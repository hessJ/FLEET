

### Author: Jonathan L. Hess (2017)
### Affiliation: SUNY Upstate Medical University, Dept. of Psychiatry & Behavioral Sciences
### Contact: hessjo@upstate.edu

### R version 3.3.0 (2016-05-03, Supposedly Educational)
### Platform: x86_64-w64-mingw32


path_to_fleet = "~/Path/to/fleet.R"

path_to_gwas = "~/Path/to/gwas.txt"


# Set working directory where fleet.R is located
setwd(paste(path_to_fleet))


## Set LD-clumping parameters for Plink

r2 = 0.1 ## minimum r2 value for a SNP for clumping
ld_window = 1000 ## minimum distance between markers for clumping
snp_field = "snpid" ## column name for SNPs
clump_field = "p" ## Column name for p-value3

## Set a p-value of enrichment (from weighted regression) to initiate permutations:

permStartThreshold = 5e-05


## Multithreaded options:
n_cores = 2 ## User specified option

## Check for multi-threaded option:
if(is.null(n_cores) == TRUE){n_cores = 1} else {n_cores = n_cores}

###

dir.create("clumped")
dir.create("pruned")
dir.create("summary_stats")
dir.create("plots")


time_stamp = format(Sys.time(), "%a-%b-%d-%Y")


# Import packages

library(plyr)
library(data.table)
library(GenomicRanges)
require(ggplot2)
require(gridExtra)
require(ggrepel)
require(foreach)
require(doParallel)
require(xtable)



# positional map of SNPs in hg19

snp137 = readRDS("g1000_CEU_bimToGranges.Rdata")
snp137 = snp137[!names(snp137) %in% "."]
snpids = as.character(names(snp137))
dup.names = snpids[duplicated(snpids)]
snp137 = snp137[!names(snp137) %in% dup.names]
snp137_df = as.data.frame(snp137)

# Pull out SNPs in the MHCx region
mhc = snp137_df[snp137_df$seqnames %in% "chr6" & snp137_df$start >= 24e6 & snp137_df$end <= 34e6, ]


## Path to 1000G files 
ref_data = list.files(path = path_to_fleet, full.names = TRUE, pattern= glob2rx("CEU*1kg_phase1*bed"))
ref_data = gsub(".bed", "", ref_data)
ref_data

### Path to GWAS file
gwas = paste(path_to_gwas, sep = "")

## Read GWAS file and extract genome-wide significant markers
read.gwas = fread(gwas, header = T)
read.gwas = as.data.frame(read.gwas)
sig.gwas = read.gwas[read.gwas[,colnames(read.gwas) %in% clump_field] < 5e-08, ]


### threshold for declaring a SNP genome-wide significant based on its -log10(P)
gwas.sig.threshold = -log10(5e-08)




##################### Module A.  Pruning 1000G CEU data for variants in high LD 
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################


### Prune 1KG reference 

for( i in 1:length(ref_data)){
  cat("\r Pruning reference data for variants in high LD (short range)",i)
  prune = paste("plink --bfile ", ref_data[[i]], " --maf .05 --indep 100 5 2 --out ",path_to_fleet,"/pruned/PRUNED_",basename(ref_data[[i]]), sep = "")
  system(prune,show.output.on.console = FALSE)
}
prune.in = list.files(path = paste(path_to_fleet,"/pruned/",sep=""), pattern = "prune.in", full.names = TRUE)
path = path[grepl("1kg_phase1_", path)]
if(length(prune.in) != 22){stop("Unexpected number of LD-pruned reference data files detected! Please check ~/pruned/ directory for errors.")}
prune.in = lapply(prune.in, function(x) read.table(x, header = F))
sum(unlist(lapply(prune.in, nrow))) ## total number of SNPs retained after high-LD pruning step


prune.in = as.data.frame(do.call(rbind, prune.in))
prune.in.plus = c(as.character(prune.in[,1]), sig.gwas[,colnames(sig.gwas) %in% clump_field])
write.table(prune.in.plus, file = paste(path_to_fleet,"/pruned/PRUNE.IN",sep = ""), quote = F, row.names = F, col.names = F)




##################### Module B.  Clumping GWAS summary statistics file using LD-pruned 1000G data
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################



## Clump summary statistics for GWAS provided in the list 'gwas' (integrates genome-wide significant markers that were dropped from LD-prune)
for( i in 1:length(gwas)){
  
  
  ## Enable multi-threading for clumping GWAS
  cl = makeCluster(n_cores)
  registerDoParallel(cl)
  
  
  print(paste("Clumping GWAS results:", basename(gwas[[i]])))
  
  
  for( j in 1:length(ref_data)){
    
    cmd = paste("plink --bfile ", ref_data[[j]]," --exclude ",path_to_fleet,"/dupvars.txt --freq --clump ", gwas[[i]], " --clump-p1 1.0 --clump-p2 1.0 --clump-r2 ", r2," --clump-kb ",ld_window," --extract ",path_to_fleet,"/pruned/PRUNE.IN --clump-snp-field ",snp_field," --clump-field ",clump_field," --out F:/ref_data/g1000/qc/clumped/CLUMP_",basename(gwas[[i]]),"_",j, sep = "")
    
    system(cmd,show.output.on.console = FALSE)
    
  }
  stopCluster(cl)
  registerDoSEQ()
  
  
}





##################### Module C.  Annotating hg19 common variants with functional/regulatory elements
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################




### Path to functional regulatory annotations
paths = list.files(path = path_to_fleet, full.names = T, pattern = "DHS|CODING|PROMOTER|SPLICESITE|UTR|FANTOM|INTRON|HISTONE|RBP|TFBS")
paths = paths[grepl(".Rdata", paths)]


## Genomic annotations by position, type, and tissue of origin
names(paths) = gsub(".Rdata", "", basename(paths))

annot.cat = list(c("CODING", "INTRON", "FIVEUTR", "THREEUTR", "PROMOTER","SPLICESITE"),
                 "DHS", c("HISTONE"), c("TFBS"), c("RBP"), c("FANTOM"))


## Overlap annotations onto SNP positions 

annot.list = list()
for( i in 1:length(paths)){
  
  cat("\r Mapping genomic annotations to SNPs",i)
  
  annot.gr = readRDS(paths[[i]])
  
  ## For splice junctions, add a 20 base window 
  if(names(paths)[[i]] == "SPLICESITE"){
    
    start(annot.gr) = start(annot.gr) - 20
    end(annot.gr) = end(annot.gr) + 20
    
  }
  
  
  ## For cell/tissue enhancers add a 500 base window 
  if(names(paths)[[i]] == "FANTOM"){
    
    start(annot.gr) = start(annot.gr) - 500
    end(annot.gr) = end(annot.gr) + 500
    
  }
  
  ## For histone/DHS, subset brain and immune cell/tissue types
  if(names(paths)[[i]] == "HISTONE" | names(paths)[[i]] == "DHS"){
    # annot.gr = annot.gr[grepl(brain.tissue, elementMetadata(annot.gr)$tissue)]
    
    cell.grab.histone = c("brain|lobe|cortex|astro|cerebellum|pyramidal|lymph|B cell|frontal|hippo|cingulate|gyrus|neuro|pericyte|neutro|monocyte|dendritic|helper|lymphocyte|spleen|thymus|killer")
    
    annot.gr = annot.gr[grepl(cell.grab.histone, annot.gr$tissue),]
    annot.gr = annot.gr[!grepl("renal|liver|kidney", annot.gr$tissue), ]
    
    elementMetadata(annot.gr)$id = paste(elementMetadata(annot.gr)$id,"_", elementMetadata(annot.gr)$tissue, sep = "" )}
  
  overlaps = findOverlaps(snp137, annot.gr) ## SNPs overlapping annotations
  
  snp.hits = names(snp137)[queryHits(overlaps)]
  
  features = elementMetadata(annot.gr)$id[subjectHits(overlaps)]
  
  
  snp.comb = data.table(SNP = as.character(snp.hits), ID =  gsub("-human", "", as.character(features)))
  
  annot.list[[i]] = snp.comb
  
}



################## 







##################### Module D.  Genome-wide enrichment test for functional annotations across LD-clumped loci  
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################


ld.df.list <- list()
glm_combine = list()
for( i in 1:length(gwas)){
  
  clumped = list.files(path =paste(path_to_fleet,"/clumped/",sep=""), full.names = T, pattern = glob2rx("CLUMP_*.clumped"))
  clumped = clumped[grepl(basename(gwas[[i]]), clumped)]
  
  
  clumped.df = lapply(clumped, function(x) fread(x, header = T, stringsAsFactors = FALSE))
  sum(unlist(lapply(clumped.df, nrow)))
  clumped.df = rbindlist(clumped.df)
  
  clumped.df$SP2[clumped.df$SP2 %in% "NONE"] <- clumped.df$SNP[clumped.df$SP2 %in% "NONE"]
  
  
  snp_list = clumped.df$SP2
  snp_list = sapply(snp_list, function(x) strsplit(x, ","))
  names(snp_list) = clumped.df$SNP
  
  snp_friends = sapply(snp_list, function(x) length(x))
  
  rep_function = function(x,y) {
    rep(x, y)
  }
  snp_reps = mapply(rep_function, x = names(snp_friends),y = snp_friends)
  ld.df = data.frame(INDEX = unlist(snp_reps), SNP= unlist(snp_list))
  split_char = strsplit(as.character(ld.df$SNP), "[(]")
  ld.df$SNP = unlist(lapply(split_char, function(x) x[[1]]))
  
  
  ## freqs in clump
  
  freqs = list.files(path = "G:/PGC_Results/clumped", full.names = T, pattern = glob2rx("CLUMP_*.frq"))
  freqs = freqs[grepl(basename(gwas[[i]]), freqs)]
  
  
  freqs.df = lapply(freqs, function(x) fread(x, header = T))
  sum(unlist(lapply(freqs.df, nrow)))
  freqs.df = rbindlist(freqs.df)
  cols = c("SNP", "MAF")
  freqs.df = freqs.df[,cols,with = FALSE]
  
  
  
  ### Calculate LD between index markers and surrounding markers:
  write.table(clumped.df$SNP, file = paste(path_to_fleet,"/clumped/INDEX.snplist",sep=""), quote = F, col.names = F)
  
  print(paste("Calculating LD for index markers"))
  for( j in 1:length(ref_data)){
    
    cmd = paste("G:/1KG/plink --bfile ", ref_data[[j]]," --ld-snp-list ", path_to_fleet,"/clumped/INDEX.snplist --ld-window-kb ", ld_window," --r2 --ld-window-r2 0.8 --out ", path_to_fleet,"/clumped/LDINT_",basename(gwas[[i]]),"_",j, sep = "")
    
    cmd = gsub("\n", "", cmd)
    cmd
    
    system(cmd, show.output.on.console = FALSE)
    
  }
  
  ld = list.files(path = "G:/PGC_Results/clumped", full.names = T, pattern = paste("LDINT_",basename(gwas[[i]]), sep = ""))
  ld = ld[!grepl(".log", ld)]
  
  ld.int = lapply(ld, function(x) fread(x, header = T))
  ld.int = rbindlist(ld.int)
  ld.int = as.data.frame(ld.int)
  ld.int = ld.int[,colnames(ld.int) %in% c("SNP_A", "CHR_A", "BP_A", "SNP_B", "R2")]
  colnames(ld.int) = c("CHR","BP","INDEX","SNP","R2")
  

  
  ld.maf = merge(ld.df , freqs.df , by = "SNP")
  
  
  maf.bin = aggregate(ld.maf$MAF, by=list(ld.maf$INDEX), function(x) sum(-log(x)))
  colnames(maf.bin) = c("INDEX", "MAF")
  
  ld.df = merge(ld.df, maf.bin, by = "INDEX")
  ld.df = ld.df[!duplicated(ld.df), ]
  
  names(clumped.df)[names(clumped.df) %in% "SNP"] = "INDEX"
  
  ## add P-value to LD interval
  clumped.df = as.data.frame(clumped.df)
  clumped.df = clumped.df[match(ld.df$INDEX, clumped.df$INDEX), ]
  ld.df$P = clumped.df$P
  
  ## Remove MHC from LD intervals 
  ld.df = ld.df[!ld.df$INDEX %in% rownames(mhc), ]
  ld.df$INDEX = as.factor(ld.df$INDEX)
  
  
  ## calculate number of LD friends
  ld_friends = unlist(table(ld.df$INDEX))
  ld_friends = as.data.frame(ld_friends)
  ld_friends = ld_friends[match(ld.df$INDEX, ld_friends$Var1), ]
  ld.df$LD_FRIENDS = ld_friends$Freq
  ld.df$LD_FRIENDS[is.na(ld.df$LD_FRIENDS)] = 0
  
  ## association scores 
  ld.df$LOGP = -log10(ld.df$P)
  ld.df$ZSCORE = qnorm(1-ld.df$P/2)
  
  ld.df.list[[i]] = ld.df
  
}





for(n in 1:length(ld.df.list)){
  
  ld.df <- ld.df.list[[n]]
  
  coef_trait = list()
  glm_trait  = list()
  
  for( a in 1:length(paths)){
    
    cat("\r Annotations ", a);flush.console()
    
    
    annot.gr = readRDS(paths[[a]])
    
    if(names(paths)[[a]] == "MF" | names(paths)[[a]] == "BP" | names(paths)[[a]] =="CC" | names(paths)[[a]] == "PANTHER"){
      
      start(annot.gr) = start(annot.gr)-20e3
      end(annot.gr) = end(annot.gr)+20e3
    }
    
    
    if(names(paths)[[a]] == "SPLICESITE"){
      
      start(annot.gr) = start(annot.gr) - 20
      end(annot.gr) = end(annot.gr) + 20
      
    }
    
    ## For cell/tissue enhancers add a 500 base window 
    if(names(paths)[[a]] == "FANTOM"){
      
      start(annot.gr) = start(annot.gr) - 200
      end(annot.gr) = end(annot.gr) + 200
      
    }
    
    ## For histone/DHS, subset brain and immune cell/tissue types
    if(names(paths)[[a]] == "HISTONE" | names(paths)[[a]] == "DHS" ){
      # annot.gr = annot.gr[grepl(brain.tissue, elementMetadata(annot.gr)$tissue)]
      
      # cell.grab.histone = c("brain|lobe|cortex|astrocy|cerebellum|pyramidal|lymph|B cell|frontal|hippo|cingulate|gyrus|neuro|pericyte|neutro|monocyte|dendritic|helper|lymphocyte|spleen|thymus|killer")
      # 
      # annot.gr = annot.gr[grepl(cell.grab.histone, annot.gr$tissue),]
      # annot.gr = annot.gr[!grepl("renal|liver|kidney", annot.gr$tissue), ]
      # 
      elementMetadata(annot.gr)$id = paste(elementMetadata(annot.gr)$id,"_", elementMetadata(annot.gr)$tissue, sep = "" )
    }
    
    overlaps = findOverlaps(snp137, annot.gr) ## SNPs overlapping annotations
    
    snp.hits = names(snp137)[queryHits(overlaps)]
    
    features = elementMetadata(annot.gr)$id[subjectHits(overlaps)]
    widths = width(annot.gr)[subjectHits(overlaps)]
    
    snp.comb = data.table(SNP = as.character(snp.hits), ID =  gsub("-human", "", as.character(features)))
    
    
    tmp = snp.comb
    tmp$widths = widths
    tmp = tmp[!duplicated(tmp)]
    
    
    # mean.widths = aggregate(tmp$widths, by=list(tmp$ID), mean)
    # mean.widths = mean.widths[order(mean.widths$x, decreasing = T), ]
    
    
    gwas.sig.threshold = -log10(5e-08)
    
    tracks = unique(tmp$ID)
    
    tmp = tmp[tmp$SNP %in% ld.df$SNP,]
    
    order.vars = as.character(unique(ld.df$INDEX))
    
    
    df.hold = ld.df[!duplicated(ld.df$INDEX), ]
    df.hold = df.hold[match(order.vars, df.hold$INDEX), ]
    
    q_plot = list()
    c_list = list()
    glm_list = list()
    keep.track.names= list()
    
    for(b in 1:length(tracks)){
      
      cat("\r Enrichment analysis in LD intervals (", names(paths)[[a]],")", b, "of", length(tracks))
      
      sub = tmp[tmp$ID %in% tracks[[b]],]
      
      if(nrow(sub) < 1) next
      
      index.hits = unique(ld.df$INDEX[ld.df$SNP %in% sub$SNP])
      sub = data.frame(INDEX = index.hits, ID = tracks[[b]], hits = 1)
      
      sub = sub[match(order.vars, sub$INDEX), ]
      sub$INDEX = order.vars
      sub$hits[is.na(sub$hits)] = 0 
      hit.counts = sub$hits
      
      
      
      df.stats = as.data.frame(do.call(cbind, list(df.hold, hit.counts)))
      
      
      new.tracks = sub("[[:punct:]]", "_", tracks[[b]])
      new.tracks = gsub("[-|,/.;:'() ]", "_", new.tracks)
      # new.tracks = paste("ID_", new.tracks, sep="")
      colnames(df.stats)[ncol(df.stats)] = "hits"
      
      keep.track.names[[b]] = new.tracks
      
      freq_annot = table(df.stats[,ncol(df.stats)])/nrow(df.stats)
      
      if(length(freq_annot) < 2) next
      
      ## Enrichment analysis
      
      
      gwas.sig.annot.hits = df.stats$INDEX[df.stats$hits  == 1 & df.stats$LOGP >= gwas.sig.threshold]
      n.gwas.hits = length(gwas.sig.annot.hits)
      
      
      gw.vars = paste(gwas.sig.annot.hits, collapse = "|", sep = "")
      
      ### threshold SNPs by GWAS significance 
      
      group_counts = table(df.stats[,ncol(df.stats)])
      min_group = min(group_counts)
      # formula  = formula(paste(tracks[[b]], " ~ ZSCORE + MAF + LD_FRIENDS"))
      # fit = speedglm(formula, data = df,  family=binomial('logit'))
      # df = df.stats[!is.infinite(df.stats$ZSCORE), ]
      
      
      
      lmform = formula(paste("ZSCORE ~ LD_FRIENDS  + MAF + hits"))
      fit = lm(lmform, data = df.stats, weights = (1/df.stats$MAF))
      # fit = logistf::logistf(binary_hits  ~ PT + MAF, data = df)
      coef = summary(fit)$coefficients
      coef = as.data.frame(coef)
      # coef$OR = exp(coef$Estimate)
      coef$pred = rownames(coef)
      coef$df = paste(summary(fit)$df, sep ="", collapse= "//")
      coef$resp = tracks[[b]]
      coef$rsq = summary(fit)$r.squared
      coef$GWS_VAR_HITS = gw.vars
      coef$freq = freq_annot[[2]]
      
      ## Rapid permutation (GW-SIG):
      pval.check = coef$`Pr(>|t|)`[coef$pred %in% "hits"] 
      
      if(pval.check < permStartThreshold) {
        
        gws.hits.1 = df.stats[df.stats$P < 5e-08, ]
        gws.hits.2 = df.stats[df.stats$P < 1e-06, ]
        gws.hits.3 = df.stats[df.stats$P < 1e-04, ]
        
        gw.sig.list = list(gws.hits.1, gws.hits.2, gws.hits.3)
        names(gw.sig.list) = c(5e-08, 1e-06, 1e-04)
        
        perm_p = list(slot1 = NA, slot2= NA, slot3=NA)
        
        for( s in 1:length(gw.sig.list)){ 
          
          set = gw.sig.list[[s]]
          
          set.hits = sum(set$hits)
          
          cat("\r Permutations for GW-significance level:", names(gw.sig.list)[[s]])
          
          cat("\n")
          
          if( nrow(set) > 0){
            
            
            set.dist = set$MAF
            
            gwas.maf = mean(set.dist)
            gwas.maf.sd = 1.96 * (sd(set.dist)/sqrt(length(set.dist)))
            
            maf.bounds = c(gwas.maf - gwas.maf.sd, gwas.maf + gwas.maf.sd)
            gwas.sig.annot.hits = paste(gwas.sig.annot.hits, collapse = "|")
            
            
            perm_hits = list()
            
            for(p in 1:1000){
              
              cat("\r MAF-matched permutation",p)
              repeat{
                
                grab = df.stats[df.stats$MAF > maf.bounds[[1]] & df.stats$MAF < maf.bounds[[2]], ]
                compare = t.test(grab$MAF, set.dist)
                
                if(compare$p.value > .05) break}
              
              grab = grab[sample(nrow(grab)), ]
              grab = grab$hits[1:nrow(set)]
              
              perm_hits[[p]] = sum(grab > 0)/length(grab)
              
            }
            
            perm_p[[s]] = sum((set.hits/nrow(set)) < unlist(perm_hits)) / 1000
            
          }
        }
      }
      
      perm_df <- data.frame(GW1_5e08_perm = perm_p[[1]], 
                            GW2_1e06_perm = perm_p[[2]],
                            GW3_1e04_perm = perm_p[[3]])
      
      coef = cbind(coef, perm_df)
      
      
      
      glm_list[[b]] = coef
      
      
      
    }
    
    
    
    cat("\n")
    
    
    glm_all = ldply(glm_list)
    glm_all = glm_all[glm_all$pred %in% "hits",]
    # glm_all$pred = unlist(keep.track.names)[[1]]
    glm_all$fdr = p.adjust(glm_all[,4], "fdr")
    glm_all$bonf = p.adjust(glm_all[,4], "bonferroni")
    glm_all = glm_all[order(glm_all[,4], decreasing= F), ]
    glm_all$annot = names(paths)[[a]]
    glm_all$trait = basename(gwas[[n]])
    
    
    glm_all = glm_all[order(glm_all$`Pr(>|t|)`,decreasing = T),]
    
    glm_trait[[a]] = glm_all
    
    
    # write.csv(coef_all, file = paste("C:/Users/HessJo/Dropbox/encode/summary_stats/", names(paths)[[a]], "_", paste("CLUMP_",basename(gwas[[i]]), sep = ""), "_REGRESSION.csv", sep = "" ), quote = F, row.names = F)
    write.csv(glm_all, file = paste("C:/Users/HessJo/Dropbox/encode/summary_stats/LM-PERM_", names(paths)[[a]], "_", paste("CLUMP_",basename(gwas[[n]]), sep = ""), "_REGRESSION.csv", sep = "" ), quote = F, row.names = F)
    
    gc()
  }
  
  glm_combine[[n]] = glm_trait
  
}



##################### Module E.  Plotting results and exporting summary statistics to an html file
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################




  
  names(annot.cat) = c("Gene regions", "DHS", "Chromatin modifications",
                       "TF binding sites", "RBP binding sites",
                       "Enhancers")
  
  
  
  for (g in 1:length(gwas)){
    
    basename = paste("*",basename(gwas[[g]]),"*_REGRESSION.csv", sep = "")
    results = list.files(path = paste(path_to_fleet, "/summary_stats/", sep = ""), pattern = glob2rx(basename), full.names = TRUE)
    results = lapply(results, function(x) read.csv(x, header = T))
    results = plyr::ldply(results)
    results = results[!is.na(results$annot), ]
    
    # results = results[results$freq > .01, ]
    # results = results[results$annot %in% unlist(annot.cat), ]
    
    
    # results = lapply(glm_combine, function(x) x[[1]])
    # results = ldply(results)
    #
    results$P = results$Pr...t..
    results$enrich = results$t.value
    
    results$logP = -log10(results$P)
    results$bonf = p.adjust(results$P, "bonferroni")
    results$fdr = p.adjust(results$P, "fdr")
    results = results[order(results$annot), ]
    annot.type = unique(results$annot)
    results$label.top = NA
    results$resp = gsub("_differentially_expressed_enhancers", "", results$resp)
    results$resp = toupper(results$resp)
    
    split = split(results, as.character(results$annot))
    split = lapply(split, function(x) x[order(x$P, decreasing = F), ])
    
    for(m in 1:length(split)){
      
      if(names(split)[[m]] == "HISTONE" | names(split)[[m]] == "DHS"){
        
        cell.grab.histone = c("brain|lobe|cortex|astro|cerebellum|pyramidal|lymph|B cell|frontal|hippo|cingulate|gyrus|neuro|pericyte|neutro|monocyte|dendritic|helper|lymphocyte|spleen|thymus|killer")
        split[[m]] = split[[m]][grepl(cell.grab.histone,ignore.case = T, split[[m]]$resp),]
        split[[m]] = split[[m]][!grepl("muscle|renal|liver|medalis|gastro",ignore.case = T, split[[m]]$resp),]
        
      }
      
      if(nrow(split[[m]]) < 2){ split[[m]]$label.top[[1]] = as.character(split[[m]]$resp[[1]])}
      if(nrow(split[[m]]) > 2){ split[[m]]$label.top[1:2] = as.character(split[[m]]$resp[1:2])}
    }
    
    results = as.data.frame(do.call(rbind, split))
    # results$label.top[results$bonf > .05] <- NA
    
    col = rainbow(length(annot.cat))
    names(col) = names(annot.cat)
    col = as.data.frame(col)
    col$features = as.character(rownames(col))
    
    n_tests = nrow(results)
    alpha = .05/n_tests
    alpha.z = qnorm(alpha)
    
    ## group annot.cat
    
    groups =  lapply(annot.cat, data.frame)
    groups = ldply(groups)
    colnames(groups) = c("features","annot")
    
    
    results = merge(results, groups, by = "annot")
    results = merge(results, col , by ="features")
    
    results = results[order(results$features), ]
    
    results$labels = results$resp
    
    results$label.top[results$P > .05 & results$enrich < 0] <- NA
    
    
    
    # results$z.value[results$z.value > 8] <- 8
    # results$z.value[results$z.value < -8] <- -8
    
    # results = results[results$enrich > 0, ]
    
    trait.names = unique(results$trait)
    
    g = ggplot(results, aes(x = results$features, y= results$enrich, col = results$features)) + geom_point(position="jitter") + scale_colour_brewer(palette="Dark2") + theme_bw() +
      xlab(NULL) + guides(col = FALSE) + geom_hline(yintercept = c(-alpha.z, 0, alpha.z), col = c("red", "black", "red"), linetype = c("dotted", "solid", "dotted"), size = 0.5)  +
      ylim(min = min(-10, min(results$enrich*1.5)), max = max(10, max(results$enrich*1.5))) +
      ylab("Annotation Enrichment (z-score)")  + geom_label_repel(label = results$label.top, size = 2, segment.color = NA)
    
    
    
    
    
    png(paste(path_to_fleet,"/broad_plots/JITTERPLOT_TRAIT=",trait.names,".png", sep = ""), res = 300, units = "in",  width= 8, height = 7)
    print(g)
    dev.off()
    
    
    sig = results
    
    box = ggplot(sig, aes(x = sig$features, y = sig$logP, fill = sig$features))+ theme_bw() + geom_hline(yintercept = -log10(alpha), col = "red", linetype="dashed") +
      geom_boxplot(width = 0.5) + xlab(NULL) + ylab("Significance (OR > 0), -log10(P)") +scale_fill_brewer(palette="Dark2") +guides(fill = FALSE) 
    
    png(paste(path_to_fleet,"/broad_plots/BOXPLOT_TRAIT=",trait.names,".png", sep = ""), res =300, units = "in", height = 5, width = 7)
    print(box)
    dev.off()
    
    ### X table 
    
    
    require(xtable)
    
    sig = results[results$enrich > 0 & results$bonf < .05, ]
    sig$PVAL = format(sig$P, scientific=  T, digits = 3)
    sig$BONF= format(sig$bonf, scientific = T, digits = 3)
    sig$FREQ = format(sig$freq*100, digits= 1)
    sig$BETA = sig$Estimate
    sig$SE = sig$`Std. Error.`
    xtab = xtable(sig[,colnames(sig) %in% c("features", "resp", "BETA", "SE", "BONF", "FREQ", "PVAL")])
    print.xtable(xtab, type="html", file= paste(path_to_fleet,"/enrichmentreport_",trait.names,".html", sep = ""))
    
    
  }
  
  
  
  
  
