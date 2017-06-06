

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
paths = list.files(path = path_to_fleet, full.names = T, pattern = "CC|MF|PANTHER|DHS|CODING|PROMOTER|SPLICESITE|UTR|FANTOM|INTRON|HISTONE|RBP|TFBS")
paths = paths[grepl(".Rdata", paths)]


## Genomic annotations by position, type, and tissue of origin
names(paths) = gsub(".Rdata", "", basename(paths))

annot.cat = list(c("CODING", "INTRON", "FIVEUTR", "THREEUTR", "PROMOTER","SPLICESITE"), c("PANTHER", "CC","MF"),
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
  
  clumped = list.files(path = "F:/ref_data/g1000/qc/clumped/", full.names = T, pattern = glob2rx("CLUMP_*.clumped"))
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
  
  freqs = list.files(path = "F:/ref_data/g1000/qc/clumped", full.names = T, pattern = glob2rx("CLUMP_*.frq"))
  freqs = freqs[grepl(basename(gwas[[i]]), freqs)]
  
  
  freqs.df = lapply(freqs, function(x) fread(x, header = T))
  sum(unlist(lapply(freqs.df, nrow)))
  freqs.df = rbindlist(freqs.df)
  cols = c("SNP", "MAF")
  freqs.df = freqs.df[,cols,with = FALSE]
  
  
  
  ### Calculate LD between index markers and surrounding markers:
  write.table(clumped.df$SNP, file = "G:/PGC_Results/clumped/INDEX.snplist", quote = F, col.names = F)
  
  print(paste("Calculating LD for index markers"))
  for( j in 1:length(ref_data)){
    
    # cmd = paste("G:/1KG/plink --bfile ", ref_data[[j]]," --exclude F:/ref_data/g1000/dupvars.txt --ld-snp-list G:/PGC_Results/clumped/INDEX.snplist --ld-window-kb ", ld_window," --r2 --ld-window-r2 0.6 --out G:/PGC_results/clumped/LDINT_",basename(gwas[[i]]),"_",j, sep = "")
    cmd = paste("G:/1KG/plink --bfile ", ref_data[[j]]," --ld-window 99999 --exclude F:/ref_data/g1000/dupvars.txt --ld-snp-list G:/PGC_Results/clumped/INDEX.snplist --ld-window-kb ", ld_window," --r2 --ld-window-r2 0.1 --out G:/PGC_results/clumped/LDINT_",basename(gwas[[i]]),"_",j, sep = "")
    
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
  
  ### subset ld.df by SNPs in high ld (r>0.6) to index marker
  # ld.df = ld.df[ld.df$SNP %in% ld.int$SNP, ]
  
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
  
  ## calculate LD scores for index SNPs 
  
  ld.score = aggregate(ld.int$R2, by = list(ld.int$INDEX), sum)
  colnames(ld.score) = c("INDEX", "LDSCORE")
  
  ## association scores 
  ld.df$LOGP = -log10(ld.df$P)
  ld.df$P[ld.df$P < 1e-15] <- 1e-15
  ld.df$ZSCORE = qnorm(1-ld.df$P/2)
  ld.df = merge(ld.df, ld.score, by="INDEX")
  
  ld.df.list[[i]] = ld.df
  
}

## FLEET enrichment analysis


start_clock = proc.time()

for(n in 1:length(ld.df.list)){
  
  ld.df <- ld.df.list[[n]]
  
  coef_trait = list()
  glm_trait  = list()
  
  for( a in 1:length(paths)){
    
    cat("\r    Loading annotations:",a,"...")

    annot.gr = readRDS(paths[[a]])
    
    
    # merge data table with ld.df stats
    tmp = merge(data.table(ld.df[,c("SNP","INDEX")]), annot.gr, by = "SNP", all.x = T)
    
    tmp = as.data.frame(tmp)
    tmp = tmp[,!colnames(tmp) %in% "SNP"]
    tmp = tmp[!duplicated(tmp$INDEX), ]
    tmp[is.na(tmp)] = 0 
    
    # merge index information
    all_df = merge(data.table(ld.df[,!colnames(ld.df) %in% "SNP"]), tmp, by = "INDEX")
    all_df = all_df[!duplicated(all_df)]
    all_df = as.data.frame(all_df)
    
    
    new.tracks = sub("[[:punct:]]", "", colnames(all_df))
    new.tracks = gsub("[-|,/.;:'() ]", "", colnames(all_df))
    new.tracks = gsub("[+]", "_", new.tracks)
    colnames(all_df) = new.tracks
    
    # threshold for declaring genome-wide significace 
    gwas.sig.threshold = -log10(5e-08) # threshold for calling INDEX marker genome-wide significant
    
    # grab unique track names (for looping)
    tracks = colnames(all_df)[-c(1:7)]
    names(all_df)[names(all_df) %in% tracks] = paste("ID_", tracks, sep = "")
    tracks = paste("ID_", tracks, sep = "")
    
    q_plot = list()
    c_list = list()
    glm_list = list()
    keep.track.names= list()
    bounds_df = list()
    
    for(b in 1:length(tracks)){
      
      cat("\r Enrichment analysis in LD intervals (", names(paths)[[a]],")", b, "of", length(tracks))
      # 
      freq_annot = table(all_df[,colnames(all_df) %in% tracks[[b]]])/nrow(all_df)
      
      if(length(freq_annot) < 2) next
      
      ## Enrichment analysis
      
      
      gwas.sig.annot.hits = all_df[all_df$LOGP >= gwas.sig.threshold & all_df[,colnames(all_df) %in% tracks[[b]]] == 1, ]
      n.gwas.hits = nrow(gwas.sig.annot.hits)
      
      
      gw.vars = paste(gwas.sig.annot.hits$INDEX, collapse = "|", sep = "")
      
      
      lmform = formula(paste("ZSCORE ~ LDSCORE  + MAF  + ", tracks[[b]]))
  
      
      fit = lm(lmform, data = all_df, weights = (1/all_df$MAF))
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
      coef$annot = basename(lapply(strsplit(paths[[a]], "[+]"), function(x) x[[1]])[[1]])
      
      ## Rapid permutation (GW-SIG):
      pval.check = coef$`Pr(>|t|)`[coef$pred %in% tracks[[b]]] 
      pval.check

      perm_df <- data.frame(GW1_5e08_perm = NA,
                            GW2_1e06_perm = NA,
                            GW3_1e04_perm = NA)

      if(pval.check < permStartThreshold & coef$Estimate[coef$pred %in% tracks[[b]]] > 0) {

        gws.hits.1 = all_df[all_df$P < 5e-08, ]
        gws.hits.2 = all_df[all_df$P < 1e-06, ]
        gws.hits.3 = all_df[all_df$P < 1e-04, ]

        gw.sig.list = list(gws.hits.1, gws.hits.2,gws.hits.3)
        names(gw.sig.list) = c(5e-08, 1e-06, 1e-04)

        perm_p = list(slot1 = NA, slot2= NA, slot3=NA)
        
        bounds = list()
        for( s in 1:length(gw.sig.list)){

          set = gw.sig.list[[s]]

          if(length(tracks) > 1){set.hits = colSums(set[,colnames(set) %in% tracks])} else {set.hits = sum(set[,colnames(set) %in% tracks])}
         
          # iterate over annotations, calculate summary statistics
          premelt = set[,colnames(set) %in% c("INDEX", tracks)]
          premelt = melt(premelt)
          hit_melt = premelt[premelt$value == 1, ]
          
          snp_var_names = paste(premelt$INDEX, collapse = "|", sep = "")
          
          if(nrow(set) < 1) next
          
          premelt.counts = table(hit_melt$variable)
          premelt = merge(set[,!colnames(set) %in% tracks], premelt, by ="INDEX")

          if( nrow(premelt) > 0){

            gwas.maf = aggregate(premelt$MAF, by=list(premelt$variable), mean)
            gwas.maf.sd = aggregate(premelt$MAF, by=list(premelt$variable), function(x) c(sd(x), sqrt(length(x))))
            gwas.maf.sd$ci = 1.96 * gwas.maf.sd$x[,1]/gwas.maf.sd$x[,2]
            
            gwas.ld = aggregate(premelt$LDSCORE, by=list(premelt$variable), mean)
            gwas.ld.sd = aggregate(premelt$LDSCORE, by=list(premelt$variable), function(x) c(sd(x), sqrt(length(x))))
            gwas.ld.sd$ci = 1.96 * gwas.ld.sd$x[,1]/gwas.ld.sd$x[,2]
            
            gwas.counts = nrow(set)
            gwas.hits = table(as.character(premelt$variable))
            
            
            t.out =  data.frame(id = gwas.maf$Group.1,
                                     pthres = names(gw.sig.list)[[s]], 
                                     n_gwas_var = gwas.counts,
                                     n_gwas_var_hits = as.vector(gwas.hits), 
                                     mean_maf = gwas.maf$x, 
                                     ci_maf = gwas.maf.sd$ci,
                                     mean_ldscore = gwas.ld$x, 
                                     ci_ldscore = gwas.ld.sd$ci,
                                    index = snp_var_names)
            
            t.out$ci_maf[is.na(t.out$ci_maf)] = t.out$mean_maf[is.na(t.out$ci_maf)]
            t.out$ci_ldscore[is.na(t.out$ci_ldscore)] = t.out$mean_ldscore[is.na(t.out$ci_ldscore)]

            bounds[[s]] = t.out
            

          }
          
        }
        
        bounds_df[[b]] = do.call(rbind, bounds)
        

        # perm_df <- data.frame(GW1_5e08_perm = perm_p[[1]],
                              # GW2_1e06_perm = perm_p[[2]])
      }
      
      # coef = cbind(coef, bounds_df)
      
      
      glm_list[[b]] = coef
      
      
      
    }
    
    all_bounds = ldply(bounds_df)
    all_bounds = all_bounds[!duplicated(all_bounds), ]
    
    split_bounds = split(all_bounds, all_bounds$id)
    
    no_rows = lapply(split_bounds, nrow)
    no_rows = which(unlist(no_rows) < 1)
    
    if(length(no_rows) > 0){split_bounds = split_bounds[-no_rows]}
    
    all_df = data.table(all_df)
    
    ## Permutations over set 
    
    if(length(split_bounds) > 0){

    for(p in 1:length(split_bounds)){
      
      cat("\rPermutation of set:", p)
      splim = split_bounds[[p]]
      
      if(nrow(splim) < 1) next
      
      
      grab = lapply(1:nrow(splim), function(x) all_df[all_df$MAF > splim$mean_maf[[x]] - splim$ci_maf[[x]] &  
                                                        all_df$MAF < splim$mean_maf[[x]] + splim$ci_maf[[x]]])
      
    
      sizes = lapply(grab, nrow)
      
      
      
      grab = lapply(grab, data.frame)
      size_sample = splim$n_gwas_var
     
      size_error = which(unlist(sizes) < unlist(size_sample))
  
      
      sample_sum = list()
      for( perms in 1:nPerms){
       
      # cat("\r Sampling SNPs from LD-score matched distribution:",perms)
      
      sampled = lapply(1:length(size_sample), function(x) grab[[x]][sample(nrow(grab[[x]]), size_sample[[x]]), ])
      sampled = lapply(sampled, function(x) x[,colnames(x) %in% names(split_bounds)[[p]]])
      sampled = lapply(sampled, sum)
      names(sampled) = splim$pthres
      sample_sum[[perms]] = sampled
      
      }
      
      if(nrow(splim) == 1){sdf = unlist(sample_sum); splim$perm_pval = sum(sdf > splim$n_gwas_var_hits)/nPerms}

      
      if(nrow(splim ) > 1){
        
        sdf = as.data.frame(do.call(rbind, sample_sum))
        sdf_sum = lapply(sdf, function(x) unlist(x))
        
        sdf_pval = list()
        for( z in 1:nrow(splim)){
          sdf_pval[[z]] = sum(sdf_sum[[z]]/size_sample > splim$n_gwas_var_hits[[z]]/splim$n_gwas_var[[z]])/nPerms
        }
        
        splim$perm_pval = unlist(sdf_pval)
      
        
      }
      
      
      split_bounds[[p]] = splim
      
    }

    perm_annot = as.data.frame(do.call(rbind, split_bounds))
    perm_annot$perm_pval[perm_annot$perm_pval == 0] = 1/nPerms
    
    write.csv(perm_annot, file = paste(summary_stats, names(paths)[[a]], "_", paste("CLUMP_",basename(gwas[[n]]), sep = ""), "_PERMUTATION.csv", sep = "" ), quote = F, row.names = F)
    }
    
    cat("\n")
    
    
    glm_all = ldply(glm_list)
    glm_all = glm_all[glm_all$pred %in% tracks,]
    glm_all = glm_all[order(glm_all[,4], decreasing= F), ]
    glm_all$trait = basename(gwas[[n]])
    glm_all = glm_all[order(glm_all$`Pr(>|t|)`,decreasing = T),]
    
    glm_trait[[a]] = glm_all
    
    
    write.csv(glm_all, file = paste(summary_stats, names(paths)[[a]], "_", paste("CLUMP_",basename(gwas[[n]]), sep = ""), "_REGRESSION.csv", sep = "" ), quote = F, row.names = F)
   
    gc()
  }
  
  glm_combine[[n]] = glm_trait
  
}




##################### Module E.  Plotting results and exporting summary statistics to an html file
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################




  
  names(annot.cat) = c("Gene regions",
                       "Pathways",
                       "DHS", "Chromatin modifications",
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
    sig$FDR= format(sig$fdr, scientific = T, digits = 3)
    sig$FREQ = format(sig$freq*100, digits= 1)
    sig$BETA = sig$Estimate
    sig$SE = sig$`Std. Error.`
    xtab = xtable(sig[,colnames(sig) %in% c("features", "FDR", "resp", "BETA", "SE", "BONF", "FREQ", "PVAL")])
    print.xtable(xtab, type="html", file= paste(path_to_fleet,"/enrichmentreport_",trait.names,".html", sep = ""))
    
    
  }
  
  
  
  
  
