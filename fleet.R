#!/usr/bin/env Rscript


# Import packages
suppressMessages(library(optparse))
suppressMessages(require(argparse))

## Auto-detect the directory name for fleet.R
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "fleet.R")

# create log file
srcPath = dirname(script.name)
#log_con <- file(paste(srcPath,"fleet.log",sep=""),open="a")

# set masthead

masthead = as.character("
==========================================================
*
* Functional LD-interval EnrichmEnt Test (FLEET)
*
* Jonathan L. Hess, PhD and Stephen J. Glatt, PhD (c) 2017
*
* SUNY Upstate Medical University, PsychGENe Lab
*
* Contact: hessjo@upstate.edu
*
* https://github.com/hessJ/FLEET
*
* GNU GENERAL PUBLIC LICENSE v3
===========================================================\n")

#cat(masthead, file = log_con)
cat(masthead) # display masthead 


# ========================================================================
started = Sys.time()
cat("\n")
cat(paste("Start time:", started))
cat("\n")

cat("\nLocation of fleet:",script.name,"\n")

# create parser object
parser <- ArgumentParser()

option_list = list(
  make_option(c("-G", "--gwas"), type="character", default="", 
              help="Path to GWAS summary statistics. Column headers are required. Allowed delim = sep, tab, or comma", metavar="character"),
  
  make_option(c("-O", "--out"), type="character", default="fleetOut", 
              help="output file name [default = %default]", metavar="character"),
  
  make_option(c("-R", "--r2"), type="double", default=0.6, 
              help="R-squared threshold for linkage disequilibrium calculations [default = %default]", metavar="double"),
  
  make_option(c("-W", "--ld-window"), type="integer", default=1000, 
              help="Size of window (kilobases) for calculating linkage disequilibrium [default = %default]", metavar="integer"),
  
  make_option(c("-S", "--snp-field"), type="character", default="", 
              help="SNP column header in GWAS file", metavar="character"),
  
  make_option(c("-P", "--pcol"), type="character", default="", 
              help="P-value column header in GWAS file", metavar="character"),

  make_option(c("--robust"), type="logical", default=TRUE,
              help="Computing robust standard errors using White method (via vcovHC function in sandwich pkg) [default = %default]", metavar="double"),
  
  make_option(c("-N", "--nPerms"), type="integer", default=1000, 
              help="Number of permutations to perform [default = %default]", metavar="integer"),
  
  make_option(c("-T", "--threads"), type="integer", default=1, 
              help="Number of cores for parallel operations [default = %default]", metavar="integer"),
  
  make_option(c("-L", "--label-annotations"), type="character", default="", 
              help="Path to file containining labels for annotation sources", metavar="character"),
  
  make_option(c("-D", "--rd-annots"), type="character", default=paste(srcPath,"/annotations/",sep=""),
              help="Path to .Rdata annotations [default = %default]", metavar="character"),
  
  make_option(c("-F", "--annot-cnt"), type="double", default=10, 
              help="Minimum annotation count observed across LD-clumps [default = %default]", metavar="double"),
  
  ## Run options....
  
  make_option(c("-M", "--fleet-prune-ref"), type="logical", default = TRUE, 
              help="Initiate LD-pruning step of 1KG reference data. Only needs to be run once. [default = %default]", metavar="logical"),
  
  make_option(c("-A", "--fleet-annotate"), type="logical", default=TRUE, 
              help="Annotate LD-clumps with bedtools [default = %default]", metavar="logical"),
  
  make_option(c("-E", "--fleet-enrichment"), type="logical", default=TRUE, 
              help="Perform enrichment analysis with weighted linear models [default = %default]", metavar="logical"),
  
  make_option("--fleet-permutation", type="logical", default=TRUE, 
              help="Perform enrichment analysis with permutation (randomizing annotations) [default = %default]", metavar="logical"),
  
  make_option("--fast-permutation", type="logical", default=FALSE, 
              help="Simple permutation analysis [default = %default]", metavar="logical"),
  
  make_option("--robust-permutation", type="logical", default=FALSE, 
              help="Permutation analysis that will sample variants from the MAF bin of target SNPs [default = %default]", metavar="logical"),
  
  make_option("--speed", type="character", default='fast', 
              help="Change behavior of linear models (fast mode: SET becomes response variable, slow mode: Z-score becomes response variable) [default = %default]", metavar="character"),
  
  make_option("--pthres", type="character", default="",
              help="Table with P-value threshold(s) for SNP bins [default P-values < 5e-08, 1e-07, and 1e-06]", metavar="character"),
  
  make_option("--plots", type="logical", default=FALSE,
              help="Turning this on will produce multiple plots to display summary statistics [default = %default]", metavar="logical")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# handling empty arguments
if (opt$gwas == ""){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

if (opt$`snp-field` == ""){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

if (opt$pcol == ""){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}
# ========================================================================

flags_in_use = as.data.frame(t(as.data.frame(opt)))
rownames(flags_in_use) = paste("--", rownames(flags_in_use), sep = "")
flags_in_use$id = rownames(flags_in_use)
flags_in_use = flags_in_use[!flags_in_use$V1  == FALSE & !flags_in_use$V1 %in% "", ]

flags_in_use$stdout = paste("    ", gsub("[.]", "-", rownames(flags_in_use)), flags_in_use$V1, "\n", sep = " ")

cat("\nFlags in use:\n", flags_in_use$stdout)


# //// Dependencies
cat("\nPackages loading...\n")
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doParallel)))
suppressWarnings(suppressMessages(library(lmtest)))
suppressWarnings(suppressMessages(library(sandwich)))
suppressWarnings(suppressMessages(library(doParallel)))
suppressWarnings(suppressMessages(library(broom)))
cat("Successfully loaded packages!\n")



# //// Functions 

# ========================================================================
# Import SNP coordinates from 1KG reference:

readSNPcoords  = function(){
  cat("Reading SNP coordinates from 1KG bim files")
  read_in = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
  read_in = lapply(read_in, function(x) fread(x, header = FALSE, showProgress=FALSE))
  read_df = rbindlist(read_in)
  dups = as.data.frame(read_df[duplicated(read_df$V2), ])
  dups = unique(dups$V2)
  write.table(dups, file = paste(srcPath,'/qc/ref/duplicatevars.txt',sep=""), quote = F, row.names = F, col.names = F)
  colnames(read_df) = c("seqnames", "SNP", "empty", "start", "A1", "A2")
  read_df = as.data.frame(read_df)
  read_df$seqnames = paste("chr", read_df$seqnames, sep = "")
  read_df = read_df[!read_df$SNP %in% dups, ]
  snp137_df <<- read_df
}


# pull out snps in the xMHC region
xMHCregion = function(){
  mhc =  snp137_df[snp137_df$seqnames %in% "chr6" & snp137_df$start >= 24e6 & snp137_df$start <= 35e6, ]
  snp137_df_nomhc <<- snp137_df[!snp137_df$SNP %in% mhc$SNP, ]
  mhc <<- mhc
}

# Set annotations for FLEET
readAnnots = function(){
  read_annots = read.table(opt$`label-annotations`, header = T, sep = "\t")
  annot.cat <<- read_annots
  
  # path to annotation .Rdata files
  rdata_files = list.files(path = opt$`rd-annots`, full.names = T, pattern = ".Rdata")
  
  if(nrow(read_annots) != length(rdata_files)){stop("Number of .Rdata files does not match table with annotation labels!")}
}



# read GWAS file:
fleetReadGWAS = function(){
  
  # simple delim detection
  readlines = readLines(opt$gwas, n = 1)
  tab = strsplit(readlines, "\t")[[1]]
  space = strsplit(readlines, "[ ]")[[1]]
  comma = strsplit(readlines, "[,]")[[1]]
  delim_size = lapply(list(tab, space, comma), length)
  
  # infer delim:
  infer = which(unlist(delim_size) == max(unlist(delim_size)))
  if(infer == 1){cat("\nGuessing that GWAS file is tab-delimited\n"); read.gwas = fread(opt$gwas, header = T, sep = "\t", showProgress=FALSE)}
  if(infer == 2){cat("\nGuessing that GWAS file is space-delimited\n"); read.gwas = fread(opt$gwas, header = T, sep = " ", showProgress=FALSE)}
  if(infer == 3){cat("\nGuessing that GWAS file is comma-delimited\n"); read.gwas = fread(opt$gwas, header = T, sep = ",", showProgress=FALSE)}
  read.gwas = as.data.frame(read.gwas)
  
  cat("Writing SNP list...\n")
  snp.list = read.gwas[,colnames(read.gwas) %in% opt$`snp-field`]
  # snp.list = snp.list[snp.list %in% snp137_df$SNP]
  snp.list = data.table(SNP = snp.list)
  fwrite(snp.list, nThread = opt$threads, file = paste(srcPath, "/qc/pruned/gwas.snp.list", sep = ""), quote = F, row.names = F, col.names = F)
  
  sumstats <<- read.gwas
  
}


### Prune 1KG reference 
refPrune = function(){
  
  ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
  ref_data = gsub(".bim", "", ref_data)
  
  for( i in 1:length(ref_data)){
    indep.cmd = paste('--indep-pairwise 100 1 0.1')
    cat("\nPruning variants to approximate linkage equilibrium ( Plink command:",indep.cmd,"): ",i)
    prune = paste("plink --bfile ", ref_data[[i]], " --maf .05 --geno .05 --mind .05 --hwe 1e-06 --freq --extract ", srcPath,"/qc/pruned/gwas.snp.list ",indep.cmd," --out ", srcPath,"/qc/pruned/PRUNED_",basename(ref_data[[i]]), sep = "")
    
    system(prune, ignore.stdout = TRUE)
    
  }
  
}


## Prune GWAS file
fleetPrune = function() { 
  
  prune.in = list.files(path = paste(srcPath,"/qc/pruned/",sep=""), pattern = "prune.in", full.names = TRUE)
  prune.in = lapply(prune.in, function(x) read.table(x, header = F, stringsAsFactors = FALSE))
  total_indep = sum(unlist(lapply(prune.in, nrow)))
  cat("\nNumber of variants retained in reference panel after LD-pruning:", total_indep, "\n")
  
  prune.in = as.data.frame(do.call(rbind, prune.in))
  prune.in = prune.in[!duplicated(prune.in$V1), ]
  
  read.gwas = sumstats
  sig.gwas = read.gwas[read.gwas[,colnames(read.gwas) %in% opt$pcol] < 5e-08, ]
  
  prune.in.plus = c(as.character(prune.in), sig.gwas[,colnames(sig.gwas) %in% opt$`snp-field`])
  prune.in.plus = prune.in.plus[!duplicated(prune.in.plus)]
  
  write.table(prune.in.plus, file = paste(srcPath,"/qc/pruned/",opt$out, "_PRUNE.IN",sep=""), quote = F, row.names = F, col.names = F)
  
  names(read.gwas)[names(read.gwas) %in% c(opt$`snp-field`, opt$pcol)] = c("SNP", "P")
  sumstats <<- read.gwas
  
}


## Clump GWAS data 

# fleetClump = function() {
#   
#   
#   ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
#   ref_data = gsub(".bim", "", ref_data)
#   
#   
#   ## Enable multi-threading for clumping GWAS
#   cat("\nUsing", opt$threads, "threads for parallel processes\n")
#   
#   # cl = makeCluster(opt$threads)
#   # registerDoParallel(cl)
#   
#   cat("\nEnrichment analysis of trait:", basename(opt$gwas), "\n")
#   
#   
#   cat("\nClumping GWAS results\n")
#   
#   # foreach( j = 1:length(ref_data), .packages=c("optparse", "argparse")) %dopar% {
#   
#   cmd = list()
#   for(j in 1:length(ref_data)){
#     
#     # ensure that spaces are escaped in GWAS file path
#     gwas_file = gsub(" ", "\\ ", opt$gwas, fixed = T)
#     
#     cmd[[j]] = paste("plink --bfile ", ref_data[[j]]," --exclude ", srcPath,"/qc/ref/duplicatevars.txt --freq --clump ", gwas_file, " --clump-p1 1.0 --clump-p2 1.0 --clump-r2 ", opt$r2," --clump-kb ",opt$`ld-window`," --extract ",srcPath,"/qc/pruned/",opt$out,"_PRUNE.IN --clump-snp-field ", opt$`snp-field`," --clump-field ", opt$`clump-field`," --out ",srcPath,"/qc/clumped/CLUMP_",opt$out,"_",j, sep = "")
#     
#     
#   }
#   
#   mclapply(cmd, function(x) system(x, ignore.stdout = TRUE), mc.cores  = opt$threads)
#   gc()
#   
#   # stopCluster(cl)
#   # registerDoSEQ()
#   
#   cat("\nCompleted clumping!\n")
#   
# }


## Combine LD-clumps with summary statistics
fleetLD = function() {
  
  
  # read 1KG reference data
  ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
  ref_data = gsub(".bim", "", ref_data)
  ref_data
  
  cat("\nReading in pruned variants\n")
  # clumped = list.files(path = paste(srcPath,"/qc/clumped/",sep=""), full.names = T, pattern = glob2rx("CLUMP*clumped"))
  clumped = list.files(path = paste(srcPath,"/qc/pruned/",sep=""), full.names = T, pattern = glob2rx("*prune.in"))
  # clumped = clumped[grepl(opt$out, clumped)]
 
  # Read in clumps
  # clumped.df = lapply(clumped, function(x) fread(x, header = T, stringsAsFactors = FALSE))
  clumped.df = lapply(clumped, function(x) fread(x, header = F, stringsAsFactors = FALSE, showProgress=FALSE))
  sum(unlist(lapply(clumped.df, nrow)))
  clumped.df = rbindlist(clumped.df)
  
  # read GWAS results
  gwasResults = sumstats
  names(gwasResults)[names(gwasResults) %in% opt$`snp-field`] = "INDEX"
  names(gwasResults)[names(gwasResults) %in% opt$pcol] = "P"
  cset = as.data.frame(clumped.df)
  colnames(cset) = "INDEX"
  cset = as.data.frame(merge(setDT(cset), setDT(gwasResults), by ='INDEX'))
  
  # Remove rows containing no clumped markers
  # clumped.df$SP2[clumped.df$SP2 %in% "NONE"] <- clumped.df$SNP[clumped.df$SP2 %in% "NONE"]
  # 
  # # Munge clumps
  # snp_list = clumped.df$SP2
  # snp_list = sapply(snp_list, function(x) strsplit(x, ","))
  # names(snp_list) = clumped.df$SNP
  # 
  # snp_friends = sapply(snp_list, function(x) length(x))
  # 
  # rep_function = function(x,y) {
  #   rep(x, y)
  # }
  # snp_reps = mapply(rep_function, x = names(snp_friends),y = snp_friends)
  # ld.df = data.frame(INDEX = unlist(snp_reps), SNP= unlist(snp_list))
  # split_char = strsplit(as.character(ld.df$SNP), "[(]")
  # ld.df$SNP = unlist(lapply(split_char, function(x) x[[1]]))
  # 
  # cat("....Calculating MAF\n")
  # 
  # freqs = list.files(path = paste(srcPath,"/qc/clumped", sep = ""), full.names = T, pattern = glob2rx("CLUMP*frq"))
  freqs = list.files(path = paste(srcPath,"/qc/pruned", sep = ""), full.names = T, pattern = glob2rx("*frq"))
  # freqs = freqs[grepl(opt$out, freqs)]
  # 
  freqs.df = lapply(freqs, function(x) fread(x, header = T, showProgress=FALSE))
  sum(unlist(lapply(freqs.df, nrow)))
  freqs.df = rbindlist(freqs.df)
  cols = c("SNP", "MAF")
  freqs.df = freqs.df[,cols,with = FALSE]
  # 
  # 
  cat("....Calculating LD scores for index markers\n")
  # 
  # write.table(clumped.df$SNP, file = paste(srcPath,"/qc/clumped/INDEX.snplist",sep=""), quote = F, col.names = F)
  fwrite(data.table(cset$INDEX), nThread = opt$threads, file = paste(srcPath,"/qc/clumped/INDEX.snplist",sep=""), quote = F, col.names = F, row.names = F)
  # 
  # print(paste("Calculating LD for index markers"))

  for( j in 1:length(ref_data)){

    cmd = paste("plink --bfile ", ref_data[[j]]," --ld-window 99999 --exclude ", srcPath,"/qc/ref/duplicatevars.txt --ld-snp-list ",srcPath,"/qc/clumped/INDEX.snplist --ld-window-kb ", opt$`ld-window`," --r2 --ld-window-r2 ", opt$r2, " --out ",srcPath,"/qc/clumped/LDINT_",opt$out,"_",j, sep = "")

    cmd = gsub("\n", "", cmd)
    cmd

    system(cmd, ignore.stdout = TRUE)

  }
  
  ld = list.files(path = paste(srcPath,"/qc/clumped",sep=""), full.names = T, pattern = paste("LDINT_",opt$out, sep = ""))
  ld = ld[!grepl(".log", ld)]
  
  ld.int = lapply(ld, function(x) fread(x, header = T, showProgress=FALSE))
  ld.int = rbindlist(ld.int)
  ld.int = as.data.frame(ld.int)
  
  # === calculate interval width (in base pairs)
  IntervalStart = ld.int[order(ld.int$BP_B, decreasing = F), ]
  IntervalStart = IntervalStart[!duplicated(IntervalStart$SNP_A), ]
  IntervalEnd = ld.int[order(ld.int$BP_B, decreasing =T), ]
  IntervalEnd = IntervalEnd[!duplicated(IntervalEnd$SNP_A), ]
  IntervalStart = IntervalStart[match(IntervalEnd$SNP_A, IntervalEnd$SNP_A), ]
  WIDTH = data.table(INDEX = IntervalStart$SNP_A, INT_WIDTH_KB = abs(IntervalStart$BP_B - IntervalEnd$BP_B)/1e3)
  
  ld.int = ld.int[,colnames(ld.int) %in% c("SNP_A", "CHR_A", "BP_A", "SNP_B", "R2")]
  colnames(ld.int) = c("CHR","BP","INDEX","SNP","R2")
  
  
  # ld.int = merge(cset[,c("INDEX", "SNP")], ld.int , by = "INDEX")

  ld.maf = merge(ld.int , freqs.df , by = "SNP")
  # 
  maf.bin = aggregate(ld.maf$MAF, by=list(ld.maf$INDEX), function(x) sum(-log(x)))
  colnames(maf.bin) = c("INDEX", "MAF")
  # 
  ld.df = merge(cset, maf.bin, by = "INDEX")
  ld.df = ld.df[!duplicated(ld.df), ]
  
  cat("....Removing xMHC from LD clumps\n")
  ld.df = ld.df[!ld.df$INDEX %in% rownames(mhc), ]
  ld.df$INDEX = as.factor(ld.df$INDEX)
  
  cat("....Adding interval width\n")
  ld.df = merge(WIDTH, ld.df, by = "INDEX")
  
  # 
  cat("....Finding LD friends for index variants\n")
  ld_friends = unlist(table(ld.int$INDEX))
  ld_friends = as.data.frame(ld_friends)
  ld_friends = ld_friends[match(ld.df$INDEX, ld_friends$Var1), ]
  ld.df$LD_FRIENDS = ld_friends$Freq
  ld.df$LD_FRIENDS[is.na(ld.df$LD_FRIENDS)] = 0
  
  cat("....Calculating LD scores for index variants\n")
  # 
  ld.score = aggregate(ld.int$R2, by = list(ld.int$INDEX), sum)
  colnames(ld.score) = c("INDEX", "LDSCORE")
  # ld.df$LOGP = -log10(ld.df$P)
  # ld.df$P[ld.df$P < 1e-15] <- 1e-15
  # ld.df$ZSCORE = qnorm(1-ld.df$P/2)
  # ld.df = merge(ld.df, ld.score, by="INDEX")
  
  # inverse-rank normalization
  my.invnorm = function(x){
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}
  
  ld.df = merge(ld.df, ld.score, by= "INDEX")
  ld.df$LOGP = -log10(ld.df$P)
  ld.df$P[ld.df$P < 1e-15] <- 1e-15
  ld.df$ZSCORE = qnorm(1-ld.df$P)
  ld.df$ZSCORE = my.invnorm(ld.df$ZSCORE)
  
  # 
  fwrite(setDT(ld.df), file = paste(srcPath,"/out/",opt$out,".sumstats",sep=""),quote=F,row.names=F,sep="\t",showProgress = FALSE, nThread = opt$threads)
  
  gwasMunged <<- ld.df 
  intervals <<- ld.int
  
  fwrite(setDT(ld.int), nThread = opt$threads , file = paste(srcPath,"/out/",opt$out,".intervals",sep=""),quote=F,row.names=F,sep="\t",showProgress = FALSE)
  
}


## FLEET annotation (via bedtools)

fleetAnnot = function(x) {
  
  # order snp into bed
  snp.bed = snp137_df[snp137_df$SNP %in% intervals$SNP, ]
  match.order = c('seqnames', 'start', 'start', 'SNP')
  snp.bed = snp.bed[,colnames(snp.bed) %in% match.order]
  snp.bed = snp.bed[,match(match.order, colnames(snp.bed))]
  snp.bed$start.1 = snp.bed$start + 1
  snp.bed  = snp.bed[order(snp.bed$seqnames, snp.bed$start), ]
  # snp.bed = data.frame(snp.bed, empty = NA, empty = "+")
  write.table(snp.bed, file = paste(srcPath,"/annotations/bedtools/snp.bed",sep=""), quote = F, row.names = F, col.names = F, sep = "\t")
  
  # === load annotations, run bedtools: 
  rdataAnnots = list.files(path = paste(srcPath,"/annotations/", sep=  ""), full.names = T, pattern = '.Rdata')
  rdataAnnots
  
  for( a in 1:length(rdataAnnots)){
    
    cat("\nAnnotating variants with bedtools:",rdataAnnots[[a]],"...")
    
    # Convert annotation file to bedtools format:
    annot.gr = readRDS(rdataAnnots[[a]])
    names(annot.gr) = NULL # ensure rownames are NULL (this can create problems if duplicates exist, which is highly likely)
    annot.df = as.data.frame(annot.gr)
    annot.bed = annot.df[,c("seqnames", "start", "end", "id")]
    annot.bed = annot.bed[order(annot.bed$seqnames, annot.bed$start), ]
    annot.bed = data.frame(annot.bed, empty = NA, empty = "+" )
    
    
    # === size of matrix
    castMb = round( nrow(annot.bed) * ncol(annot.bed) * 8/1e6, 3)
    nAnnots = length(unique(annot.bed$id))
    
    if(nAnnots >= 500){ cat("\n   Annotation matrix is approximately",castMb, "Mb in size and contains",nAnnots,"terms. Writing into chunks to reduce memory strain...\n")
      
      unique.tracks = unique(annot.bed$id)
      
      chunkSize = ceiling(length(unique(annot.bed$id))/10)
      
      trackChunks = split(unique.tracks, ceiling(seq_along(unique.tracks)/ chunkSize ))
      
      
      # === memory efficient operation for large annotation files: 
      for ( tc in 1:length(trackChunks)){
        
        write.table(annot.bed[annot.bed$id %in% trackChunks[[tc]], ], file = paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed",sep=""), quote= F, row.names = F, col.names = F, sep = "\t")
    
        # Run bedtools variant annotation:
        cmd = paste("bedtools intersect -a ",srcPath,"/annotations/bedtools/snp.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed -wb > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,"_",opt$out, ".txt" , sep = "")
        system(cmd, ignore.stdout = FALSE)
         }
       } else { write.table(annot.bed, file = paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed",sep=""), quote= F, row.names = F, col.names = F, sep = "\t")
        
        # Run bedtools variant annotation:
        cmd = paste("bedtools intersect -a ",srcPath,"/annotations/bedtools/snp.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed -wb > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])), "_",opt$out, ".txt" , sep = "")
        system(cmd, ignore.stdout = FALSE)
        
        }
      
    
    
    
  }
  
}
  
  
### Enrichment analysis
fleetTest = function() {
    
    
    search = ls()[grepl("gwasMunged", ls())]
    
    if(length(search) > 0) {ld.df <- gwasMunged} else{cat("\nReading summary statistics file...\n"); ld.df = as.data.frame(fread(paste(srcPath,"/out/",opt$out,".sumstats",sep=""), header = T, showProgress=FALSE)) }
    
    coef_trait = list()
    glm_trait  = list()
    
    rdataAnnots = list.files(path = paste(srcPath,'/annotations/bedtools/',sep=""), full.names = T, pattern = '.txt')
    rdataAnnots = rdataAnnots[grepl(paste(opt$out,".txt",sep=""), rdataAnnots)]
    rdataAnnots
    
    for( a in 1:length(rdataAnnots)){
      
      cat("\rLoading variant annotations:",rdataAnnots[[a]],"...")
      
      # read in and convert to wide format 
      
      bedtooled = fread(rdataAnnots[[a]], header = FALSE)
      bedtooled = bedtooled[,c(1,2,4,8)]
      bedtooled = bedtooled[!duplicated(bedtooled)]
      bedtooled$value = 1
      
      chr_annots = paste("chr",1:22,sep = "")
      
      cast_list = list()
      for(chr in 1:length(chr_annots)){
      cat("\n....Gathering annotations for chr:",chr)
      subed = bedtooled[bedtooled$V1 %in% chr_annots[[chr]]]
      if(nrow(subed) < 1) next
      cast = suppressMessages(dcast.data.table(subed, V4 ~ V8, value.var="value"))
      setnames(setDT(cast),"V4","SNP") # rename SNP column
      cast[is.na(cast)] = 0
      
      # === merge data table with ld.df stats
      tmp = merge(setDT(intervals[,c("SNP","INDEX")]), cast, by = "SNP")
      # 
      tmp = as.data.frame(tmp)
      tmp = tmp[,!colnames(tmp) %in% "SNP"]
      tmp = tmp[!duplicated(tmp$INDEX), ]
      
      cat("\n    Number of unique annotations detected:", ncol(cast)-1)
      # tmp[is.na(tmp)] = 0
      cast_list[[chr]] = setDT(tmp)
      }
      
      cast = rbindlist(cast_list, fill = TRUE)
      
      
      # === append missing snps 
      missing_snps = unique(ld.df$INDEX)
      missing_snps = missing_snps[!missing_snps %in% cast$INDEX]
      
      empty = matrix(0, nrow = length(missing_snps), ncol = ncol(cast)-1)
      empty_df = data.frame(SNP = missing_snps, empty)
      colnames(empty_df) = colnames(cast)
      
      cat("\n    Resolving missingness")
      full_annot = rbindlist(list(setDT(empty_df), setDT(cast)))
      full_annot = as.data.frame(full_annot)
      colnames(full_annot)[1] = "INDEX"
      
      # === merge data table with ld.df stats
      tmp = merge(ld.df, full_annot, by = "INDEX")
      tmp = as.data.frame(tmp)
      tmp[is.na(tmp)] = 0
     
      # === fix format of annotation labels 
      all_df = tmp
      
      new.tracks = sub("[[:punct:]]", "", colnames(all_df)[!colnames(all_df) %in% colnames(ld.df)])
      new.tracks = gsub("[-|,/;:'() ]", "_", colnames(all_df)[!colnames(all_df) %in% colnames(ld.df)])
      new.tracks = gsub("[+]", "_", new.tracks)
      new.tracks = paste("ID_", new.tracks, sep = "")
      colnames(all_df)[!colnames(all_df) %in% colnames(ld.df)] = new.tracks
      
      # === threshold for declaring genome-wide significace 
      gwas.sig.threshold = -log10(5e-08) # threshold for calling INDEX marker genome-wide significant
      
      
      # === check for genome-wide significant variants per annotation
      minColAnnots = ncol(ld.df)+1
      
      gw.sig.vars = all_df[all_df$LOGP >= gwas.sig.threshold, ]
      if(ncol(gw.sig.vars) > minColAnnots){gw.sig.vars.hits = colSums(gw.sig.vars[,minColAnnots:ncol(gw.sig.vars)])} else {gw.sig.vars.hits = sum(gw.sig.vars[,minColAnnots]); names(gw.sig.vars.hits) = colnames(gw.sig.vars)[minColAnnots]}
      
      # ===  calculate frequency of 1s in columns  
      
      if( ncol(all_df) > minColAnnots){freq_annot = colSums(all_df[, minColAnnots:ncol(all_df)])} else{freq_annot = sum(all_df[,minColAnnots]); names(freq_annot) = colnames(all_df)[[minColAnnots]]}
      
      # minimum annotation frequency
      minAnnotFreq = opt$`annot-cnt`
      
      low_freq = freq_annot[freq_annot < minAnnotFreq]
      
      cat("\nTotal number of annotations:", length(freq_annot))
      
      # === annotation QC step:
      removed_annots = length(low_freq)
      cat("\nRemoved",removed_annots,"annotations due to low frequency: (<",minAnnotFreq,"counts)\n")
      
      ## Enrichment analysis
      all_df = as.data.frame(all_df)
      all_df = all_df[,!colnames(all_df) %in% names(low_freq)]
      all_df = data.frame(inv_maf = 1/all_df$MAF, inv_ldscore = 1/all_df$LDSCORE, all_df)
      all_df$inv_maf = log(all_df$MAF)
      all_df$inv_ldscore = log(all_df$LDSCORE)
      
          tracks = colnames(all_df)[!colnames(all_df) %in% colnames(ld.df)]
          tracks = tracks[!tracks %in% c("inv_ldscore", "inv_maf")]
          
          cat("\nTotal number of LD-intervals available for enrichment tests:", nrow(all_df),"\n")
          
          if(length(tracks) < 1) next
          
          # === FLEET enrichment analysis via linear regression
          fits = list()
          rsq.fits = list()
          
          
          SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
          SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
          
                
                ## === New feature (October 3, 2017)
                if(opt$speed == 'fast'){
                  
                    dt_all_df = data.table(all_df)
                    
                    misCol = colSums(dt_all_df[,colnames(dt_all_df) %in% tracks,with=F] == 1)
                    
                    misCol = misCol[misCol < opt$`annot-cnt`]
                    
                    if(length(misCol) > 0){
                      
                      dt_all_df = dt_all_df[,!colnames(dt_all_df) %in% names(misCol),with=F]}
                    
                    response_matrix = as.matrix(dt_all_df[,colnames(dt_all_df) %in% tracks,with=F])
                    mod = model.matrix( ~ ZSCORE + LDSCORE + MAF + inv_ldscore + inv_maf + INT_WIDTH_KB, dt_all_df)
                    
                    cat("\r...Running linear regression in fast mode!\n")
                    start_reg = proc.time()
                    lmFast = lm(response_matrix ~ mod)
                    cat("\r...Time to complete regression:", (proc.time() - start_reg)[3], 'seconds')
                    lmFastCoef = summary(lmFast)
                    
                    if(ncol(response_matrix) > 1){
                    lmCoef = lapply(lmFastCoef, function(x) broom::tidy(x$coefficients))
                    names(lmCoef) = tracks
                    lmDf = ldply(lmCoef)} else{
                    lmCoef = broom::tidy(summary(lmFast)$coefficients)
                    lmDf = data.frame(SET_NAME = tracks, lmCoef)}
                    
                    
                    lmDf = lmDf[lmDf$.rownames %in% "modZSCORE", ]
                    colnames(lmDf) = c("SET_NAME", "Term", "Beta", "SE", "T","P")
                    lmDf = lmDf[,!colnames(lmDf) %in% "Term"]
                    
                    lmDf$SET_SOURCE = SETNAME
                    
                    gwVars = gw.sig.vars[,colnames(gw.sig.vars) %in% c("INDEX",tracks)]
                    gwVars = suppressMessages(melt(gwVars))
                    
                    if(nrow(gwVars) > 0){gwVars = gwVars[!gwVars$value == 0, ]}
                    
                    
                    # === annotate lm stats for annotations with genome-wide significant variants
                    if(nrow(gwVars) > 0){
                      
                    gwVars_agg = aggregate(gwVars$INDEX, by=list(gwVars$variable), function(x) paste(x, collapse="|",sep=""))
                    colnames(gwVars_agg) = c("SET_NAME","GWS_VAR_HITS")
                    gwVars_agg = gwVars_agg[match(lmDf$SET_NAME, gwVars_agg$SET_NAME), ]
                    gwVars_agg$SET_NAME = lmDf$SET_NAME
                    gwVars_agg$GWS_VAR_HITS[is.na(gwVars_agg$GWS_VAR_HITS)] = "*"
                    lmDf$GWS_VAR_HITS = gwVars_agg$GWS_VAR_HITS 
                    } else {lmDf$GWS_VAR_HITS = "*"}
                    
                    lmDf$SET_NAME = gsub("ID_", "", lmDf$SET_NAME)
                    lmDf$SET_NAME = gsub("_",".", lmDf$SET_NAME)
                    coef_df = lmDf
                    
                    }
                
          if(opt$speed == "slow"){
          for( x in 1:length(tracks)){
            
          cat('\r....Fitting linear models:',x)
          #model = lm(ZSCORE ~ as.matrix(all_df[,colnames(all_df) %in% c("LDSCORE",  "MAF", "inv_ldscore", "inv_maf", tracks[[x]])]), data = all_df, weights = (all_df$MAF^-1))
          model = lm(ZSCORE ~ as.matrix(all_df[,colnames(all_df) %in% c("LDSCORE",  "MAF", "inv_ldscore", "INT_WIDTH_KB", "inv_maf", tracks[[x]])]), data = all_df)
          res = summary(model)$coefficients 
          tidy = tidy(res)
          names = strsplit(tidy$.rownames, "])")
          names = lapply(names , function(x) x[[length(x)]])
          tidy = data.frame(predictors = unlist(names), tidy)
          
          gwas_vars = gw.sig.vars[,c(tracks[[x]], "INDEX", "P")]
          gwas_vars = ifelse(gwas_vars[,1] == 1, as.character(gwas_vars$INDEX), "")
          gwas_vars = gwas_vars[!gwas_vars %in% ""]
          gwas_vars = paste(gwas_vars, collapse = "|", sep = "")
          
          SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
          SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
          
          tidy = tidy[tidy$predictors %in% tracks[[x]], ]
          tidy = tidy[,!colnames(tidy) %in% c(".rownames")]
          tidy$GWSVARNAMES = ifelse( gwas_vars == "", "*", gwas_vars)
          tidy$SET_SOURCE = SETNAME
          
          fits[[x]] = tidy
          
          rsq.fits[[x]] = summary(model)$r.squared
          
          if(opt$robust == TRUE){
            
          hcmodel = vcovHC(model, omega = NULL, type = "HC4")
          hcFit = coeftest(model, df = Inf, hcmodel) # adjust for heteroskedasticity w/ robust standard errors
          hcTidy = tidy(hcFit)
          
          names = strsplit(hcTidy$term, "])")
          names = lapply(names , function(x) x[[length(x)]])
          hcTidy = data.frame(predictors = unlist(names), hcTidy)
          
          hcTidy = hcTidy[hcTidy$predictors %in% tracks[[x]], ]
          hcTidy = hcTidy[,!colnames(hcTidy) %in% c("term")]
          hcTidy$GWASVARNAMES = tidy$GWSVARNAMES
          hcTidy$SET_SOURCE = SETNAME
          colnames(hcTidy) = colnames(tidy)
          
          fits[[x]] = hcTidy # overwrite slot with robust se
          }
          
        
          }
          
          # Combine model results
          coef_df = ldply(fits)
          coef_df$rsq = unlist(rsq.fits) # Multiple R2 
          coef_df$annot_cnt = freq_annot[names(freq_annot) %in% coef_df$predictors]
          coef_df$annot_frq = coef_df$annot_cnt/nrow(all_df)
          coef_df$GWS_VAR_HITS = gw.sig.vars.hits[names(gw.sig.vars.hits) %in% coef_df$predictors]
          }
          
          # Permutation enrichment analysis
          all_df = as.data.frame(all_df)
          
          if(opt$`fleet-permutation` == TRUE){
      
          if(opt$`fast-permutation` == TRUE){
            # Randomize annotations, count number of variants below P-value threshold with annotation
            cat("\nRunning permutation enrichment analysis in fast mode\n")
    
              gws.hits.1 = all_df$INDEX[all_df$P < 5e-08]
              gws.hits.2 = all_df$INDEX[all_df$P < 1e-06]
              gws.hits.3 = all_df$INDEX[all_df$P < 1e-04]

              # gw.sig.list = list(gws.hits.1, gws.hits.2,gws.hits.3)
              # names(gw.sig.list) = c(5e-08, 1e-06, 1e-04)
              
              # Number of variants with corresponding P-threshold in annotation
              if(ncol(all_df) > minColAnnots){
              pt1_hits = colSums(all_df[all_df$P < 5e-08, minColAnnots:ncol(all_df)])
              pt2_hits = colSums(all_df[all_df$P < 1e-06, minColAnnots:ncol(all_df)])
              pt3_hits = colSums(all_df[all_df$P < 1e-04, minColAnnots:ncol(all_df)])
              } else {
              pt1_hits = sum(all_df[all_df$P < 5e-08, ncol(all_df)])
              pt2_hits = sum(all_df[all_df$P < 1e-06, ncol(all_df)])
              pt3_hits = sum(all_df[all_df$P < 1e-04, ncol(all_df)])
              }
              
              # randomize all_df nPerm times
              cat("\n")
              
              random_pt1_hits = list()
              random_pt2_hits = list()
              random_pt3_hits = list()
              
              for(x in 1:opt$nPerms){
                cat("\r    Permutation:",x)
                
                randomize = as.data.frame(all_df)
                randomize$LOGP = sample(randomize$LOGP)
                
                if(ncol(randomize) > 2){
                random_pt1_hits[[x]] = colSums(randomize[randomize$LOGP >= gwas.sig.threshold, minColAnnots:ncol(randomize)])
                random_pt2_hits[[x]] = colSums(randomize[randomize$LOGP >= gwas.sig.threshold, minColAnnots:ncol(randomize)])
                random_pt3_hits[[x]] = colSums(randomize[randomize$LOGP >= gwas.sig.threshold, minColAnnots:ncol(randomize)])} else {
                  random_pt1_hits[[x]] = sum(randomize[randomize$LOGP >= gwas.sig.threshold, ncol(randomize)])
                  random_pt2_hits[[x]] = sum(randomize[randomize$LOGP >= gwas.sig.threshold, ncol(randomize)])
                  random_pt3_hits[[x]] = sum(randomize[randomize$LOGP >= gwas.sig.threshold, ncol(randomize)])
                }
                
              }
              random_pt1_hits = do.call(cbind, random_pt1_hits)
              random_pt2_hits = do.call(cbind, random_pt2_hits)
              random_pt3_hits = do.call(cbind, random_pt3_hits)
              
              
              cat("\r    Calculating empirical p-values")
              perm_out = list()
              for(x in 1:length(pt1_hits)){
                
                pval_perm1 = sum(random_pt1_hits[x, ] >= pt1_hits[[x]])/opt$nPerms
                if(pval_perm1 == 0){pval_perm1 = 1/opt$nPerms}
                pval_perm2 = sum(random_pt2_hits[x, ] >= pt2_hits[[x]])/opt$nPerms
                if(pval_perm2 == 0){pval_perm2 = 1/opt$nPerms}
                pval_perm3 = sum(random_pt3_hits[x, ] >= pt3_hits[[x]])/opt$nPerms
                if(pval_perm3 == 0){pval_perm3 = 1/opt$nPerms}
                
                perm_out[[x]] = data.frame(SET_NAME = names(pt1_hits)[[x]], EmpP_pt1 = pval_perm1, EmpP_pt2 = pval_perm2, EmpP_pt3 = pval_perm3, OBSCNT = pt1_hits[[x]], MEANOBSP = mean(random_pt1_hits[x,]))
                
              }
              perm_out_df = ldply(perm_out)
              
              
              cat("\n...Compiling enrichment statistics")
              coef_df = merge(coef_df, perm_out_df, by='SET_NAME',all=T)
              
            }
            
            
              if(opt$`robust-permutation` == TRUE){
                
                cat("\nRunning permutation analysis in robust mode!\n")
                
                ## Run permutation (GW-SIG):
                # pval.check = which(coef_df$`Pr(>|t|)` <= opt$permStartThreshold & coef_df$Estimate > 0)
                
                dt = setDT(all_df)
                
                gc()
    
                  # subset SNPs by significance thresold
                  if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(5e-08, 1e-07, 1e-06)}
                  
                  # how many permutations to run per annotation
                  nPerms = opt$nPerms
                  
                  ptPermDf = list()
                  for(pt in 1:length(pthreshold)){
                    
                  gws.hits.1 = dt[dt$P < pthreshold[[pt]]]
                 
                  
                  # calculate number of "hits" among top SNPs
                  if(length(tracks) > 1){hits1 = colSums(gws.hits.1[,tracks,with=FALSE])}
                  if(length(tracks) == 1){hits1 = sum(gws.hits.1[,tracks,with=FALSE]); names(hits1) = tracks}
                  
                  hits1 = hits1[hits1 > 0]
                  
                  if(length(hits1) >= 1){
                  
                  sub = as.data.frame(gws.hits.1[,c(names(hits1), "INDEX"), with=FALSE])
                  sub = suppressMessages(melt(sub))
                  sub = sub[sub$value == 1, ]
                  sub = merge(setDT(sub), gws.hits.1[,c("INDEX", "MAF"),with=FALSE], by= "INDEX")
                  
                  # calculate avg minor allele frequency for SNPs below Pt in annotation
                  avg.maf = aggregate(sub$MAF, by = list(sub$variable), mean)
                  colnames(avg.maf) = c("variable", "mean.maf")
                  sd.maf = aggregate(sub$MAF, by = list(sub$variable), sd)
                  size.bin = aggregate(sub$MAF, by = list(sub$variable), length)
                  avg.maf$SE = sd.maf$x/sqrt(size.bin$x)
                  avg.maf$SE[is.na(avg.maf$SE)] = 0
                  avg.maf$CI = 1.96 * sd.maf$x/sqrt(size.bin$x)
                  avg.maf$CI[is.na(avg.maf$CI)] = 0
                  avg.maf$CI_high = avg.maf$mean.maf + avg.maf$CI
                  avg.maf$CI_low = abs(avg.maf$mean.maf - avg.maf$CI)
                  avg.maf$variable = as.character(avg.maf$variable)
                  
              
                  gc()
                  
                  permPval = list()
                  for(v in 1:nrow(avg.maf)){
                    start.clock = proc.time()
                    cat("\r   Permutation analysis with threshold P-value < ", pthreshold[[pt]] ,":",round(v/nrow(avg.maf)*100, 2),"% complete")
                    
                    
                    permutable = dt[!dt$INDEX %in% gws.hits.1$INDEX & dt$MAF >= floor(avg.maf$CI_low[[v]]) & dt$MAF <= ceiling(avg.maf$CI_high[[v]]), c("LOGP", "INDEX", "MAF", avg.maf$variable[[v]]),with=FALSE]
                    
                    total_size = nrow(gws.hits.1)
                    hit_ratio = (hits1[names(hits1) %in% avg.maf$variable[[v]]]/total_size)[[1]]
                    
                    if(nrow(permutable) > total_size){ 
                    
                    permSum = list(rep(NA, nPerms))
                    for(p in 1:nPerms){
                    tmp = permutable[sample(nrow(permutable), total_size)]
                    permSum[[p]] = sum(tmp[,4])/nrow(tmp)
                    rm(tmp) # remove vector
                    }
                    
                    pval = sum(unlist(permSum) >= hit_ratio)/nPerms} else {pval = NA}
                    
                    
                    permPval[[v]] = pval
                    
                    secs = proc.time() - start.clock
                  }
                 
                  colHeaderA = paste("Emp_P",pt,sep="")
                  permDf = data.frame(SET_NAME = avg.maf$variable, Pval = unlist(permPval))
                  colnames(permDf)[2] = colHeaderA
                  
                  cat("\n")
                  
                  # ptPermDf[[pt]] = permDf
                  
                  coef_df = merge(coef_df, permDf, by = "SET_NAME", all.x = TRUE)
                  
                  }
                  
                  }
              
                  

                  # perm_sets_df = Reduce(function(x,y) merge(x,y,by='predictors',all=TRUE), ptPermDf)
          
                  # if(nrow(perm_sets_df) > 0){coef_df = merge(coef_df, perm_sets_df, by = "predictors", all = TRUE)}

          
          } # close if (pval check)

        } # close if (slow permutation)
            
        cat("\nWriting results...\n")
        # coef_df = data.frame(predictors = coef_df$predictors, coef_df[,!colnames(coef_df) %in% "predictors"])
        # names(coef_df)[names(coef_df) %in% c("estimate", "std.error", "statistic", "p.value", "rsq")] = c("Beta", "SE", "Tstat", "Pval", "ModelR2")
        # coef_df$predictors = substring(coef_df$predictors, 4, nchar(as.character(coef_df$predictors))) # Munge summary statistics (enrichment tests)
        
        if(a == 1) { write.table(coef_df, file = paste(srcPath,"/out/",opt$out,"_fleetEnrichment.txt",sep=""), quote = F, row.names = F, sep = "\t")} else{ write.table(coef_df, file = paste(srcPath,"/out/",opt$out,"_fleetEnrichment.txt",sep=""), quote = F, row.names = F, col.names = F, sep = "\t", append = TRUE)}
          
      
    } # close for loop (run next annotation set)
    
} # close fleetTest()
  

### Plotting mechanisms:

fleetPlot = function(x) {
  
####### start qq plot ########
par(mfrow=c(1,1))

coef_df = read.table(paste(srcPath,"/out/",opt$out,"_fleetEnrichment.txt",sep=""), header = T, fill = NA)

coef_df  = coef_df[coef_df$Estimate > 0 & !is.na(coef_df$Estimate), ]

observed <- sort(coef_df$Pr...t..)
lobs <- -(log10(as.numeric(observed)))
lobs[lobs > 8] <- 8

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

pvalues = coef_df$Pr...t..
nn<-sapply(pvalues, length)
rs<-cumsum(nn)
re<-rs-nn+1
n<-min(nn)

conf.alpha = .05
conf.points = length(pvalues)
# conf.points = 10e3
n = length(pvalues);
# n = 1e03
mpts<-matrix(nrow=conf.points*2, ncol=2)
for(i in seq(from=1, to=conf.points)) {
  mpts[i,1]<- -log10((i-.5)/n)
  mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
  mpts[conf.points*2+1-i,1] <- -log10((i-.5)/n)
  mpts[conf.points*2+1-i,2] <- -log10(qbeta(conf.alpha/2, i, n-i))
}

title= "Quantile-Quantile Plot"

medianchi = median(qchisq(coef_df$Pr...t.., df = 1,lower.tail = F), na.rm = T)
lambda = medianchi/0.454

plot (0,0,type="n", xlim=c(0,max(lexp)), ylim=c(0,max(lobs)), xlab="Expected -log10 (P)",
      main = title,
      cex.lab =1.0, cex.axis=1.15,
      ylab = "Observed -log10 (P)", las=1 )
polygon(x = mpts[,1], y=mpts[,2], col="gray90", border = "grey", lwd = 3)
lines(c(0,8), c(0,8), col="black", lwd=1.5)

points(x = lexp, y = lobs, pch = 8, cex = .5, col = "black")

text(x = 2, y= 0.5, cex = 1, labels = substitute(paste("Median ", chi^2, " = ", MC), list(MC=format(medianchi, digits = 3))))
text(x = 2,   y = .25, cex = 1, labels = substitute(paste(lambda, " = ", LM), list(LM=format(lambda,digits = 4))))

####### end qq plot ########
}
  



### Run operations

total_clock_start = proc.time()

readSNPcoords() # always run (from 1000G)

xMHCregion() # always run

fleetReadGWAS()

# readAnnots() # always run

if(opt$`fleet-prune-ref` == TRUE){
  refPrune()
}


if(opt$`fleet-annotate` == TRUE){
  fleetLD() 
  fleetAnnot()
}

if(opt$`fleet-enrichment` == TRUE & opt$`fleet-annotate` == FALSE){
  
  fleetLD() 
  fleetTest()
  
}

if(opt$`fleet-enrichment` == TRUE & opt$`fleet-annotate` == TRUE){
  
  fleetTest()
  
}

if(opt$plots == TRUE){
  fleetPlots()
}

time_elap = proc.time()
time_elap = abs(total_clock_start[[1]] - time_elap[[1]])

cat("\nTotal time elapsed:", time_elap[[1]], "seconds\n")
 
#close(log_con)
