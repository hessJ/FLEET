#!/usr/bin/env Rscript

# global options
options(digits = 5)

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
||
|| Functional LD-interval EnrichmEnt Test (FLEET)
||
|| Jonathan L. Hess, PhD and Stephen J. Glatt, PhD (c) 2017
||
|| SUNY Upstate Medical University, PsychGENe Lab
||
|| Contact: hessjo@upstate.edu
||
|| https://github.com/hessJ/FLEET
||
|| GNU GENERAL PUBLIC LICENSE v3
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
  
  make_option(c("-R", "--r2"), type="double", default=0.5, 
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
  
  make_option(c("-F", "--annot-cnt"), type="double", default=100, 
              help="Minimum annotation count observed across LD-clumps [default = %default]", metavar="double"),
  
  make_option(c("-B", "--boot"), type="double", default=10, 
              help="Percent of intervals sampled for bootstrap regressions [default = %default%]", metavar="double"),
  
  ## Run options....
  
  # make_option(c("-M", "--fleet-prune-ref"), type="logical", default = TRUE, 
  #             help="Initiate LD-pruning step of 1KG reference data. Only needs to be run once. [default = %default]", metavar="logical"),
  
  make_option(c("-A", "--fleet-annotate"), type="logical", default=TRUE, 
              help="Annotate LD-clumps with bedtools [default = %default]", metavar="logical"),
  
  make_option(c("-E", "--fleet-enrichment"), type="logical", default=TRUE, 
              help="Perform enrichment analysis with weighted linear models [default = %default]", metavar="logical"),
  
  make_option("--fleet-permutation", type="logical", default=TRUE, 
              help="Perform enrichment analysis with permutation (randomizing annotations) [default = %default]", metavar="logical"),
  
  # make_option("--fast-permutation", type="logical", default=FALSE, 
  #             help="Simple permutation analysis [default = %default]", metavar="logical"),
  
  make_option("--robust-permutation", type="logical", default=TRUE, 
              help="Permutation analysis that will sample variants from the MAF bin of target SNPs [default = %default]", metavar="logical"),
  
  make_option("--speed", type="character", default='fast', 
              help="Change behavior of linear models (fast mode: SET becomes response variable, slow mode: Z-score becomes response variable) [default = %default]", metavar="character"),
  
  make_option("--pthres", type="character", default="",
              help="Table with P-value threshold(s) for SNP bins [default P-values < 5e-08, 1e-07, and 1e-06]", metavar="character"),
  
  make_option(c("--plink"), type="character", default = "", 
              help="Path to plink executable (v1.90) file [default = %default]", metavar="character"),
  
  make_option(c("--bedtools"), type="character", default = "", 
              help="Path to bedtools executable [default = %default]", metavar="character"),
  
  make_option("--plots", type="logical", default=FALSE,
              help="Option to auto-generate QQ-plots for genome-wide enrichment tests [default = %default]", metavar="logical")
  
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
suppressWarnings(suppressMessages(library(sfsmisc)))
suppressWarnings(suppressMessages(library(locfdr)))
suppressWarnings(suppressMessages(library(SumVg)))
cat("Successfully loaded packages!\n")

srcPath="~/Documents/FLEET/"
opt$gwas="~/Google Drive/mac_storage/PGC_results/scz2014_full.txt"
opt$out="scz2014"
opt$`snp-field`="snpid"
opt$pcol="p"
opt$bedtools="~/Documents/bin/bedtools2/bin/"
opt$plink ="~/Documents/plink_mac/"

# QQPlot method

qqPlot = function(x){
  
  chisq1 <- qchisq(1-x,1)
  
  medianchi = median(chisq1,na.rm=T)
  lambdaOut = medianchi/.454
  
  observed <- sort(x)
  lobs <- -(log10(as.numeric(observed)))
  
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  ci = .95
  N = length(x)
  observed = -log10(sort(x))
  expected = -log10(1:N / N)
  clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1))
  cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
  
  set = data.frame(expected = expected, observed = observed, clower = clower, cupper = cupper)
  
  qplot = ggplot(set, aes(x = lexp, y = lobs)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, col = 'black', lwd = 0.3) +
    theme_classic() + 
    ylab(expression(paste("Observed -",italic(log)[10]," P-value"))) +
    xlab(expression(paste("Expected -",italic(log)[10]," P-value"))) +
    scale_fill_discrete(NULL) +
    theme(text = element_text(size = 20)) +
    geom_line(aes(x = expected, y = clower), colour ='grey' ,lwd = 0.75) +
    geom_line(aes(x = expected, y = cupper), colour = 'grey', lwd = 0.75) +
    geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), fill="grey", alpha="0.2") 
  
  png(file=paste(srcPath,"/plots/QQplot.png",sep=""),res=300,units="in",height=7,width=7)
  print(qplot)
  dev.off()
  
}

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
  if(infer == 1){cat("\nGuessing that GWAS file is tab-delimited\n"); read.gwas = fread(opt$gwas, sep = "\t", showProgress=FALSE)}
  if(infer == 2){cat("\nGuessing that GWAS file is space-delimited\n"); read.gwas = fread(opt$gwas,  sep = " ", showProgress=FALSE)}
  if(infer == 3){cat("\nGuessing that GWAS file is comma-delimited\n"); read.gwas = fread(opt$gwas,  sep = ",", showProgress=FALSE)}
  read.gwas = as.data.frame(read.gwas)
  
  # freq from reference
  freq = list.files(path = paste(srcPath,"/qc/ref/",sep=""), pattern = ".frq", full.names = T)
  freq = lapply(freq , function(x) fread(x))
  freq = rbindlist(freq)
  
  # Align to reference panel
  cat("\rAssuming OR > 1 --> Risk (phenotype) increasing")
  cat("\rMatching alleles to reference panel")
  read.gwas = data.table(read.gwas)
  # rename headers
  
  names(read.gwas)[names(read.gwas) %in% c("or","OR","odds","Odd")] = "OR"
  names(read.gwas)[names(read.gwas) %in% c("se","SE","StdErr","SDERR","STDERROR")] = "se"
  
  names(read.gwas)[names(read.gwas) %in% c(opt$`snp-field`, opt$pcol)] = c("SNP","P")
  # retain only markers found in reference panel
  read.gwas = read.gwas[read.gwas$SNP %in% freq$SNP]
  # recode effect sizes and a1
  read.gwas$ORc = ifelse(read.gwas$OR < 1, 1/read.gwas$OR, read.gwas$OR)
  read.gwas$A1c = ifelse(read.gwas$OR > 1, read.gwas$a1, read.gwas$a2)
  
  freq = freq[match(read.gwas$SNP, freq$SNP)]
  match0 = freq$A1 == read.gwas$A1c
  match1 = freq$A2 == read.gwas$A1c
  
  co = data.frame(match0,match1)
  bad = which(rowSums(co) == 0)
  
  cat("\rRemoving",length(bad),"variants with mis-matched alleles between reference panel and GWAS")
  if(length(bad) > 0){
    freq = freq[-bad];
    read.gwas = read.gwas[-bad]
  }
  
  # match alleles
  match2 = freq$A1 == read.gwas$A1c
  
  freq$A1c = ifelse(match2 == F, freq$A2, freq$A1)
  freq$MAFc = ifelse(match2 == F, 1-freq$MAF, freq$MAF)
  
  read.gwas$MAFc = freq$MAFc
  
  infoCol = which(grepl("info|INFO|imputation", colnames(read.gwas)))
  
  if(length(infoCol) > 0){
    cat("\rDetected imputation qualtiy column as:",colnames(read.gwas)[infoCol])
    preinfo = nrow(read.gwas)
    postinfo = read.gwas[,colnames(read.gwas) %in% colnames(read.gwas)[infoCol],with=F]
    keepinfo = which(postinfo > .8)
    gwas.info = read.gwas[keepinfo]
    cat("\rRemoved",abs(nrow(gwas.info)-preinfo),"variants with imputation quality < 0.9")
  }
  
  # remove SNPs with chi-square > 30
  chisq = (log(gwas.info$OR)/gwas.info$se)**2
  outlier = which(chisq > 30)
  
  if(length(outlier) > 0){
    cat("\rRemoved",length(outlier),"variants with chi-square > 30")
    gwas.info = gwas.info[-outlier]
  }
  
  # set population prevalence
  K = 0.01
  RR1 = gwas.info$ORc
  RR2 = gwas.info$ORc^2
  
  PA = gwas.info$MAFc
  Paa = (1-PA)^2
  PAa = 2*PA*(1-PA)
  PAA = PA^2
  muaa=0
  faa= K/(Paa + PAa*RR1 + PAA*RR2)
  fAa= RR1*faa
  fAA= RR2*faa 
  T = qnorm(1-faa) 
  muAa = T-qnorm(1-fAa)
  muAA = T-qnorm(1-fAA)
  mean.all= PAa*muAa+ PAA*muAA
  Vg = Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
  actual.Vg =  Vg/(1+Vg) 
  
  N = 70e3
  chisq  = (log(gwas.info$OR)/gwas.info$se)**2
  VgM = chisq/(N-2+chisq)
  
  gwas.info$Vg = actual.Vg
  gwas.info$R2m = VgM
  
  fwrite(data.table(gwas.info$SNP),
         file = paste(srcPath,"/qc/gwas.snplist",sep=""),
         quote = F, row.names = F, col.names=F, sep = "\t")
  
  sumstats <<- data.frame(gwas.info)
  
}


### Prune 1KG reference 
refPrune = function(){

  ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
  ref_data = gsub(".bim", "", ref_data)

  if(opt$plink != ""){PathToPlink = paste(opt$plink,sep="/")} else{PathToPlink = ""}
  
  for( i in 1:length(ref_data)){
    cat("\rPruning progress:",round(i/length(ref_data)*100,3),"%")
    indep.cmd = paste('--indep-pairwise 50 5 0.2')
    # cat("\nPruning variants to approximate linkage equilibrium ( Plink command:",indep.cmd,"): ",i)
    prune = paste(PathToPlink,"plink --bfile ", ref_data[[i]], " --seed 1234 --maf .05 --freq --extract ", srcPath,"/qc/gwas.snplist ",indep.cmd," --out ", srcPath,"/qc/pruned/PRUNED_",basename(ref_data[[i]]), sep = "")
    system(prune, ignore.stdout = TRUE)
  }

}


## Clump GWAS data 

fleetClump = function() {

  ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
  ref_data = gsub(".bim", "", ref_data)
  
  
  PRUNEIN = list.files(path = paste(srcPath,'/qc/pruned/',sep=""), full.names = T, pattern = '.prune.in')
  PRUNEIN
  
  readPruned = lapply(PRUNEIN, function(x) fread(x, header = FALSE))
  readPruned = rbindlist(readPruned)
  
  # keepPrune = list.files(path = paste(srcPath,'/qc/pruned/',sep=""), full.names = T, pattern = 'prune.in')

  cat("\n##### Creating LD-intervals from GWAS results #####\n")
  # cat("\nUsing", opt$threads, "threads\n")
  
  fwrite(data.table(readPruned),
         file=paste(srcPath,"/qc/INDEX.txt", sep=""),
         quote  = F, row.names = F, sep ="\t", col.names = FALSE)
  
  if(opt$plink != ""){PathToPlink = paste(opt$plink,sep="/")} else{PathToPlink = ""}
  
  for( j in 1:length(ref_data)){
    cat("\rClumping progress:",round(j/length(ref_data),3)*100,"%")
    
    cmd = paste(PathToPlink,"plink --bfile ", ref_data[[j]]," --ld-window 99999 --exclude ", srcPath,"/qc/ref/duplicatevars.txt --ld-snp-list ",srcPath,"/qc/INDEX.txt --ld-window-kb ", opt$`ld-window`," --r2 --ld-window-r2 ", opt$r2, " --out ",srcPath,"/qc/clumped/LDINT_",opt$out,"_",j, sep = "")
    # run plink command
    system(cmd, ignore.stdout = TRUE)
    # remove log files from directory
    suppressMessages(file.remove(paste(srcPath,"/qc/clumped/LDINT_",opt$out,"_",j,".log", sep = "")))
    
  }
  
  
  ld = list.files(path = paste(srcPath,"/qc/clumped",sep=""), full.names = T, pattern = paste("LDINT_",opt$out, sep = ""))
  ld = ld[!grepl(".log", ld)]
  
  ld.int = lapply(ld, function(x) fread(x, header = TRUE, showProgress=FALSE))
  ld.int = rbindlist(ld.int)
  names(ld.int)[names(ld.int) %in% "CHR_A"] = "CHR"
  
  # == ADD P-values to intervals 
  sumstats = data.table(sumstats)
  sumstats = sumstats[match(ld.int$SNP_B, sumstats$SNP)]
  ld.int$Pval = sumstats$P
  
  # == Position 
  indexPos = ld.int[,c("SNP_A","CHR"),with=F]
  indexPos = indexPos[!duplicated(indexPos)]
  
  
  # === calculate interval width (in base pairs)
  IntervalRange = ld.int[,list(IntervalStart=min(BP_B), IntervalEnd = max(BP_B)),by=(SNP_A)]
  IntervalRange$LOCSIZE = abs(IntervalRange$IntervalStart - IntervalRange$IntervalEnd)
  LDFRIENDS = ld.int[,list(NSNP = length(SNP_B)),by=(SNP_A)]
  LDSCORE = ld.int[,list(LDSCORE = sum(R2)),by=(SNP_A)]
  
  # == Calculate mean P-value for interval
  IntervalPval = ld.int[,list(meanPval = mean(Pval,na.rm=TRUE)),by=(SNP_A)]
  IntervalPval$MeanZSCORE = qnorm(1-IntervalPval$meanPval)
  
  # === Combine Interval stats with sumstats 
  int.stats = data.table(IntervalRange, LDFRIENDS,LDSCORE, IntervalPval, indexPos)
  int.stats = int.stats[,!duplicated(colnames(int.stats)),with=F]
  colnames(int.stats)[1] = "SNP"
  int.stats$LOCSIZE[int.stats$LOCSIZE == 0] = 1
  int.stats$SNPDENSITY = int.stats$NSNP/(int.stats$LOCSIZE/1e3)
  
  ld.int = merge(data.table(sumstats), int.stats, by= "SNP")
  ld.int$ZSCORE = qnorm(1-ld.int$P)
  
  ld.int = ld.int[!ld.int$NSNP == 1]
  
  cat("\rDetected",nrow(ld.int),"intervals for analysis")
  cat("\rDetected",nrow(ld.int[ld.int$P < 5e-08]),"genome-wide significant intervals! Some may have been removed by pruning")
  
  
  # save interval sumstats
  fwrite(data.table(ld.int),
         file = paste(srcPath,"/out/",opt$out,".sumstats",sep=""),
         quote=F,row.names=F,
         sep="\t",showProgress = FALSE, nThread = opt$threads)
  
  # make sumstats available outside function
  gwasMunged <<- ld.int 
  
# 
#   if(opt$plink != ""){PathToPlink = paste(opt$plink,sep="/")} else{PathToPlink = ""}
#   
#   if(opt$threads > 1){
# 
#     foreach(j = 1:length(ref_data)) %dopar% {
#       # ensure that spaces are escaped in GWAS file path
#       gwas_file = gsub(" ", "\\ ", opt$gwas, fixed = T)
#       cmd = paste(PathToPlink,"plink --bfile  ", ref_data[[j]]," --maf .05 --freq --extract ",PRUNEIN[[j]]," --exclude ", srcPath,"/qc/ref/duplicatevars.txt --clump ", gwas_file, " --clump-range-border 20 --clump-range ~/Documents/FLEET/entrezhg19.gene.map --clump-p1 1.0 --clump-p2 1.0 --clump-r2 ", opt$r2," --clump-kb ",opt$`ld-window`," --clump-snp-field ", opt$`snp-field`," --clump-field ", opt$pcol," --out ",srcPath,"/qc/clumped/CLUMP_",opt$out,"_",j, sep = "")
#       system(cmd, ignore.stdout=F)
#     }
#     
#     } else{
#   
#     for(j in 1:length(ref_data)){
#       
#       cat("\r    Progress...:",j,"of",length(ref_data))
#       # ensure that spaces are escaped in GWAS file path
#       gwas_file = gsub(" ", "\\ ", opt$gwas, fixed = T)
#       cmd = paste(PathToPlink,"plink --bfile  ", ref_data[[j]]," --maf .05 --freq --exclude ", srcPath,"/qc/ref/duplicatevars.txt --clump ", gwas_file, " --clump-range-border 20 --clump-range ~/Documents/FLEET/entrezhg19.gene.map --clump-p1 1.0 --clump-p2 1.0 --clump-r2 ", opt$r2," --clump-kb ",opt$`ld-window`," --clump-snp-field ", opt$`snp-field`," --clump-field ", opt$pcol," --out ",srcPath,"/qc/clumped/CLUMP_",opt$out,"_",j, sep = "")
#       system(cmd, ignore.stdout = F)
#       
#     }
#     
#   }
  cat("\rIntervals created!")
}


## Combine LD-clumps with summary statistics
# fleetLD = function() {
#   
#   # read 1KG reference data
#   ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
#   ref_data = gsub(".bim", "", ref_data)
#   ref_data
#   
#   cat("##### Reading in clumped loci #####")
#   clumped = list.files(path = paste(srcPath,"/qc/clumped/",sep=""), full.names = T, pattern = glob2rx("CLUMP*clumped.ranges"))
#   # clumped = list.files(path = paste(srcPath,"/qc/pruned/",sep=""), full.names = T, pattern = glob2rx("*prune.in"))
#   clumped = clumped[grepl(opt$out, clumped)]
#   
#   # Read in clumps
#   clumped.df = lapply(clumped, function(x) fread(x,  stringsAsFactors = FALSE, showProgress = FALSE))
#   # clumped.df = lapply(clumped, function(x) fread(x, header = F, stringsAsFactors = FALSE, showProgress=FALSE))
#   sum(unlist(lapply(clumped.df, nrow)))
#   clumped.df = rbindlist(clumped.df)
#   clumped.df$POS = unlist(lapply(strsplit(clumped.df$POS, "[:]"), function(x) x[[2]]))
#   rangeCut = strsplit(clumped.df$POS, "[..]")
#   clumped.df$START = as.integer(unlist(lapply(rangeCut, function(x) x[[1]])))
#   clumped.df$END = as.integer(unlist(lapply(rangeCut, function(x) x[[3]])))
#   clumped.df$LOCSIZE = abs(clumped.df$END - clumped.df$START)
#   # clumped.df = clumped.df[clumped.df$N > 0, ] # intervals overlapping genes
#   
#   clumped.df = merge(clumped.df, sumstats[,c("SNP","ORc","MAFc","Vg")], by="SNP")
#   
#   
#   # gene counter
#   
#   splitter = strsplit(clumped.df$RANGES, "[,]")
#   counts = lapply(splitter, length)
#   counts[which(splitter %in% "[]")] = 0
#   clumped.df$NGENES = unlist(counts)
#   
#   freqs = list.files(path = paste(srcPath,"/qc/clumped", sep = ""), full.names = T, pattern = glob2rx("CLUMP*frq"))
#   # freqs = list.files(path = paste(srcPath,"/qc/pruned", sep = ""), full.names = T, pattern = glob2rx("*frq"))
#   freqs = freqs[grepl(opt$out, freqs)]
#   
#   freqs.df = lapply(freqs, function(x) fread(x,  showProgress=FALSE))
#   sum(unlist(lapply(freqs.df, nrow)))
#   freqs.df = rbindlist(freqs.df)
#   cols = c("SNP", "MAF")
#   freqs.df = freqs.df[,cols,with = FALSE]
#   # 
#   
#   # == average maf for interval
#   ld.maf = merge(clumped.df, freqs.df,by="SNP")
#  
# 
#   cat("....Removing LD intervals in xMHC region")
#   ld.maf = ld.maf[!ld.maf$SNP %in% rownames(mhc), ]
#   
#   ld.maf = data.table(ld.maf)
#   
#   ld.maf$LOGP = -log10(ld.maf$P)
# 
#   # ld.df$ZSCORE = my.invnorm(qnorm(1-ld.df$P))
#   ld.maf$ZSCORE = qnorm(1-ld.maf$P)
#  
#   ld.maf = ld.maf[!duplicated(ld.maf)]
#   
#   cat("##### Calculating average P-value for intervals #####")
#   clumped = list.files(path = paste(srcPath,"/qc/clumped/",sep=""), full.names = T, pattern = glob2rx("CLUMP*clumped"))
#   clumped = clumped[grepl(opt$out, clumped)]
# 
#   tagged = lapply(clumped, function(x) fread(x))
#   tagged = rbindlist(tagged)
#   tagged$SP2 = ifelse(tagged$SP2 == "NONE", tagged$SNP, tagged$SP2)
#   splitter = strsplit(tagged$SP2, ",")
#   sizesplit = unlist(lapply(splitter, function(x) length(x)))
# 
#   repNames = mapply(rep, tagged$SNP, sizesplit)
#   repNames = unlist(repNames, use.names = F)
# 
#   tag = data.frame(INDEX = repNames)
#   tag$SNP = unlist(splitter,use.names = FALSE)
#   prt = unlist(lapply(strsplit(tag$SNP, "[(]"),function(x) x[[1]]))
#   tag$SNP  = unlist(prt)
#   
#   names(sumstats)[names(sumstats) %in% c(opt$`snp-field`, opt$pcol)] = c("SNP","P")
#   sumstats = data.table(sumstats)
#   
#   tag = merge(data.table(tag), sumstats[,c("SNP","P"),with=F], by = "SNP")
#   
#   intervalMeanP = tag[,list(MeanP = mean(P)),by=(INDEX)]
#   colnames(intervalMeanP)[1] = "SNP"
#   
#   ld.maf = merge(ld.maf, intervalMeanP, by = "SNP")
#   
#   ld.maf$meanZSCORE = qnorm(1-ld.maf$MeanP)
#   
#   # remove intervals with one marker 
#   ld.maf = ld.maf[!ld.maf$N == 1]
#   
#   # remove outliers 
#   ld.maf = ld.maf[!is.infinite(ld.maf$meanZSCORE)]
#   
#   cat("\rDetected",nrow(ld.maf),"intervals for analysis")
#   cat("\rDetected",nrow(ld.maf[ld.maf$P < 5e-08]),"genome-wide significant intervals!")
#   
#   # calculate number of SNPs per kilobase of interval
#   ld.maf$LOCSIZE[ld.maf$LOCSIZE == 0] = 1
#   ld.maf$snp_density = ld.maf$N/(ld.maf$LOCSIZE/1e3)
#   
#   # plot z-score (meanP) vs z-score (minP)
#   grab = ld.maf[sample(nrow(ld.maf), 5e3)]
#   plot(grab$ZSCORE, grab$meanZSCORE,
#        col='dodgerblue', cex = 0.5, pch = 19,
#        xlab="z-score (min P)", ylab="z-score (mean P)")
#   abline(lm(grab$meanZSCORE ~ grab$ZSCORE),col='orange',lwd=2)
#   abline(v = 0, col = 'grey', lty = 2)
#   abline(h = 0, col = 'grey', lty = 2)
#   
#   
#   # Add variance attributed to SNPs (intervals), repeat across j-bootstraps
#   probsample = c(0.01, .02, .03, .04, .05, .06,.07,.08,.1,.11,.12,.125,.15)
#   
#   VGstats = list()
#   for(x in  1:length(probsample)){
#     
#   Nperm = 1e2
#   bootStrapVG = rep(NA, Nperm)
#   for(i in 1:Nperm){
#     cat("\rBootstraps:",i)
#     rows = sample(nrow(ld.maf), round(probsample[[x]]*nrow(ld.maf)) )
#     sub = clumped.df[rows]
#     bootStrapVG[[i]] = sum(sub$Vg)
#   }
#   Vg = mean(unlist(bootStrapVG))
#   bootVal = unlist(bootStrapVG)
#   SE = (((Nperm-1)/Nperm) * sum((bootVal - Vg)^2))^0.5
#   cat("\rEstimated h2g:",round(Vg,3),"(SE:",round(SE,3),")")
#   VGstats[[x]] = data.frame(Vg=Vg, SE=SE, Prop = probsample[[x]])
#   }
#   VGstats = ldply(VGstats)
#   
#   ggplot(VGstats, aes(x = factor(Prop), y = Vg*100, fill = Prop)) +
#     geom_bar(stat='identity',col='black') +
#     theme_classic()+
#     guides(fill=F)+
#     xlab("Proportion of intervals in bootstraps") +
#     ylab("Total variance (Vg)") +
#     geom_errorbar(aes(ymin = Vg*100 - SE, ymax = Vg*100 + SE), width = 0.1)
#   
#   # save interval sumstats
#   fwrite(data.table(ld.maf),
#          file = paste(srcPath,"/out/",opt$out,".sumstats",sep=""),
#          quote=F,row.names=F,
#          sep="\t",showProgress = FALSE, nThread = opt$threads)
#   
#   # make sumstats available outside function
#   gwasMunged <<- ld.maf 
#   
# }


## FLEET annotation (via bedtools)

fleetAnnot = function(x) {
  
  cat("\n##### Bedtools module ######")
  # order snp into bed
  snp.bed = gwasMunged
  colOrder = c("CHR","IntervalStart","IntervalEnd","SNP")
  snp.bed = snp.bed[,colnames(snp.bed) %in% c("SNP","CHR","IntervalStart","IntervalEnd"),with=F]
  snp.bed$CHR = paste("chr",snp.bed$CHR,sep="")
  snp.bed = snp.bed[,match(colOrder, colnames(snp.bed)),with=F]


  write.table(snp.bed,
              file = paste(srcPath,"/annotations/bedtools/snp.bed",sep=""),
              quote = F, row.names = F, col.names = F, sep = "\t")
  
  
  # === load annotations, run bedtools: 
  rdataAnnots = list.files(path = paste(srcPath,"/annotations/", sep=  ""), full.names = T, pattern = '.Rdata')
  rdataAnnots

  for( a in 1:length(rdataAnnots)){
    
    
    cat("\nLoading annotation coordinates from file:",rdataAnnots[[a]])
    
    # Convert annotation file to bedtools format:
    annot.gr = readRDS(rdataAnnots[[a]])
    mhcCoords = GRanges(seqnames="chr6", IRanges(start = 24e6, end = 35e6))
    MHCannotations = findOverlaps(annot.gr, mhcCoords)
    if(length(queryHits(MHCannotations)) > 0){ annot.gr = annot.gr[-queryHits(MHCannotations)]}
    
    names(annot.gr) = NULL # ensure rownames are NULL (this can create problems if duplicates exist, which is highly likely)
    annot.df = as.data.frame(annot.gr)
    annot.bed = annot.df[,c("seqnames", "start", "end", "id")]
    annot.bed = annot.bed[order(annot.bed$seqnames, annot.bed$start), ]
    annot.bed = data.frame(annot.bed, empty = NA, empty = "+" )
    
    # === size of matrix
    castMb = round( nrow(annot.bed) * ncol(annot.bed) * 8/1e6, 3)
    nAnnots = length(unique(annot.bed$id))
    castMb
    
    if(opt$bedtools != ""){PathToBedtools = paste(opt$bedtools,"/",sep="")} else{PathToBedtools = ""}
    
    if( castMb > 25 & nAnnots > 1 ){ 
      
      cat("\n   Annotation matrix is approximately",castMb, "Mb in size and contains",nAnnots,"terms. Writing into chunks to reduce memory strain...\n")
      
      unique.tracks = unique(annot.bed$id)
      
      size = min(20, nAnnots)
      chunkSize = ceiling(length(unique(annot.bed$id))/size)
      
      trackChunks = split(unique.tracks, ceiling(seq_along(unique.tracks)/ chunkSize ))
      
      
      # === memory efficient operation for large annotation files: 
      for ( tc in 1:length(trackChunks)){
        
        fwrite(annot.bed[annot.bed$id %in% trackChunks[[tc]], ], 
                    file = paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed",sep=""), 
                    quote= F, row.names = F, col.names = F, sep = "\t")
        
        # Run bedtools variant annotation:
        cmd = paste(PathToBedtools,"bedtools intersect -a ",srcPath,"/annotations/bedtools/snp.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed -wo > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,"_",opt$out, ".txt" , sep = "")
        system(cmd, ignore.stdout = FALSE)
        
        suppressMessages(file.remove(paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed",sep="")))
        
      }
    } else {
      
      fwrite(annot.bed, file = paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed",sep=""), quote= F, row.names = F, col.names = F, sep = "\t")
      # Run bedtools variant annotation:
      cmd = paste(PathToBedtools,"bedtools intersect -a ",srcPath,"/annotations/bedtools/snp.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed -wo > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])), "_",opt$out, ".txt" , sep = "")
      system(cmd, ignore.stdout = FALSE)
      suppressMessages(file.remove(paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed",sep="")))
      
    }
    
}

}




### Enrichment analysis
fleetTest = function(){

  search = ls()[grepl("gwasMunged", ls())]

  if(length(search) > 0) {ld.df <- gwasMunged} else{cat("\nReading summary statistics file...\n"); ld.df = as.data.frame(fread(paste(srcPath,"/out/",opt$out,".sumstats",sep=""), header = TRUE, showProgress=FALSE)) }

  
  
  coef_trait = list()
  glm_trait  = list()

  rdataAnnots = list.files(path = paste(srcPath,'/annotations/bedtools/',sep=""), full.names = TRUE, pattern = '.txt')
  rdataAnnots = rdataAnnots[grepl(paste(opt$out,".txt",sep=""), rdataAnnots)]
  rdataAnnots

  for( a in 1:length(rdataAnnots)){

    
    cat("\rLoading variant annotations:",rdataAnnots[[a]],"...")
    
    # read in and convert to wide format 
    
    bedtooled = fread(rdataAnnots[[a]], header = FALSE, stringsAsFactors = FALSE)
    bedtooled = bedtooled[,c(1,2,3,4,6,7,8,11)]
    bedtooled = bedtooled[!duplicated(bedtooled)] # remove duplicate rows 
    bedtooled$V2 = as.integer(bedtooled$V2)
    bedtooled$V3 = as.integer(bedtooled$V3)+1L
    bedtooled$V11 = as.integer(bedtooled$V11)
    bedtooled$percent_overlap = bedtooled$V11/abs(bedtooled$V6 - bedtooled$V7)
    
    cat("\rRemoving annotations that overlap < 60% of interval!")
    bedtooled = bedtooled[bedtooled$percent_overlap > 0.6] 
    bedtooled$value = 1
    
    if(nrow(bedtooled) < 1) next
    
    chr_annots = paste("chr",1:22,sep = "")
    
    cast_list = list()
    for(chr in 1:length(chr_annots)){
      
      cat("\n....Gathering annotations for chr:",chr)
      subed = bedtooled[bedtooled$V1 %in% chr_annots[[chr]]]
      if(nrow(subed) < 1) next
      cast = suppressMessages(dcast.data.table(subed, V4 ~ V8, value.var="value"))
      # setnames(setDT(cast),"V4","SNP") # rename SNP column
      colnames(cast)[1]="INDEX"
      cast[is.na(cast)] = 0
      
      counts = ifelse(cast[,2:ncol(cast)] > 0, 1,0)
      cast = data.table(cast$INDEX, counts)
      
     
      cat("\n    Number of unique annotations detected:", ncol(cast)-1)
      # tmp[is.na(tmp)] = 0
      # cast_list[[chr]] = setDT(tmp)
      cast_list[[chr]] = cast
    }
    
    cast = rbindlist(cast_list, fill = TRUE)
    
    names(ld.df)[names(ld.df) %in% "SNP"] = "INDEX"
    # === append missing snps 
    missing_snps = ld.df$INDEX[!ld.df$INDEX %in% cast$V1]
    
    empty = matrix(0, nrow = length(missing_snps), ncol = ncol(cast)-1)
    empty_df = data.frame(SNP = missing_snps, empty)
    colnames(empty_df) = colnames(cast)
    
    cat("\n    Resolving missingness")
    full_annot = rbindlist(list(setDT(empty_df), setDT(cast)))
    full_annot = as.data.frame(full_annot)
    colnames(full_annot)[1] = "INDEX"
    
    
    full_annot = full_annot[match(ld.df$INDEX, full_annot$INDEX), ]
    colnames(full_annot) = gsub(" ", "_", colnames(full_annot)) # replace spaces in ID with underscore
    colnames(full_annot) = gsub("[,]", "", colnames(full_annot)) # replace commas in ID with underscore
    
  #  ### == write a flat file containing intervals and all annotations
  #   if( a > 1){
  #     refRead = fread(file = paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.txt",sep=""), h=TRUE);
  #     refRead = data.table(refRead, full_annot)
  #     refRead = refRead[,!duplicated(colnames(refRead)),with=F]
  # 
  #   fwrite(data.table(refRead),
  #          file = paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.txt",sep=""),
  #          row.names = F, quote = F, nThread  = opt$threads, sep = "\t")} else{ fwrite(data.table(full_annot),
  #          file = paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.txt",sep=""),
  #          row.names = F, quote = F, nThread  = opt$threads, sep = "\t")}
  # 
  # }

  # # create PED file
  # readFlat = fread( paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.txt",sep=""),h=T)
  # ped = readFlat[,!colnames(readFlat) %in% "INDEX",with=F]
  # ped = ifelse(ped == 1, "I", "O")
  # 
  # colNames = colnames(ped)
  # 
  # colnames(ped) = paste("ID",1:ncol(ped),sep="")
  # orA = seq(from = 1, by = 2, length.out = ncol(ped))
  # orB = seq(from = 2, by = 2, length.out = ncol(ped))
  # order = c(orA,orB)
  # 
  # # ped2 = matrix(0, nrow = nrow(ped), ncol = ncol(ped))
  # ped2 = data.table(ped)
  # colnames(ped2) = paste("ID_B",1:ncol(ped2),sep="")
  # 
  # names(order) = c(colnames(ped), colnames(ped2))
  # order = order[order(order,decreasing=F)]
  # 
  # slap = data.table(ped2,ped)
  # setcolorder(slap, names(order))
  # 
  # # Write MAP file
  # map = data.frame(CHR = 1, ID = colnames(readFlat)[-1], POS = 0, BP = 1:(ncol(readFlat)-1))
  # fwrite(map,  file = paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.map",sep=""), quote = F, row.names = F, sep=  "\t", col.names=F)
  # # create FAM file
  # fam = data.frame(FID = 1:nrow(readFlat), IID = readFlat$INDEX, PAT = 0, MAT = 0, SEX = 0, PHEN = 0)
  # ped = data.table(fam, slap)
  # # write ped file
  # fwrite(ped,  file = paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.ped",sep=""), quote = F, row.names = F, sep=  " ", col.names=F)
  # # Write PHEN file
  # phen = data.frame(FID = 1:nrow(readFlat), IID = readFlat$INDEX, PHEN = ld.df$ZSCORE, NGENES = ld.df$NGENES, LOCSIZE = ld.df$LOCSIZE, NSNP = ld.df$N, MAF = ld.df$MAF)
  # fwrite(phen,  file = paste(srcPath,"/out/",opt$out,"_IntervalAnnotations.phen",sep=""), quote = F, row.names = F, sep=  "\t", col.names=T)


# }


    
    # full_annot[,-1] = apply(full_annot[,-1], 2, function(x) my.invnorm(x))
    
    # === merge data table with ld.df stats
    
    grab = ld.df
    full_annot = data.table(full_annot)
    full_annot = full_annot[match(grab$INDEX, full_annot$INDEX)]
    tmp = data.table(grab, full_annot)
    tmp = tmp[,!duplicated(colnames(tmp)),with=F]
    tmp = as.data.frame(tmp)
    tmp[is.na(tmp)] = 0
    
    # === fix format of annotation labels 
    all_df = tmp
  
    # 
    
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
    
    
    tracks = colnames(all_df)[!colnames(all_df) %in% colnames(ld.df)]
    
    
    if(length(tracks) < 1) next
    
    # === FLEET enrichment analysis via linear regression
    fits = list()
    rsq.fits = list()
    
    
    SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
    SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
    
    all_df = all_df[!is.infinite(all_df$ZSCORE), ]
    # all_df = all_df[!duplicated(all_df), ]
    
    cat("\nTotal number of LD-intervals available for enrichment tests:", nrow(all_df),"\n")
    
    # all_df = all_df[all_df$ZSCORE < 4,]
    # all_df$ZSINV = my.invnorm(all_df$ZSCORE)
    
    # opt$speed = "fast"
    ## === linear model in fast mode
    if(opt$speed == 'fast'){
      
      dt_all_df = data.table(all_df)
      
      misCol = colSums(dt_all_df[,colnames(dt_all_df) %in% tracks,with=F] == 1)

      misCol = misCol[misCol < opt$`annot-cnt`]

      if(length(misCol) > 0){
        dt_all_df = dt_all_df[,!colnames(dt_all_df) %in% names(misCol),with=F]
        }
      
      keep_tracks = tracks[tracks %in% colnames(dt_all_df)]
      
      response_matrix = as.matrix(dt_all_df[,colnames(dt_all_df) %in% keep_tracks,with=F])
      response = dt_all_df$MeanZSCORE
      
      
      cat("\rFitting linear regressions in fast mode (Zscore ~ Set (multiple) + covariates)!\n")
      start_reg = proc.time()
      mod0 = model.matrix( ~. , data = dt_all_df[,colnames(dt_all_df) %in% c(keep_tracks, "MAFc", "LDSCORE", "SNPDENSITY"),with=F])
      lmFast = lm(response ~ -1 + mod0)
      cat("\r...Time elapsed for regression:", (proc.time() - start_reg)[3], 'seconds\n')
      lmFastCoef = summary(lmFast)
      
     
          lmCoef = broom::tidy(summary(lmFast)$coefficients)
          lmCoef = lmCoef[!grepl("Intercept|MAFc|LDSCORE|SNPDENSITY", lmCoef$.rownames), ]
          lmCoef$.rownames = gsub("mod0", "", lmCoef$.rownames)
         
      
      lmDf = lmCoef
      
      lmDf = lmDf[,!colnames(lmDf) %in% ".id"]
      lmDf = lmDf[,!colnames(lmDf) %in% "Term"]
      colnames(lmDf) = c("SET_NAME",  "Beta", "SE", "T","P")
      
      lmDf$SET_SOURCE = SETNAME
      
      gwVars = gw.sig.vars[,colnames(gw.sig.vars) %in% c("INDEX",keep_tracks)]
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
        splitter = strsplit(gwVars_agg$GWS_VAR_HITS, "[|]")
        lmDf$N = unlist(lapply(splitter,length))
      } else {lmDf$GWS_VAR_HITS = "*"}
      
      # == Calculate Vg for SNPs in annotation category
      
      VgCalc_list = list()
      for(x in 1:length(keep_tracks)){
        cat("\rCalculating Vg for SNPs in annotation category:",round(x/length(keep_tracks)*100,2),"%")
      sumVgInAnnot = which(dt_all_df[,colnames(dt_all_df) %in% keep_tracks[[x]],with=F] == 1)
      
      zall = log(dt_all_df$ORc)/dt_all_df$se
      zall = zall[sumVgInAnnot]
      
      uniqueZval = length(unique(zall))
      
      if(uniqueZval > 20){
        VgInAnnot = try(SumVg.binary(zall=zall, 
                                     method="jack",
                                     d=length(zall)/5,
                                     repl=5, 
                                     out="unconditional",
                                     caseNo=35476, 
                                     ctrlNo=46839, 
                                     K=0.01))
        if(VgInAnnot != "error"){
          VgCalc = data.frame(t(unlist(VgInAnnot)))
          VgCalc$VgZ = VgCalc[,1]/VgCalc[,2]
          VgCalc$VgP = pnorm(-abs(VgCalc$VgZ))
          VgCalc$VgProportion = VgCalc$Est.SumVg * (length(sumVgInAnnot)/nrow(dt_all_df))
          VgCalc$N = length(sumVgInAnnot)
          VgCalc$SET_NAME = keep_tracks[[x]]
          VgCalc_list[[x]] = VgCalc
          } 
        
      
      }
      
      }
      VgCalc_df = ldply(VgCalc_list)
      VgCalc_df = VgCalc_df[match(lmDf$SET_NAME, VgCalc_df$SET_NAME), ]
      
      lmDf = data.frame(lmDf, VgCalc_df[,!colnames(VgCalc_df) %in% "SET_NAME"])
      
      coef_df = lmDf
      }
    
    
    # opt$speed = "slow"
    # === linear model in slow mode
    if(opt$speed == "slow"){
      
      
      dt = setDT(all_df)
      
      dt_all_df = data.table(all_df)
      
      misCol = colSums(dt_all_df[,colnames(dt_all_df) %in% tracks,with=F] == 1)
      
      misCol = misCol[misCol < opt$`annot-cnt`]
      
      if(length(misCol) > 0){
        dt_all_df = dt_all_df[,!colnames(dt_all_df) %in% names(misCol),with=F]
      }
      
      phenotypes = dt_all_df[,colnames(dt_all_df) %in% tracks,with=F]
      phenotypes = data.frame(phenotypes)
      ZSCORE = dt_all_df$MeanZSCORE
      MAF = dt_all_df$MAFc
      LDSCORE = dt_all_df$LDSCORE
      SNPD = dt_all_df$SNPDENSITY
      
      
      cat("\rFitting linear regressions (Zscore ~ Set + covariates)\n")
      startclock = proc.time()
      lmFit = list()
      
      tracks = tracks[tracks %in% colnames(dt_all_df)]
      
      lmFit = list()
      for(x  in 1:length(tracks)){
     
         cat("\rRegressions:",round(x/length(tracks),3)*100,"%")
        
      if( length(tracks) > 0) {
      mod0 = model.matrix( ~ phenotypes[,x] + MAF + LDSCORE + SNPD) } else{
      mod0 = model.matrix( ~ phenotypes + MAF + LDSCORE + SNPD) }
      
      fit = broom::tidy(summary(lm(ZSCORE ~ -1 + mod0))$coefficients)
      rm(mod0)
      fit$SET_NAME = tracks[[x]]
      fit = fit[grepl("phenotype", fit$.rownames), ]
      
      sumVgInAnnot = which(dt_all_df[,colnames(dt_all_df) %in% tracks[[x]],with=F] == 1)
      
      zall = log(dt_all_df$ORc)/dt_all_df$se
      zall = zall[sumVgInAnnot]
      
      uniqueZval = length(unique(zall))
      
      if(uniqueZval > 20){
      VgInAnnot = try(SumVg.binary(zall=zall, 
                               method="jack",
                               d=length(zall)/5,
                               repl=5, 
                               out="unconditional",
                               caseNo=35476, 
                               ctrlNo=46839, 
                               K=0.01))
      if(VgIntAnnot != "error"){
      VgCalc = data.frame(t(unlist(VgInAnnot)))
      VgCalc$VgZ = VgCalc[,1]/VgCalc[,2]
      VgCalc$VgP = pnorm(-abs(VgCalc$VgZ))
      VgCalc$VgPercent = VgCalc$Est.SumVg/(length(sumVgInAnnot)/nrow(dt_all_df)*100)
      VgCalc$N = length(sumVgInAnnot)
      
      fit = data.frame(SET_NAME = fit$SET_NAME, 
                        Beta = fit$Estimate,
                        SE = fit$Std..Error,
                        T = fit$t.value,
                        P = fit$Pr...t.., VgCalc)}
      } 
      
     if(uniqueZval <= 20){ fit = data.frame(SET_NAME = fit$SET_NAME, 
                                  Beta = fit$Estimate,
                                  SE = fit$Std..Error,
                                  T = fit$t.value,
                                  P = fit$Pr...t.., N = length(sumVgInAnnot))}
      
      
      lmFit[[x]] = fit
      }
      
      tidy = ldply(lmFit)
      
      set_names = tidy$SET_NAME
      bsStats = list()
      for(x in 1:length(set_names)){
        
        gwas_vars = gw.sig.vars[,c(tracks[[x]], "INDEX", "P")]
        gwas_vars = ifelse(gwas_vars[,1] == 1, as.character(gwas_vars$INDEX), "")
        gwas_vars = gwas_vars[!gwas_vars %in% ""]
        gwas_vars = paste(gwas_vars, collapse = "|", sep = "")
        
        SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
        SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
       
        tidy$GWSVARNAMES[[x]] = ifelse( gwas_vars == "", "*", gwas_vars)
        
        # === binomial enrichment analysis
        
        expectedProp = table(phenotypes[,x])[[2]]
        nSize = nrow(phenotypes)
        baseline = expectedProp/nSize
        
        if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(5e-08, 1e-06, 1e-05, 1e-04, 5e-02, 1e-1, 2e-1, 5e-1)}
        
        cat("\rBinomial tests:",round(x/length(set_names),3)*100,"%");flush.console()
        binomstats = list()
        for( z in 1:length(pthreshold)){
          
          gws.hits.1 = dt[dt$P < pthreshold[[z]]]
          hits1 = colSums(gws.hits.1[,tracks[[x]],with=FALSE])[[1]]
          
          observedProp = hits1/nrow(gws.hits.1)
          binomPval = binom.test(x = hits1, n =  nrow(gws.hits.1), p = baseline, "greater")$p.value
          bfse = data.frame(FSE = observedProp/baseline, PVAL = binomPval)
          colnames(bfse) = paste(paste("BINOM_",pthreshold[[z]],sep=""),"_",colnames(bfse),sep="")
          binomstats[[z]] = bfse
        }
        binomstats_df = as.data.frame(do.call(cbind, binomstats))
        binomstats_df$SET_NAME = set_names[[x]]
        bsStats[[x]] = binomstats_df
      }
      bsStats = ldply(bsStats)
      
        cat("Time elapsed:",(proc.time() - startclock)[[3]],"seconds")
        
        tidy = merge(tidy, bsStats, by= "SET_NAME")
        # Combine model results
        coef_df = tidy
  
    }
    
    
    
    # Permutation enrichment analysis
    all_df = as.data.frame(all_df)
    
    if(opt$`fleet-permutation` == TRUE){
      
      
      if(opt$`robust-permutation` == TRUE){
        
        cat("\nRunning permutation analysis in robust mode!\n")
        
        ## Run permutation (GW-SIG):
        # pval.check = which(coef_df$`Pr(>|t|)` <= opt$permStartThreshold & coef_df$Estimate > 0)
        
        dt = setDT(all_df)
       
        gc()
        
        # subset SNPs by significance thresold
        if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(1e-06, 1e-04, 1e-03)}
        
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
            sub = sub[!duplicated(sub), ]
            gws.hits.1 = gws.hits.1[match(sub$INDEX, gws.hits.1$INDEX)]
            sub = data.table(sub, gws.hits.1[,c("INDEX", "SNPDENSITY"),with=FALSE])
            sub = sub[!duplicated(sub)]
            
            # === Calculate snp density per interval
            snpdensity = sub[,list(meanSD = mean(SNPDENSITY), varSD = sd(SNPDENSITY), lengthSD = length(SNPDENSITY)),by=list(variable)]
            snpdensity$varSD[is.na(snpdensity$varSD)] = 0
            snpdensity$lb95 =  snpdensity$meanSD - (snpdensity$varSD/sqrt(snpdensity$lengthSD))*2.58
            snpdensity$ub95 =  snpdensity$meanSD + (snpdensity$varSD/sqrt(snpdensity$lengthSD))*2.58
            
            permPval = list()
            for(v in 1:nrow(snpdensity)){
              
              # permutable = dt[!dt$INDEX %in% sub$INDEX & dt$snp_density >= snpdensity$lb95[[v]] & dt$MAF <= snpdensity$ub95[[v]], c("SNPDENSITY", as.character(snpdensity$variable[[v]])),with=FALSE]
              permutable = dt[,c("SNPDENSITY", as.character(snpdensity$variable[[v]])),with=FALSE]
              
              total_size = nrow(gws.hits.1)
              hit_ratio = (snpdensity$lengthSD[[v]]/total_size)[[1]]
              
              if(nrow(permutable) > total_size){ 
                
                permSum = list(rep(NA, nPerms))
                for(p in 1:nPerms){
                  tmp = permutable[sample(nrow(permutable), total_size)]
                  permSum[[p]] = sum(tmp[,2])/nrow(tmp)
                  rm(tmp) # remove vector
                }
                cat("\r   Permutation analysis with threshold P-value < ", pthreshold[[pt]] ,":",round(v/nrow(snpdensity)*100, 2),"% complete")
                pval = sum(unlist(permSum) >= hit_ratio)/nPerms} else {pval = NA}
                  
                  
                pval = ifelse(pval == 0, 1/nPerms, pval)
              
              permPval[[v]] = pval
            }
              
             
            }
            
            colHeaderA = paste("Emp_P",pt,sep="")
            permDf = data.frame(SET_NAME = snpdensity$variable, Pval = unlist(permPval))
            colnames(permDf)[2] = colHeaderA
            
            cat("\n")
            
            # ptPermDf[[pt]] = permDf
            
            coef_df = merge(coef_df, permDf, by = "SET_NAME", all.x = TRUE)
            
          }
          
        }
        
        
        
        
      } # close permutation
      
    
    cat("\nWriting results...\n")
    
    if(a == 1) { write.table(coef_df, file = paste(srcPath,"/out/",opt$out,"_fleetEnrichment.txt",sep=""), quote = F, row.names = F, sep = "\t")} else{ write.table(coef_df, file = paste(srcPath,"/out/",opt$out,"_fleetEnrichment.txt",sep=""), quote = F, row.names = F, col.names = F, sep = "\t", append = TRUE)}

  } # close for loop (run next annotation set)
  
} # close fleetTest()


### Plotting mechanisms:

fleetPlot = function(x) {
  
  cat("\n#####Creating a QQ-plot######\n")
  par(mfrow=c(1,1))
  
  coef_df = read.table(paste(srcPath,"/out/",opt$out,"_fleetEnrichment.txt",sep=""), header = T, fill = NA)
  
  coef_df  = coef_df[coef_df$Estimate > 0 & !is.na(coef_df$Estimate), ] # enrichment QQplot
  
  observed <- sort(coef_df$Pr...t..)
  lobs <- -(log10(as.numeric(observed)))
  lobs[lobs > 12] <- 12
  
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
  lines(c(0,12), c(0,12), col="black", lwd=1.5)
  
  points(x = lexp, y = lobs, pch = 8, cex = .5, col = "black")
  
  text(x = 2, y = 0.5, cex = 1, labels = substitute(paste("Median ", chi^2, " = ", MC), list(MC=format(medianchi, digits = 3))))
  text(x = 2, y = .25, cex = 1, labels = substitute(paste(lambda, " = ", LM), list(LM=format(lambda,digits = 4))))
  
  ####### end qq plot ########
}




### Run operations

total_clock_start = proc.time()

readSNPcoords() # always run (from 1000G)

xMHCregion() # always run

fleetReadGWAS()

# readAnnots() # always run

# if(opt$`fleet-prune-ref` == TRUE){
  # refPrune()
# }


if(opt$`fleet-annotate` == TRUE){
  refPrune()
  fleetClump()
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
