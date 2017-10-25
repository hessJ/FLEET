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
|| (C) 2017 Jonathan L. Hess, PhD and Stephen J. Glatt, PhD
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
  
  make_option(c("-R", "--r2"), type="double", default=0.6, 
              help="R-squared threshold for linkage disequilibrium calculations [default = %default]", metavar="double"),
  
  make_option(c("-W", "--ld-window"), type="integer", default=1000, 
              help="Size of window (kilobases) for calculating linkage disequilibrium [default = %default]", metavar="integer"),
  
  make_option(c("-S", "--snp-field"), type="character", default="", 
              help="SNP column header in GWAS file", metavar="character"),
  
  make_option(c("-P", "--pcol"), type="character", default="", 
              help="P-value column header in GWAS file", metavar="character"),
  
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
  
  make_option(c("--SumVg"), type="logical", default=FALSE, 
              help="Estimate heritability based on SNP z-scores in annotation category with SumVg method [default = %default]", metavar="logical"),
  
  make_option(c("--binomial"), type="logical", default=FALSE, 
              help="Perform binomial tests of enrichment across p-value bins [default = %default]", metavar="logical"),
  
  ## Run options....
  
  make_option(c("--Nca"), type="numeric", default = "",
              help="Number of cases in GWAS.  Required for estimation of Vg for SNPs in category [default = %default]", metavar="numeric"),
  
  make_option(c("--Nco"), type="numeric", default = "",
              help="Number of controls in GWAS. Required for estimation of Vg for SNPs in category [default = %default]", metavar="numeric"),
  
  make_option("--refa1", type="character", default=TRUE, 
              help="Header for column containing reference allele", metavar="character"),
  
  make_option("--popprev", type="numeric", default="", 
              help="Population prevalence for phenotype (assumed as binary)", metavar="numeric"),
  
  make_option(c("-M", "--fleet-prune-ref"), type="logical", default = TRUE,
              help="Initiate LD-pruning step of 1KG reference data. Only needs to be run once. [default = %default]", metavar="logical"),
  
  make_option(c("-A", "--fleet-annotate"), type="logical", default=TRUE, 
              help="Annotate LD-clumps with bedtools [default = %default]", metavar="logical"),
  
  make_option(c("-E", "--fleet-enrichment"), type="logical", default=TRUE, 
              help="Perform enrichment analysis with weighted linear models [default = %default]", metavar="logical"),
  
  make_option("--fleet-permutation", type="logical", default=FALSE, 
              help="Perform enrichment analysis with permutation (randomizing annotations) [default = %default]", metavar="logical"),
  
  # make_option("--fast-permutation", type="logical", default=FALSE, 
  #             help="Simple permutation analysis [default = %default]", metavar="logical"),
  
  make_option("--robust-permutation", type="logical", default=TRUE, 
              help="Permutation analysis that will sample variants from the MAF bin of target SNPs [default = %default]", metavar="logical"),
  
  make_option("--speed", type="character", default='fast', 
              help="Change behavior of linear models (fast mode; multiple sets added as predictors to regression, slow mode; each set regressed onto Z-scores separately) [default = %default]", metavar="character"),
  
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
  # print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

if (opt$`snp-field` == ""){
  # print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

if (opt$pcol == ""){
  # print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}


if (opt$refa1 == ""){
  # print_help(opt_parser)
  stop("Must provide header for reference allele\n", call.=FALSE)
}


if (opt$popprev == ""){
  # print_help(opt_parser)
  stop("Must provide population prevalence for phenotype\n", call.=FALSE)
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
opt$Nca = 35476
opt$Nco = 46839
opt$refa1 = "a1"
opt$popprev=0.01

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
    theme(text = element_text(size = 10)) +
    geom_line(aes(x = expected, y = clower), colour ='grey' ,lwd = 0.75) +
    geom_line(aes(x = expected, y = cupper), colour = 'grey', lwd = 0.75) +
    geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), fill="grey", alpha="0.2") 
  
  # png(file=paste(srcPath,"/plots/QQplot.png",sep=""),res=300,units="in",height=7,width=7)
  print(qplot)
  # dev.off()
  
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
# xMHCregion = function(){
#   mhc =  snp137_df[snp137_df$seqnames %in% "chr6" & snp137_df$start >= 24e6 & snp137_df$start <= 35e6, ]
#   snp137_df_nomhc <<- snp137_df[!snp137_df$SNP %in% mhc$SNP, ]
#   mhc <<- mhc
# }

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
  # rename headers (SNP, OR, SE, A1)
  names(read.gwas)[names(read.gwas) %in% c("or","OR","odds","Odd")] = "OR"
  names(read.gwas)[names(read.gwas) %in% c("se","SE","StdErr","SDERR","STDERROR")] = "se"
  names(read.gwas)[names(read.gwas) %in% c(opt$`snp-field`, opt$pcol)] = c("SNP","P")
  names(read.gwas)[names(read.gwas) %in% c(opt$`refa1`)] = c("A1")
  
  # retain only markers found in reference panel
  read.gwas = read.gwas[read.gwas$SNP %in% freq$SNP]
  
  freq = freq[match(read.gwas$SNP, freq$SNP)]
  match0 = freq$A1 == read.gwas$A1
  match1 = freq$A2 == read.gwas$A1
  
  co = data.frame(match0,match1)
  bad = which(rowSums(co) == 0)
  
  cat("\nRemoving",length(bad),"variants with mis-matched alleles between reference panel and GWAS\n")
  if(length(bad) > 0){
    freq = freq[-bad];
    read.gwas = read.gwas[-bad]
  }
  
  # match alleles
  match2 = freq$A1 == read.gwas$A1
  
  freq$A1c = ifelse(match2 == F, freq$A2, freq$A1)
  freq$MAFc = ifelse(match2 == F, 1-freq$MAF, freq$MAF)
  
  read.gwas$MAFc = freq$MAFc
  
  infoCol = which(grepl("info|INFO|imputation", colnames(read.gwas)))
  
  if(length(infoCol) > 0){
    cat("\nDetected imputation qualtiy column as:",colnames(read.gwas)[infoCol])
    preinfo = nrow(read.gwas)
    postinfo = read.gwas[,colnames(read.gwas) %in% colnames(read.gwas)[infoCol],with=F]
    keepinfo = which(postinfo > .6)
    gwas.info = read.gwas[keepinfo]
    cat("\nRemoved",abs(nrow(gwas.info)-preinfo),"variants with imputation quality < 0.6\n")
  }
  
  # remove SNPs with chi-square > 30
  chisq = (log(gwas.info$OR)/gwas.info$se)**2
  outlier = which(chisq > 300)
  
  if(length(outlier) > 0){
    cat("\rRemoved",length(outlier),"variants with chi-square > 30")
    gwas.info = gwas.info[-outlier]
  }
  
  # set population prevalence
  K <- opt$popprev
  RR1 = gwas.info$OR
  RR2 = gwas.info$OR^2
  
  PA = gwas.info$MAFc
  Paa = (1-PA)^2
  PAa = 2*PA*(1-PA)
  PAA = PA^2
  muaa=0
  faa= K/(Paa + PAa*RR1 + PAA*RR2)
  fAa= RR1*faa
  fAA= RR2*faa 
  Tr = qnorm(faa,lower.tail = FALSE) 
  muAa = Tr-qnorm(fAa, lower.tail = FALSE)
  muAA = Tr-qnorm(fAA, lower.tail = FALSE)
  mean.all= PAa*muAa+ PAA*muAA
  Vg = Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
  actual.Vg =  Vg/(1+Vg) 
  
  N = opt$Nca + opt$Nco
  chisq  = (log(gwas.info$OR)/gwas.info$se)**2
  VgChisq = chisq/(N-2+chisq)
  
  gwas.info$Vg = actual.Vg
  gwas.info$VgChisq = VgChisq
  
  # add position info for SNPs
  gwas.info = gwas.info[,!grepl("chr|CHR|chrom|chromosome", colnames(gwas.info)),with=F]
  gwas.info = gwas.info[,!grepl("base|BP|bp", colnames(gwas.info)),with=F]
  snp137_df = data.table(snp137_df)
  snp137_df = snp137_df[match(gwas.info$SNP, snp137_df$SNP)]
  gwas.info$CHR = snp137_df$seqnames
  gwas.info$BP = snp137_df$start
  
  fwrite(data.table(gwas.info$SNP),
         file = paste(srcPath,"/qc/gwas.snplist",sep=""),
         quote = F, row.names = F, col.names=F, sep = "\t")
  
  sumstats <<- data.frame(gwas.info)
  
}


# = read gene coordinates:
snp2gene = function(x){
  
  cat("\nCreating Granges objects for SNP and gene coordinates")
  path_to_entrez_genes = fread(paste(srcPath,"/entrezhg19.gene.map",sep=""),h=FALSE)
  generanges = GRanges(seqnames = paste("chr", path_to_entrez_genes$V1, sep = ""), IRanges(start = path_to_entrez_genes$V2, end =path_to_entrez_genes$V3), symbol = path_to_entrez_genes$V4)
  
  # = expand gene window
  # start(generanges) = start(generanges) - 100e3
  # end(generanges) = end(generanges) + 100e3
  
  snpranges = GRanges(seqnames = sumstats$CHR, 
                      IRanges(start = sumstats$BP, end = sumstats$BP), SNP = sumstats$SNP, MAF = sumstats$MAFc, OR = sumstats$OR, SE = sumstats$se, Pval = sumstats$P)
  
  cat("\nMapping SNPs to genes (hg19)")
  # == find overlaps
  fo = findOverlaps(snpranges, generanges)
  snp2gene = snpranges[queryHits(fo)]
  gene2snp = generanges[subjectHits(fo)]
  snp2gene$symbol = gene2snp$symbol
  
  # == calculate gene scores from SNP z-scores
  geneScores = data.frame(snp2gene)
  geneScores = data.table(geneScores)
  geneScores$Effect = log(geneScores$OR)/geneScores$SE
  geneScores$Zscore = qnorm(geneScores$Pval, lower.tail = FALSE)
  
  genecoord = data.frame(generanges)
  
  width = width(generanges)
  width = data.table(symbol = generanges$symbol, width)
  
  cat("\nCalculating mean SNP p-value per gene")
  
  # minP(SNP)
  minP = geneScores[order(geneScores$Pval,decreasing = F)]
  minP = minP[!duplicated(minP$symbol)]
  minP = minP[,colnames(minP) %in% c("SNP","symbol","OR","SE"),with=F]
  
  # Mean P-value 
  meanP = geneScores[,list(Pval = mean(-log10(Pval)), MAF = mean(MAF)),by=list(symbol)]
  
  meanP$Pval = 10^-(meanP$Pval)
  meanP$meanZscore = qnorm(meanP$Pval, lower.tail = F)
  
  meanP = merge(meanP, minP, by='symbol')
  meanP = meanP[!meanP$symbol %in% ""]
  
  meanP$OR= ifelse(meanP$OR < 1, 1/meanP$OR, meanP$OR)
  meanP$MAF = ifelse(meanP$MAF, 1-meanP$MAF, meanP$MAF)
  
  # sumvg on genes
  K <- opt$popprev
  RR1 = meanP$OR
  RR2 = meanP$OR^2
  PA = meanP$MAF
  Paa = (1-PA)^2
  PAa = 2*PA*(1-PA)
  PAA = PA^2
  muaa=0
  faa= K/(Paa + PAa*RR1 + PAA*RR2)
  fAa= RR1*faa
  fAA= RR2*faa 
  Tr = qnorm(faa,lower.tail = FALSE) 
  muAa = Tr-qnorm(fAa, lower.tail = FALSE)
  muAA = Tr-qnorm(fAA, lower.tail = FALSE)
  mean.all= PAa*muAa+ PAA*muAA
  Vg = Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
  actual.Vg =  Vg/(1+Vg) 
  
  # 
  meanP$Est.Vg = actual.Vg

  # remove rows with same SNPs
  nonDup = meanP[!duplicated(meanP$SNP)]
  
  # clump top SNP 
  fwrite(data.table(nonDup$SNP), 
         file = paste(srcPath,"/qc/geneindex.snplist",sep=""), quote=F, row.names = F,sep="\t",col.names = F)
  
  
  refData = list.files(path=paste(srcPath,"/qc/ref",sep=""), pattern='.bed',full.names = T)
  refData = gsub(".bed", "", refData)
  
  if(opt$plink != ""){PathToPlink = opt$plink}else{PathToPlink = ""}
  for(z in 1:length(refData)){
    cat("\rLD-clumping gene results:",z)
    cmd = paste(PathToPlink,"/plink --bfile ",refData[[z]], " --extract ",srcPath,"/qc/geneindex.snplist --indep-pairwise 100 5 0.1 --out ",srcPath,"/qc/pruned/GENEPRUNED_",z,sep="")
    system(cmd, ignore.stdout = TRUE)
  }
  # read LD-pruned snp list for genes
  geneLdPruned = list.files(path = paste(srcPath,"/qc/pruned/",sep=""), pattern="GENEPRUNED",full.names=TRUE)
  geneLdPruned = geneLdPruned[grepl(".prune.in", geneLdPruned)]
  geneLdPruned = lapply(geneLdPruned, function(x) fread(x, h=FALSE))
  geneLdPruned = rbindlist(geneLdPruned)
  
  keepLd = meanP[meanP$SNP %in% geneLdPruned$V1]
  keepLd = keepLd[order(keepLd$OR,decreasing=TRUE)]
  keepLd = keepLd[!duplicated(keepLd$SNP)]
  
  # == Calculate median Chi-sq
  zall = log(keepLd$OR)/keepLd$SE
  
  # == estimate heritability explained by intervals (SumVg method)
  
  EstVgAll = SumVg.binary(zall=zall, 
                          method="jack",
                          d=length(zall)/5,
                          repl=5, 
                          out="unconditional",
                          caseNo=opt$Nca, 
                          ctrlNo=opt$Nco, 
                          K=opt$popprev)
  # == Total heritability:
  EstVgAll = data.frame(t(unlist(EstVgAll)))
  EstVgAll$VgZ = EstVgAll[,1]/EstVgAll[,2]
  EstVgAll$VgP = pnorm(-abs(EstVgAll$VgZ))
  EstVgAll$VgP = ifelse(EstVgAll$VgP == 0, 1e-300, EstVgAll$VgP)
  if(EstVgAll[,1] > 1.0){cat("\rWarning! Heritability estimate (liability scale) from intervals was greater than 1, possibly due to insufficient LD pruning")}
  cat("\nTotal heritability from genes (liability scale):",round(EstVgAll[,1],4),"(SE:",round(EstVgAll$SE.SumVg,4),")")
  cat("\nHeritabiliy z-score: ",EstVgAll$VgZ)
  cat("\nHeritabiliy p-value: ",EstVgAll$VgP)
  
  EstVgAll$Ngene = nrow(keepLd)

  fwrite(EstVgAll, 
         file = paste(srcPath,"/out/",opt$out,"_geneH2-sumvg.txt",sep=""),
         quote = F, row.names = F, sep="\t")
  
  # QQ plot of gene scores
  png(paste(srcPath,"/plots/qqplot.png",sep=""),res=300,units="in",height=6,width=6)
  qqPlot(meanP$Pval)
  dev.off()
  
  cat("\nCreating Manhattan plot of gene p-values")
  man = merge(meanP, data.table(genecoord), by ='symbol')
  sub = man[!man$symbol %in% ""]
  names(sub)[names(sub) %in% c("seqnames","start")]=c("CHR","BP")
  sub$CHR = as.integer(gsub("chr", "", sub$CHR))
  sub = sub[order(sub$CHR, sub$BP), ]
  
  
  col = data.frame(CHR = 1:22, col = c("navy", "forestgreen"))
  
  sub = merge(sub, col, by = "CHR")
  sub$pos = NA
  sub$pos = ifelse(sub$CHR == 1, sub$BP, NA)
  
  median_pos = list()
  chr_grab = 1:22
  for( i in 1:length(chr_grab)){
    
    if(i == 1){
      prior_max = min(sub$pos[sub$CHR == 1]) - 1
    }
    
    if(i > 1){
      k = i - 1
      prior_max = max(sub$pos[sub$CHR %in% k])
    }
    
    chr_seq = sub$BP[sub$CHR %in% chr_grab[[i]]]
    true_min = min(chr_seq)
    diff = chr_seq - true_min
    new_chr_seq = diff + (prior_max + 1)
    
    sub$pos[sub$CHR == i] <- new_chr_seq
    
    median_pos[[i]] =  mean(sub$pos[sub$CHR == i])
  }
  
  names(median_pos) = 1:22
  
  median_pos = ldply(median_pos)
  
  sub$SYMBOL = ifelse(sub$Pval < .05/nrow(sub), as.character(sub$symbol), NA)
  
  mPlot = ggplot(sub, aes(x = sub$pos, y=-log10(sub$Pval))) + 
    geom_point(size = 1.1, col  = as.character(sub$col)) + 
    xlab("Genomic coordinate") + 
    theme_classic() +
    ylab(expression(paste("-log"[10],"(P-value)"))) +
    ylim(min = 0, max = 9) +
    theme(axis.text.x = element_text(size = 5)) +
    scale_x_continuous(name="Genomic coordinate", breaks=median_pos$V1, labels=median_pos$.id) +
    geom_hline(yintercept = c(-log10(5e-08), -log10(.05/nrow(sub))), col = c("red","orange"), lwd = 0.5, linetype = c("solid", "dashed")) +
    geom_hline(aes(yintercept = 0.0), lwd = 0.5, col ='grey') +
    geom_text_repel(aes(label = sub$SYMBOL), nudge_y = 0.3, segment.size = 0.2, segment.colour = "black", fontface = 'italic', size = 2, col = "black")
  
  png(paste(srcPath,"/plots/manhattan.png",sep=""),res=300,units="in",height=6.6,width=12.5)
  print(mPlot)
  dev.off()
  
  
  # SNP density 
  geneSnpDensity = geneScores[,list(NSNP = length(SNP)),by=list(symbol)]
  geneSnpDensity = merge(geneSnpDensity, width, by="symbol")
  geneSnpDensity = geneSnpDensity[!geneSnpDensity$symbol %in% ""]


  # gene scores
  finalGeneScores = merge(meanP, geneSnpDensity,by='symbol')
  finalGeneScores$meanZscore[is.infinite(finalGeneScores$meanZscore)] = qnorm(1e-300,lower.tail=FALSE)
  finalGeneScores = merge(finalGeneScores, genecoord[,!colnames(genecoord) %in% 'width'], by='symbol')
  names(finalGeneScores)[names(finalGeneScores) %in% c("seqnames")] = "CHR"
  finalGeneScores$SNPDENSITY = finalGeneScores$NSNP/(finalGeneScores$width/1e3)
  finalGeneScores = finalGeneScores[,!colnames(finalGeneScores) %in% 'strand',with=F]
  
  fwrite(finalGeneScores, 
         file = paste(srcPath,"/out/",opt$out,".fleet.geneScores.txt",sep=""),
         quote = F, row.names = F, sep=  "\t")
  
  geneScores <<- finalGeneScores
  
  cat("\nCompleted derivation of gene scores")  
}

### Prune 1KG reference 
refPrune = function(){

  ref_data = list.files(path = paste(srcPath,'/qc/ref/',sep=""), full.names = T, pattern = 'bim')
  ref_data = gsub(".bim", "", ref_data)

  if(opt$plink != ""){PathToPlink = paste(opt$plink,sep="/")} else{PathToPlink = ""}
  
  for( i in 1:length(ref_data)){
    cat("\rPruning progress:",round(i/length(ref_data)*100,3),"%")
    indep.cmd = paste('--indep-pairwise 100 1 0.1')
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
  ld.int$DIST = abs(ld.int$BP_A - ld.int$BP_B)
  
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
  IntervalPval$MeanZSCORE = qnorm(IntervalPval$meanPval, lower.tail = FALSE)
  
  # === Combine Interval stats with sumstats 
  int.stats = data.table(IntervalRange, LDFRIENDS,LDSCORE, IntervalPval, indexPos)
  int.stats = int.stats[,!duplicated(colnames(int.stats)),with=F]
  colnames(int.stats)[1] = "SNP"
  int.stats$LOCSIZE[int.stats$LOCSIZE == 0] = 1
  int.stats$SNPDENSITY = int.stats$NSNP/(int.stats$LOCSIZE/1e3)
  
  ld.int = merge(data.table(sumstats), int.stats, by= "SNP")
  ld.int$ZSCORE = qnorm(ld.int$P,lower.tail = FALSE)
  
  ld.int = ld.int[!ld.int$NSNP == 1]
  
  cat("\rDetected",nrow(ld.int),"intervals for analysis\n")
  cat("\rDetected",nrow(ld.int[ld.int$P < 5e-08]),"genome-wide significant intervals! Some may have been removed by pruning\n")
  
  
  # save interval sumstats
  fwrite(data.table(ld.int),
         file = paste(srcPath,"/out/",opt$out,".sumstats",sep=""),
         quote=F,row.names=F,
         sep="\t",showProgress = FALSE, nThread = opt$threads)
  
  # make sumstats available outside function
  gwasMunged <<- ld.int 
  
  cat("\nIntervals created!")
}



## === Annotate intervals with bedtools2; genomic coordinates required for intervals and (epi)genomic annotation(s)
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
  
  
  gene.bed = geneScores
  colOrder = c("CHR","start", "end","symbol")
  gene.bed = gene.bed[,colnames(gene.bed) %in% colOrder,with=F]
  gene.bed = gene.bed[,match(colOrder, colnames(gene.bed)),with=F]
  
  
  write.table(gene.bed,
              file = paste(srcPath,"/annotations/bedtools/gene.bed",sep=""),
              quote = F, row.names = F, col.names = F, sep = "\t")
  
  # === load annotations, run bedtools: 
  rdataAnnots = list.files(path = paste(srcPath,"/annotations/", sep=  ""), full.names = T, pattern = '.Rdata')
  rdataAnnots

  for( a in 1:length(rdataAnnots)){
    
    cat("\nLoading annotation coordinates from file:",rdataAnnots[[a]])
    
    # Convert annotation file to bedtools format:
    annot.gr = readRDS(rdataAnnots[[a]])
    # mhcCoords = GRanges(seqnames="chr6", IRanges(start = 24e6, end = 35e6))
    # MHCannotations = findOverlaps(annot.gr, mhcCoords)
    # if(length(queryHits(MHCannotations)) > 0){ annot.gr = annot.gr[-queryHits(MHCannotations)]}
    
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
        
        
        cmd = paste(PathToBedtools,"bedtools intersect -a ",srcPath,"/annotations/bedtools/gene.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed -wo > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,"_",opt$out, ".txt" , sep = "")
        system(cmd, ignore.stdout = FALSE)
        
        
        suppressMessages(file.remove(paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed",sep="")))
        
      }
    } else {
      
      fwrite(annot.bed, file = paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed",sep=""), quote= F, row.names = F, col.names = F, sep = "\t")
      # Run bedtools variant annotation:
      cmd = paste(PathToBedtools,"bedtools intersect -a ",srcPath,"/annotations/bedtools/snp.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed -wo > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])), "_",opt$out, ".txt" , sep = "")
      system(cmd, ignore.stdout = FALSE)
      cmd = paste(PathToBedtools,"bedtools intersect -a ",srcPath,"/annotations/bedtools/gene.bed -b ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,".bed -wo > ",srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),"_CHUNK_",tc,"_",opt$out, ".txt" , sep = "")
      system(cmd, ignore.stdout = FALSE)
      suppressMessages(file.remove(paste(srcPath,"/annotations/bedtools/",gsub(".Rdata", "", basename(rdataAnnots[[a]])),".bed",sep="")))
      
    }
    
}

}




##### Statistical approaches (Genome-wide linear regression of intervals stratified by annotation, permutation analysis across p-value bins, binomial enrichment tests across p-value bins, SumVg)
# === SumVg citation: Hon-Cheong SO and  Pak C. SHAM (2015) "SumVg: Total heritability explained by all variants in genome-wide association studies based on summary statistics with standard error estimates"
fleetTest = function(){

  search = ls()[grepl("gwasMunged", ls())]

  if(length(search) > 0) {ld.df <- gwasMunged} else{cat("\nReading summary statistics file...\n"); ld.df = as.data.frame(fread(paste(srcPath,"/out/",opt$out,".sumstats",sep=""), header = TRUE, showProgress=FALSE)) }

  
  
  coef_trait = list()
  glm_trait  = list()

  rdataAnnots = list.files(path = paste(srcPath,'/annotations/bedtools/',sep=""), full.names = TRUE, pattern = '.txt')
  rdataAnnots = rdataAnnots[grepl(paste(opt$out,".txt",sep=""), rdataAnnots)]
  rdataAnnots
  
  # == Calculate median Chi-sq
  zall = log(ld.df$OR)/ld.df$se
  chisq = zall**2
  meanchisq = mean(chisq)
  lambdagc = median(chisq)/0.454
  cat("\nMean chi-square:",round(meanchisq,3))
  cat("\nLambda GC:",round(lambdagc,3))
  
  # == estimate heritability explained by intervals (SumVg method)
  
  EstVgAll = SumVg.binary(zall=zall, 
               method="jack",
               d=length(zall)/5,
               repl=5, 
               out="unconditional",
               caseNo=opt$Nca, 
               ctrlNo=opt$Nco, 
               K=opt$popprev)
  # == Total heritability:
  EstVgAll = data.frame(t(unlist(EstVgAll)))
  EstVgAll$VgZ = EstVgAll[,1]/EstVgAll[,2]
  EstVgAll$VgP = pnorm(-abs(EstVgAll$VgZ))
  EstVgAll$VgP = ifelse(EstVgAll$VgP == 0, 1e-300, EstVgAll$VgP)
  if(EstVgAll[,1] > 1.0){cat("\rWarning! Heritability estimate (liability scale) from intervals was greater than 1, possibly due to insufficient LD pruning")}
  cat("\nTotal heritability (liability scale):",round(EstVgAll[,1],4),"(SE:",round(EstVgAll$SE.SumVg,4),")")
  cat("\nHeritabiliy z-score: ",EstVgAll$VgZ)
  cat("\nHeritabiliy p-value: ",EstVgAll$VgP)
  
  
  for( a in 1:length(rdataAnnots)){

    cat("\nLoading variant annotations:",rdataAnnots[[a]],"...")
    
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
    } #!/close annotation by chromosome
    
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
    
    
    cat("\nTotal number of LD-intervals available for enrichment tests:", nrow(all_df),"\n")
    

    ## === linear model in fast mode
    if(opt$speed == 'fast'){
      
      dt_all_df = data.table(all_df)
      
      misCol = colSums(dt_all_df[,colnames(dt_all_df) %in% tracks,with=F] == 1)

      misCol = misCol[misCol < opt$`annot-cnt`]

      if(length(misCol) > 0){
        dt_all_df = dt_all_df[,!colnames(dt_all_df) %in% names(misCol),with=F]
        }
      
      keep_tracks = tracks[tracks %in% colnames(dt_all_df)]
      
      if(length(keep_tracks) < 1) next
      
      response_matrix = as.matrix(dt_all_df[,colnames(dt_all_df) %in% keep_tracks,with=F])
      response = dt_all_df$MeanZSCORE
      
      
      cat("\rFitting linear regressions in fast mode (Zscore ~ Set [multiple] + covariates)!\n")
      start_reg = proc.time()
      mod0 = model.matrix( ~. , data = dt_all_df[,colnames(dt_all_df) %in% c(keep_tracks, "MAFc", "LDSCORE", "LOCSIZE", "SNPDENSITY"),with=F])
      lmFast = lm(response ~ -1 + mod0)
      cat("\r...Time elapsed for regression:", (proc.time() - start_reg)[3], 'seconds\n')
      lmFastCoef = summary(lmFast)
      
          lmCoef = broom::tidy(summary(lmFast)$coefficients)
          lmCoef = lmCoef[!grepl("Intercept|MAFc|LDSCORE|SNPDENSITY|LOCSIZE", lmCoef$.rownames), ]
          lmCoef$.rownames = gsub("mod0", "", lmCoef$.rownames)
         
      lmDf = lmCoef
      
      lmDf = lmDf[,!colnames(lmDf) %in% ".id"]
      lmDf = lmDf[,!colnames(lmDf) %in% "Term"]
      colnames(lmDf) = c("SET_NAME",  "Beta", "SE", "T","P")
      
      lmDf$SET_SOURCE = SETNAME
      
      # == Jackknife delete-d 
      dt_all_df = dt_all_df[order(dt_all_df$LDSCORE), ]
      groups = dplyr::ntile(1:nrow(dt_all_df), 5)
      Ns = unique(groups)
      jack_list = list()
      for(j in Ns){
        cat("\rJackknife regression:",j)
        jack = dt_all_df[-which(groups %in% Ns[[j]])]

        response = jack$MeanZSCORE
        mod0 = model.matrix( ~. , data = jack[,colnames(jack) %in% c(keep_tracks, "MAFc", "LDSCORE", "SNPDENSITY"),with=F])
        lmFast = lm(response ~ -1 + mod0)

        lmFastCoef = summary(lmFast)

        lmCoef = broom::tidy(summary(lmFast)$coefficients)
        lmCoef = lmCoef[!grepl("Intercept|MAFc|LDSCORE|SNPDENSITY", lmCoef$.rownames), ]
        lmCoef$.rownames = gsub("mod0", "", lmCoef$.rownames)

        jackDf = lmCoef

        jackDf = jackDf[,!colnames(jackDf) %in% ".id"]
        jackDf = jackDf[,!colnames(jackDf) %in% "Term"]
        colnames(jackDf) = c("SET_NAME",  "Beta", "SE", "T","P")
        jack_list[[j]] = jackDf
      }
      jack_df = ldply(jack_list)
      jack_df = data.table(jack_df)

      meanVal = jack_df[,list(Mest = mean(Beta)),by=(SET_NAME)]
      n = length(Ns)
      for(y in 1:nrow(meanVal)){
        jackVals = jack_df$Beta[jack_df$SET_NAME %in% meanVal$SET_NAME[[y]]]
        # jackSE = (((n-1)/n) * sum((jackVals - meanVal$M[[y]])^2))^0.5
        u = jackVals
        n = length(jackVals)
        jackSE = sqrt(((n - 1)/n) * sum((u - mean(u))^2))
        meanVal$JACKSE[[y]] = round(jackSE,4)
        meanVal$JACKTval[[y]] = meanVal$M[[y]]/jackSE
        meanVal$JACKPval[[y]] = pnorm(-abs(meanVal$JACKTval[[y]]))
      }

     lmDf = merge(lmDf, meanVal,by="SET_NAME")
     #
      
      
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
      
      # lmDf = data.frame(lmDf, VgCalc_df[,!colnames(VgCalc_df) %in% "SET_NAME"])
      
      
      coef_df = lmDf
    } #!/close fast enrichment analysis
    
    
    # === linear model in iterative ("slow") mode
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
        
        lmFit[[x]] = fit
      }
      
      tidy = ldply(lmFit)
      
    
        SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
        SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
        
        tidy$GWSVARNAMES[[x]] = ifelse( gwas_vars == "", "*", gwas_vars)
        
        # Combine model results
        coef_df = tidy
        
      } #!/close slow enrichment analysis
    
      # == SumVg process
    if(opt$SumVg == T){
  
      # == Calculate Vg for SNPs in annotation category
      
      VgCalc_list = list()
      for(x in 1:length(keep_tracks)){
        
        sumVgInAnnot = which(dt_all_df[,colnames(dt_all_df) %in% keep_tracks[[x]],with=F] == 1)
        
        dt_all_df = dt_all_df[order(dt_all_df$LDSCORE,decreasing = TRUE)]
        zall = log(dt_all_df$OR)/dt_all_df$se
        zall = zall[sumVgInAnnot]
        
        uniqueZval = length(unique(zall))
        
        if(uniqueZval > 10){
          VgInAnnot = try(SumVg.binary(zall=zall, 
                                       method="jack",
                                       d=length(zall)/5,
                                       repl=5, 
                                       out="unconditional",
                                       caseNo=opt$Nca, 
                                       ctrlNo=opt$Nco, 
                                       K=0.01))
          if(class(VgInAnnot) != "try-error"){
            VgCalc = data.frame(t(unlist(VgInAnnot)))
            VgCalc$VgZ = VgCalc[,1]/VgCalc[,2]
            VgCalc$VgP = pnorm(-abs(VgCalc$VgZ))
            VgCalc$VgEnrichment = (VgCalc$Est.SumVg/EstVgAll$Est.SumVg)/(length(sumVgInAnnot)/nrow(dt_all_df))
            VgCalc$SEEnrichment = (VgCalc$SE.SumVg)/(length(sumVgInAnnot)/nrow(dt_all_df))
            VgCalc$N = length(sumVgInAnnot)
            VgCalc$SET_NAME = keep_tracks[[x]]
            VgCalc_list[[x]] = VgCalc
          } 
          if(class(VgInAnnot) == "try-error"){
            VgCalc_list[[x]] = data.frame(VgZ = NA, VgP = NA, VgEnrichment = NA, SEEnrichment = NA, N = NA, SET_NAME = keep_tracks[[x]])
          }
          cat("\rCalculating Vg for SNPs in annotation category:",round(x/length(keep_tracks)*100,2),"%")
          
        }
        
      }
      VgCalc_df = ldply(VgCalc_list)
      VgCalc_df = VgCalc_df[match(lmDf$SET_NAME, VgCalc_df$SET_NAME), ]
      
      coef_df = merge(coef_df, VgCalc_df, by="SET_NAME")
      
    } #!/close sumvg method
    
    
    
      
      # == binomial test across p-value bins 
    if(opt$binomial == TRUE){
      
      
      # === binomial enrichment analysis
      
      set_names = keep_tracks

      bsStats = list()
      for(x in 1:length(set_names)){
            
      phenotypes = dt_all_df[,colnames(dt_all_df) %in% set_names[[x]],with=F]
      expectedProp = table(phenotypes)[[2]]
      nSize = nrow(phenotypes)
      baseline = expectedProp/nSize
      
      if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(5e-08, 1e-06, 1e-05, 1e-04, 5e-02, 1e-1, 2e-1, 5e-1)}
      
      cat("\rBinomial tests:",round(x/length(set_names),3)*100,"%");flush.console()
      binomstats = list()
      for( z in 1:length(pthreshold)){
        
        gws.hits.1 = dt_all_df[dt_all_df$P < pthreshold[[z]]]
        hits1 = colSums(gws.hits.1[,set_names[[x]],with=FALSE])[[1]]
        
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
      
      coef_df = merge(coef_df, bsStats, by= "SET_NAME")
      
    } #!/close binomial method
    
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
        if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(1e-07, 1e-06, 1e-05)}
        
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



# == Fleet Test on gene scores

fleetGS = function(){
  
  search = ls()[grepl("geneScores", ls())]
  
  if(length(search) > 0) {ld.df <- geneScores} else{cat("\nReading summary statistics file...\n"); ld.df = as.data.frame(fread(paste(srcPath,"/out/",opt$out,".fleet.geneScores.txt",sep=""), header = TRUE, showProgress=FALSE)) }
  
  
  
  coef_trait = list()
  glm_trait  = list()
  
  rdataAnnots = list.files(path = paste(srcPath,'/annotations/bedtools/',sep=""), full.names = TRUE, pattern = '.txt')
  rdataAnnots = rdataAnnots[grepl(paste(opt$out,".txt",sep=""), rdataAnnots)]
  rdataAnnots
  
  for( a in 1:length(rdataAnnots)){
    
    cat("\nLoading variant annotations:",rdataAnnots[[a]],"...")
    
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
    } #!/close annotation by chromosome
    
    cast = rbindlist(cast_list, fill = TRUE)
    
    names(ld.df)[names(ld.df) %in% "symbol"] = "INDEX"
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
    all_df$LOGP = -log10(all_df$Pval)
    
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
    
    
    SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
    SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
    
    all_df = all_df[!is.infinite(all_df$ZSCORE), ]
    
    
    cat("\nTotal number of LD-intervals available for enrichment tests:", nrow(all_df),"\n")
    
    
    ## === linear model in fast mode
    if(opt$speed == 'fast'){
      
      dt_all_df = data.table(all_df)
      
      misCol = colSums(dt_all_df[,colnames(dt_all_df) %in% tracks,with=F] == 1)
      
      misCol = misCol[misCol < opt$`annot-cnt`]
      
      if(length(misCol) > 0){
        dt_all_df = dt_all_df[,!colnames(dt_all_df) %in% names(misCol),with=F]
      }
      
      keep_tracks = tracks[tracks %in% colnames(dt_all_df)]
      
      if(length(keep_tracks) < 1) next
      
      response_matrix = as.matrix(dt_all_df[,colnames(dt_all_df) %in% keep_tracks,with=F])
      response = dt_all_df$MeanZSCORE
      
      
      cat("\rFitting linear regressions in fast mode (Zscore ~ Set [multiple] + covariates)!\n")
      start_reg = proc.time()
      mod0 = model.matrix( ~. , data = dt_all_df[,colnames(dt_all_df) %in% c(keep_tracks, "MAFc", "LDSCORE", "LOCSIZE", "SNPDENSITY"),with=F])
      lmFast = lm(response ~ -1 + mod0)
      cat("\r...Time elapsed for regression:", (proc.time() - start_reg)[3], 'seconds\n')
      lmFastCoef = summary(lmFast)
      
      lmCoef = broom::tidy(summary(lmFast)$coefficients)
      lmCoef = lmCoef[!grepl("Intercept|MAFc|LDSCORE|SNPDENSITY|LOCSIZE", lmCoef$.rownames), ]
      lmCoef$.rownames = gsub("mod0", "", lmCoef$.rownames)
      
      lmDf = lmCoef
      
      lmDf = lmDf[,!colnames(lmDf) %in% ".id"]
      lmDf = lmDf[,!colnames(lmDf) %in% "Term"]
      colnames(lmDf) = c("SET_NAME",  "Beta", "SE", "T","P")
      
      lmDf$SET_SOURCE = SETNAME
      
      # == Jackknife delete-d 
      dt_all_df = dt_all_df[order(dt_all_df$LDSCORE), ]
      groups = dplyr::ntile(1:nrow(dt_all_df), 5)
      Ns = unique(groups)
      jack_list = list()
      for(j in Ns){
        cat("\rJackknife regression:",j)
        jack = dt_all_df[-which(groups %in% Ns[[j]])]
        
        response = jack$MeanZSCORE
        mod0 = model.matrix( ~. , data = jack[,colnames(jack) %in% c(keep_tracks, "MAFc", "LDSCORE", "SNPDENSITY"),with=F])
        lmFast = lm(response ~ -1 + mod0)
        
        lmFastCoef = summary(lmFast)
        
        lmCoef = broom::tidy(summary(lmFast)$coefficients)
        lmCoef = lmCoef[!grepl("Intercept|MAFc|LDSCORE|SNPDENSITY", lmCoef$.rownames), ]
        lmCoef$.rownames = gsub("mod0", "", lmCoef$.rownames)
        
        jackDf = lmCoef
        
        jackDf = jackDf[,!colnames(jackDf) %in% ".id"]
        jackDf = jackDf[,!colnames(jackDf) %in% "Term"]
        colnames(jackDf) = c("SET_NAME",  "Beta", "SE", "T","P")
        jack_list[[j]] = jackDf
      }
      jack_df = ldply(jack_list)
      jack_df = data.table(jack_df)
      
      meanVal = jack_df[,list(Mest = mean(Beta)),by=(SET_NAME)]
      n = length(Ns)
      for(y in 1:nrow(meanVal)){
        jackVals = jack_df$Beta[jack_df$SET_NAME %in% meanVal$SET_NAME[[y]]]
        # jackSE = (((n-1)/n) * sum((jackVals - meanVal$M[[y]])^2))^0.5
        u = jackVals
        n = length(jackVals)
        jackSE = sqrt(((n - 1)/n) * sum((u - mean(u))^2))
        meanVal$JACKSE[[y]] = round(jackSE,4)
        meanVal$JACKTval[[y]] = meanVal$M[[y]]/jackSE
        meanVal$JACKPval[[y]] = pnorm(-abs(meanVal$JACKTval[[y]]))
      }
      
      lmDf = merge(lmDf, meanVal,by="SET_NAME")
      #
      
      
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
      
      # lmDf = data.frame(lmDf, VgCalc_df[,!colnames(VgCalc_df) %in% "SET_NAME"])
      
      
      coef_df = lmDf
    } #!/close fast enrichment analysis
    
    
    # === linear model in iterative ("slow") mode
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
        
        lmFit[[x]] = fit
      }
      
      tidy = ldply(lmFit)
      
      
      SETNAME = strsplit(gsub(opt$out, "", basename(rdataAnnots[[a]])), "_CHUNK")[[1]][[1]]
      SETNAME = gsub("_|_.txt|.txt", "", SETNAME)
      
      tidy$GWSVARNAMES[[x]] = ifelse( gwas_vars == "", "*", gwas_vars)
      
      # Combine model results
      coef_df = tidy
      
    } #!/close slow enrichment analysis
    
    # == SumVg process
    if(opt$SumVg == T){
      
      # == Calculate Vg for SNPs in annotation category
      
      VgCalc_list = list()
      for(x in 1:length(keep_tracks)){
        
        sumVgInAnnot = which(dt_all_df[,colnames(dt_all_df) %in% keep_tracks[[x]],with=F] == 1)
        
        dt_all_df = dt_all_df[order(dt_all_df$LDSCORE,decreasing = TRUE)]
        zall = log(dt_all_df$OR)/dt_all_df$se
        zall = zall[sumVgInAnnot]
        
        uniqueZval = length(unique(zall))
        
        if(uniqueZval > 10){
          VgInAnnot = try(SumVg.binary(zall=zall, 
                                       method="jack",
                                       d=length(zall)/5,
                                       repl=5, 
                                       out="unconditional",
                                       caseNo=opt$Nca, 
                                       ctrlNo=opt$Nco, 
                                       K=0.01))
          if(class(VgInAnnot) != "try-error"){
            VgCalc = data.frame(t(unlist(VgInAnnot)))
            VgCalc$VgZ = VgCalc[,1]/VgCalc[,2]
            VgCalc$VgP = pnorm(-abs(VgCalc$VgZ))
            VgCalc$VgEnrichment = (VgCalc$Est.SumVg/EstVgAll$Est.SumVg)/(length(sumVgInAnnot)/nrow(dt_all_df))
            VgCalc$SEEnrichment = (VgCalc$SE.SumVg)/(length(sumVgInAnnot)/nrow(dt_all_df))
            VgCalc$N = length(sumVgInAnnot)
            VgCalc$SET_NAME = keep_tracks[[x]]
            VgCalc_list[[x]] = VgCalc
          } 
          if(class(VgInAnnot) == "try-error"){
            VgCalc_list[[x]] = data.frame(VgZ = NA, VgP = NA, VgEnrichment = NA, SEEnrichment = NA, N = NA, SET_NAME = keep_tracks[[x]])
          }
          cat("\rCalculating Vg for SNPs in annotation category:",round(x/length(keep_tracks)*100,2),"%")
          
        }
        
      }
      VgCalc_df = ldply(VgCalc_list)
      VgCalc_df = VgCalc_df[match(lmDf$SET_NAME, VgCalc_df$SET_NAME), ]
      
      coef_df = merge(coef_df, VgCalc_df, by="SET_NAME")
      
    } #!/close sumvg method
    
    
    
    
    # == binomial test across p-value bins 
    if(opt$binomial == TRUE){
      
      
      # === binomial enrichment analysis
      
      set_names = keep_tracks
      
      bsStats = list()
      for(x in 1:length(set_names)){
        
        phenotypes = dt_all_df[,colnames(dt_all_df) %in% set_names[[x]],with=F]
        expectedProp = table(phenotypes)[[2]]
        nSize = nrow(phenotypes)
        baseline = expectedProp/nSize
        
        if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(5e-08, 1e-06, 1e-05, 1e-04, 5e-02, 1e-1, 2e-1, 5e-1)}
        
        cat("\rBinomial tests:",round(x/length(set_names),3)*100,"%");flush.console()
        binomstats = list()
        for( z in 1:length(pthreshold)){
          
          gws.hits.1 = dt_all_df[dt_all_df$P < pthreshold[[z]]]
          hits1 = colSums(gws.hits.1[,set_names[[x]],with=FALSE])[[1]]
          
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
      
      coef_df = merge(coef_df, bsStats, by= "SET_NAME")
      
    } #!/close binomial method
    
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
        if(opt$pthres != ""){pthreshold = read.table(opt$pthres, header = F)[,2]} else{pthreshold=c(1e-07, 1e-06, 1e-05)}
        
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
  
  
} # close 


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

# xMHCregion() # always run

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
