

#!/usr/bin/env Rscript

# set masthead

masthead = as.character("
  ==========================================================
  *
  * Functional LD-clump EnrichmEnt Test (FLEET)
  *
  * Jonathan L. Hess, PhD and Stephen J. Glatt, PhD (c) 2017
  *
  * SUNY Upstate Medical University, PsychGENe Lab
  *
  * Contact: hessjo@upstate.edu
  *
  * GNU GENERAL PUBLIC LICENSE v3
  ===========================================================\n")

cat(masthead) # display masthead 


# Import packages

lib_list = list("outparse", "argparse", "plyr", "data.table", "GenomicRanges", "FDb.UCSC.snp137common.hg19",
  "TxDb.Hsapiens.UCSC.hg19.knownGene", "ggplot2", "gridExtra", "corrplot", "ggrepel", "reshape2",
  "foreach", "doParallel")


suppressMessages(library(optparse))
suppressMessages(require(argparse))

# ========================================================================
# create parser object
parser <- ArgumentParser()

option_list = list(
  make_option(c("-g", "--gwas"), type="character", default="", 
              help="Path to GWAS summary statistics", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default = %default]", metavar="character"),
  
  make_option(c("-r2", "--r2"), type="double", default=0.2, 
              help="R-squared threshold for linkage disequilibrium calculations [default = %default]", metavar="double"),
  
  make_option(c("-ldw", "--ld-window"), type="integer", default=1000, 
              help="Size of window (kilobases) for calculating linkage disequilibrium [default = %default]", metavar="integer"),
  
  make_option(c("-s", "--snp-field"), type="character", default="", 
              help="SNP column header in GWAS file", metavar="character"),
  
  make_option(c("-c", "--clump-field"), type="character", default="", 
              help="P-value column header in GWAS file", metavar="character"),
  
  
  make_option(c("-p", "--permStartThreshold"), type="double", default=5e-03, 
              help="Minimum P-value from enrichment test to initiate permutation analysis [default = %default]", metavar="double"),
  
  make_option(c("-n", "--nPerms"), type="integer", default=1000, 
              help="Number of permutations to perform [default = %default]", metavar="integer"),
  
  make_option(c("-t", "--threads"), type="integer", default=1, 
              help="Number of cores for parallelization [default = %default]", metavar="integer"),
  
  make_option(c("-l", "--label-annotations"), type="character", default="annotations/annotation.txt", 
              help="Path to annotation table [default = %default]", metavar="character"),
  
  make_option(c("-d", "--rd-annots"), type="character", default="annotations/", 
              help="Path to .Rdata annotations [default = %default]", metavar="character")
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# handling empty arguments
if (opt$gwas == ""){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}
# ========================================================================

if(opt$gwas != ""){
  cat("\nPackages loaded....")
  cat("\nPath to GWAS file:", opt$gwas, "\n")
}

# //// Dependencies

suppressMessages(library(plyr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(FDb.UCSC.snp137common.hg19))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(corrplot))
suppressMessages(library(ggrepel))
suppressMessages(library(reshape2))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))



# //// Functions 

# ========================================================================
# Import SNP coordinates from 1KG reference:


readSNPcoords  = function(){
  cat("Reading SNP coordinates from 1KG bim files")
  read_in = list.files(path = 'qc/ref/', full.names = T, pattern = 'bim')
  read_in = lapply(read_in, function(x) fread(x, header = FALSE))
  read_df = rbindlist(read_in)
  dups = read_df[duplicated(read_df$V2), ]
  write.table(dups$V2, file = 'qc/ref/duplicatevars.txt', quote = F, row.names = F, col.names = F)
  colnames(read_df) = c("seqnames", "SNP", "empty", "start", "A1", "A2")
  read_df = as.data.frame(read_df)
  read_df$seqnames = paste("chr", read_df$seqnames, sep = "")

  snp137_df <<- read_df
}


# pull out snps in the xMHC region
xMHCregion = function(){
  mhc =  snp137_df[snp137_df$seqnames %in% "chr6" & snp137_df$start >= 24e6 & snp137_df$start <= 35e6, ]
  snp137_df_nomhc <<- snp137_df[!snp137_df$SNP %in% mhc$SNP, ]
}

# Set annotations for FLEET
readAnnots = function(){
  read_annots = read.table(opt$annots, header = T, sep = "\t")
  annot.cat <<- read_annots
  
  # path to annotation .Rdata files
  rdata_files = list.files(path = opt$`rd-annots`, full.names = T, pattern = ".Rdata")
  
  if(nrow(read_annots) != length(rdata_files)){stop("Number of .Rdata files does not match table with annotation labels!")}
}



