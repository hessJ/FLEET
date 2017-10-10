#!/usr/bin/env Rscript


# Import packages
suppressMessages(library(optparse))
suppressMessages(require(argparse))
suppressMessages(require(GenomicRanges))
suppressMessages(require(data.table))

## Auto-detect the directory name for fleet.R
initial.options <- commandArgs(trailingOnly = FALSE)

# set masthead

masthead = as.character("
==========================================================
||
|| Functional LD-interval EnrichmEnt Test (FLEET)
|| *** Custom annotation extension script ***
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

# create parser object
parser <- ArgumentParser()

option_list = list(
  
  
  make_option(c( "--file"), type="character", default="", 
              help="Path to file (.txt, .csv, .tsv, etc.) containing annotations (with chromosome and range)", metavar="character"),
  
  make_option(c( "--chr"), type="character", default="", 
              help="Colname with chromosome", metavar="character"),
  
  make_option(c( "--start"), type="character", default="", 
              help="Colname for start of range (base pair)", metavar="character"),
  
  make_option(c( "--end"), type="character", default="", 
              help="Colname for end of range (base pair)", metavar="character"),
  
  make_option(c( "--id"), type="character", default="", 
              help="Colname with annotation labels", metavar="character"),
  
  make_option(c( "--out"), type="character", default="", 
              help="Output path for convert annotation table (GRanges object)", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


convertAnnot = function() {
  
  cat("\nImporting annotation table...")
  guessdelim = readLines(opt$file,1)
  tab = length(strsplit( guessdelim, "[\t]")[[1]])
  space = length(strsplit( guessdelim, " ")[[1]])
  comma = length(strsplit( guessdelim, "[,]")[[1]])
  
  set = list(tab,space,comma)
  names(set) = c("tab","space","comma")
  delim_type = c("\t"," ", ",")
  guesser = which(set %in% max(c(tab,space,comma)))
  
  cat("\nGuessing that file is:",names(set[guesser]),"delimited")
  
  read_annot = as.data.frame(fread(opt$file, h=T, sep = delim_type[guesser], fill = T))
  if(ncol(read_annot) < 2){stop("Unexpected number of fields in:",opt$file)}
  
  if(unique(grepl("chr", read_annot[,colnames(read_annot) %in% opt$chr])) == FALSE){read_annot[,colnames(read_annot) %in% opt$chr] = paste("chr",read_annot[,colnames(read_annot) %in% opt$chr],sep="")}
  
  cat("\nConverting to GRanges object...")
  
  annot_gr = GRanges(seqnames = read_annot[,colnames(read_annot) %in% opt$chr], 
                     IRanges(start = read_annot[,colnames(read_annot) %in% opt$start], 
                             end = read_annot[,colnames(read_annot) %in% opt$end]),
                     id = read_annot[,colnames(read_annot) %in% opt$id])
  

  
  cat("\nWriting output...")
  saveRDS(annot_gr, file = paste(opt$out,".Rdata",sep=""))
  
  
}

convertAnnot()
