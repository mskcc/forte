#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

usage <- function() {
    message("Usage:")
    message("rna_fillout_reformat.R --fillout_maf <fillout.maf> --maf <original.maf> --output_maf <output.maf> --rna_sample_id <rna_sample_id>")
}

#negate
'%nin%' = Negate('%in%')
opt = commandArgs(TRUE)

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))
    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) gsub('-','_',x[1])))
    parsed_args[! is.na(parsed_args)]
}

opt <- parse_args(paste(opt,collapse=" "))

required_args <- c("fillout_maf","maf","output_maf","rna_sample_id")
if (length(setdiff(required_args,names(opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

rna_fillout_formatting <- function(opt){
  fillout_maf <- fread(opt$fillout_maf, data.table=FALSE)

  fillout_maf$rna_genotyped_variant_frequency <- fillout_maf$t_variant_frequency
  #fillout_maf$t_variant_frequency <- NULL
  for(col in colnames(fillout_maf)[grepl('^t_',colnames(fillout_maf))]){
    new_col <- paste0('rna_',col)
    fillout_maf[new_col] <- fillout_maf[col]
    #fillout_maf[col] <- NULL # why are we editing the old column? isn't it just dropped anyway?
  }

  #fillout_maf <- fillout_maf  %>% mutate(loc_spec = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele1))
  merge_cols <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1")
  fillout_maf <- fillout_maf[,c(merge_cols,colnames(fillout_maf)[grepl("^rna_",colnames(fillout_maf))])]

  maf <- fread(opt$maf, data.table = FALSE)
  merged_maf <- merge(maf,fillout_maf, by = merge_cols, all = T) # merged on multiple columns instead of loc_spec. wouldn't this work?
  merged_maf$RNA_ID <- opt$rna_sample_id
  return(merged_maf)

}

fillex<- rna_fillout_formatting(opt)

write.table(fillex,opt$output_maf,quote = FALSE, row.names = FALSE,sep = "\t")
