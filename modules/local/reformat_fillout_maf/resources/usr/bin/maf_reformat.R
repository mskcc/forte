#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
})

usage <- function() {
    message("Usage:")
    message("maf_reformat.R --maf <input.maf> --output_maf <output.maf>")
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

required_args <- c("maf","output_maf")
if (length(setdiff(required_args,names(opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

maf <- fread(opt$maf, data.table = F)
maf <- maf %>%
    mutate(
        # Using the Hugo Symbol column as an ID column because GBCMS does not allow users to control which columns should be kept.
        Hugo_Symbol = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2)
    )

write.table(maf,file = opt$output_maf,quote = F,sep = '\t', row.names = F)
