#!/usr/local/bin/Rscript

# __author__      = "Alexandria Dymun"
# __maintainer__  = "Anne Marie Noronha"
# __email__       = "pintoa1@mskcc.org;noronhaa@mskcc.org"
# __version__     = "0.0.1"

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
    fillout_maf <- fread(opt$fillout_maf, data.table=FALSE) %>%
        rename(loc_spec = Hugo_Symbol)


    fillout_maf$rna_genotyped_variant_frequency <- fillout_maf$t_variant_frequency
    for(col in colnames(fillout_maf)[grepl('^t_',colnames(fillout_maf))]){
        new_col <- paste0('rna_',col)
        fillout_maf[new_col] <- fillout_maf[col]
    }


    fillout_maf <- fillout_maf[,c('loc_spec',colnames(fillout_maf)[grepl("^rna_",colnames(fillout_maf))])]

    maf <- fread(opt$maf, data.table = FALSE) %>%
        mutate(
            loc_spec = str_c(Chromosome,':',Start_Position,':',End_Position,':',Reference_Allele,':',Tumor_Seq_Allele2)
        )
    merged_maf <- merge(maf, fillout_maf, by = c("loc_spec"), all = T)
    merged_maf$RNA_ID <- opt$rna_sample_id
    merged_maf <- merged_maf %>% select(-c(loc_spec))
    return(merged_maf)

}

fillex <- rna_fillout_formatting(opt)

write.table(fillex,opt$output_maf,quote = FALSE, row.names = FALSE,sep = "\t")
