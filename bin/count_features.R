#!/usr/local/bin/Rscript

# __author__      = "Anne Marie Noronha"
# __email__       = "noronhaa@mskcc.org"
# __version__     = "0.0.1"

usage <- function() {
    message("Usage:")
    message("count_features.R --abundance <abundance.tsv> --gtf <annotation.gtf> --sample <sampleid>")
}

args = commandArgs(TRUE)

if (is.null(args) | length(args)<1) {
    usage()
    quit()
}

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

args_opt <- parse_args(paste(args,collapse=" "))

abundance <- read.table(args_opt$abundance, header=T)
gtf <- args_opt$gtf
Sample <- args_opt$sample

gtf_df <- as.data.frame(rtracklayer::import(gtf))

gtf_df <- unique(gtf_df[,c("transcript_id","gene_id")])

abundance <- merge(
    abundance,
    gtf_df,
    by.x="target_id",
    by.y="transcript_id",
    all.x = T,
    all.y =F
)
TranscriptCount <- dim(abundance[abundance$est_counts > 0,])[[1]]

abundance_gene <- aggregate(
    x = abundance$est_counts,
    by = list(abundance$gene_id),
    FUN = sum
)
GeneCount <- dim(abundance_gene[abundance_gene$x > 0,])[[1]]

out_df <- data.frame(Sample, TranscriptCount, GeneCount)
write.table(
    out_df,
    paste0(Sample,".kallisto.customsummary.txt"),
    sep="\t",
    row.names = F,
    append = F,
    quote = F
)
