#!/usr/local/bin/Rscript

# __author__      = "Anne Marie Noronha"
# __email__       = "noronhaa@mskcc.org"
# __contributor__ = ""
# __version__     = "0.1.0"
# __status__      = "Dev"


suppressPackageStartupMessages({
    library(data.table)
})

usage <- function() {
    message("Usage:")
    message("fingerprint_heatmap.R --input_table <inputtable> [--output <output>]")
}

args <- commandArgs(trailingOnly = TRUE)

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

opt <- parse_args(paste(args,collapse=" "))

required_args <- c("input_table")

if (length(setdiff(required_args,names(opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

if (! "output" %in% names(opt)){
    opt$output <- "heatmap.tsv"
}

combined_data <- list()
samplenames <- c()
file_paths <- fread(opt$input_table,header=T)

if (length(setdiff(c("sample","path"),names(file_paths))) > 0) {
    message("input to --input_table should have columns 'sample' and 'path'")
    quit()
}

for (i in 1:nrow(file_paths)) {
    sample_name <- file_paths$sample[i]
    file_path <- file_paths$path[i]

    # Check if the file exists before attempting to read it
    if (file.exists(file_path)) {
        data <- fread(file_path)
        data[, paste0(sample_name,"_Depth") := sum(as.integer(gsub(".*:","",strsplit(get(paste0(sample_name,"_Counts"))," ")[[1]]))) ]
        print(head(data))
        combined_data[[sample_name]] <- data
        samplenames <- c(samplenames,sample_name)
    } else {
        warning(paste("File not found:", file_path))
        combined_data[[sample_name]] <- NULL # or NA, or any other placeholder
    }
}


combined_data_tab <- do.call(cbind, combined_data)
#print(length(combined_data))

#print(dim(combined_data_tab))

combinations <- combn(samplenames, 2)

comparison.mtx <- matrix(1, nrow = length(samplenames), ncol = length(samplenames))
rownames(comparison.mtx) <- samplenames
colnames(comparison.mtx) <- samplenames

for (j in 1:dim(combinations)[[2]]){
    sampleA <- combinations[1,j]
    sampleB <- combinations[2,j]
    keepCols <- c(paste0(sampleA,".",sampleA,"_Genotypes"),paste0(sampleB,".",sampleB,"_Genotypes"), paste0(sampleA,".",sampleA,"_Depth"),paste0(sampleB,".",sampleB,"_Depth"))
    slice <- combined_data_tab[,..keepCols]
    slice <- data.frame(slice)
    slice <- slice[slice[,1] != "--" & slice[,2] != "--"  & slice[,3] > 20 & slice[,4] > 20,]

    total <- dim(slice)[1]
    match <- sum(slice[,1] == slice[,2])
    comparison.mtx[sampleA,sampleB] <- match/total
    comparison.mtx[sampleB,sampleA] <- match/total

}

write.table(comparison.mtx, opt$output , col.names=NA, quote = FALSE, sep="\t")
