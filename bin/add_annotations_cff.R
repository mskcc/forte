#!/usr/local/bin/Rscript
# __author__      = "Anne Marie Noronha"
# __email__       = "noronhaa@mskcc.org"
# __version__     = "0.0.1"


suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
})

usage <- function() {
    message("Usage:")
    message("add_annotations_cff.R --cff-file <file.cff> --agfusion-file <agfusion.tsv> --oncokb-file <oncokb.tsv> --out-prefix <prefix>")
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

possible_args = c("cff", "oncokb", "agfusion", "out_prefix")
if (length(setdiff(names(args_opt),possible_args)) > 0){
    message("Invalid options")
    usage()
    quit()
}

if (length(setdiff(possible_args,names(args_opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

oncokb_file = args_opt$oncokb
agfusion_file = args_opt$agfusion
cff_file = args_opt$cff
out_prefix = args_opt$out_prefix

cff = fread(cff_file)
final_cff_cols <- c(names(cff))
agfusion_tab = fread(agfusion_file) %>% select(c(`5'_transcript`,`3'_transcript`,`5'_breakpoint`,`3'_breakpoint`,Fusion_effect))
final_cff_cols <- c(final_cff_cols,"Fusion_effect")
if (!is.null(oncokb_file)){
    oncokb_tab = fread(oncokb_file) %>% select(-Fusion)
    final_cff_cols = c(final_cff_cols,names(oncokb_tab %>% select(-Tumor_Sample_Barcode)))
    cff <- merge(
        cff,
        oncokb_tab,
        by.x ="FID",
        by.y = "Tumor_Sample_Barcode",
        all.x = T,
        all.y=F
    )
}

cff <- merge(
    cff,
    agfusion_tab,
    by.x = c("gene5_transcript_id","gene3_transcript_id","gene5_breakpoint","gene3_breakpoint"),
    by.y = c("5'_transcript","3'_transcript","5'_breakpoint","3'_breakpoint"),
    all.x = T,
    all.y = T
)

cff <- as.data.frame(cff)[,c(final_cff_cols)]
#cff <- cff %>% mutate(!!final_cff_cols[34] := Fusion_effect) %>% select(-c(Fusion_effect))

write.table(
    cff,
    paste0(out_prefix, ".cff"),
    row.names = F,
    quote = F,
    sep = "\t",
    col.names = ! "V1" %in% final_cff_cols
)

