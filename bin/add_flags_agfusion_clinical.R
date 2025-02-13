#!/usr/local/bin/Rscript
# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __version__     = "0.0.1"


suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
})

usage <- function() {
    message("Usage:")
    message("add_flags_agfusion_clinical.R --cff-file <file.cff> --agfusion-file <agfusion.tsv> --transcript_allowlist <transcript_allowlist.txt> --out-prefix <prefix>")
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

possible_args = c("cff", "agfusion","transcript_allowlist", "out_prefix")
if (length(setdiff(names(args_opt),possible_args)) > 0){
    message("Invalid options")
    usage()
    quit()
}

required_args <- c("cff","agfusion","transcript_allowlist", "out_prefix")
if (length(setdiff(required_args,names(args_opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

agfusion_file = args_opt$agfusion
cff_file = args_opt$cff
out_prefix = args_opt$out_prefix
transcript_allowlist_file <- args_opt$transcript_allowlist

cff = fread(cff_file) %>%
      select(gene5_transcript_id,gene3_transcript_id,gene5_breakpoint,gene3_breakpoint,gene5_renamed_symbol,gene3_renamed_symbol,cluster) %>%
      mutate(ID = paste(gene5_transcript_id,gene3_transcript_id,gene5_breakpoint,gene3_breakpoint))
agfusion_tab = fread(agfusion_file) %>%
               mutate(ID = paste(`5'_transcript`,`3'_transcript`,`5'_breakpoint`,`3'_breakpoint`))
transcript_allowlist <- fread(transcript_allowlist_file) %>% select(ensembl_transcript,refseq_transcript)


agfusion_tab <- agfusion_tab %>% mutate(Metafusion_Transcript_Pair = ID %in% cff$ID)

agfusion_tab[,ID:=NULL]


cff <- cff %>%
       select(gene5_breakpoint,gene3_breakpoint,gene5_renamed_symbol,gene3_renamed_symbol,cluster) %>%
       distinct()

agfusion_tab <- merge(agfusion_tab,
                      cff,
                      by.x = c("5'_gene","3'_gene","5'_breakpoint","3'_breakpoint"),
                      by.y = c("gene5_renamed_symbol","gene3_renamed_symbol","gene5_breakpoint","gene3_breakpoint"),
                      all.x = T,
                      all.y = F)

agfusion_tab <- merge(agfusion_tab,
                      transcript_allowlist,
                      by.x = "5'_transcript",
                      by.y = "ensembl_transcript",
                      all.x = T,
                      all.y = F)
setnames(agfusion_tab,"refseq_transcript","refseq_transcript_id_5")

agfusion_tab <- merge(agfusion_tab,
                      transcript_allowlist,
                      by.x = "3'_transcript",
                      by.y = "ensembl_transcript",
                      all.x = T,
                      all.y = F)
setnames(agfusion_tab,"refseq_transcript","refseq_transcript_id_3")

write.table(
    agfusion_tab,
    paste0(out_prefix, ".expanded_agfusion_transcripts.tsv"),
    row.names = F,
    append = F,
    quote = F,
    sep = "\t"
)

