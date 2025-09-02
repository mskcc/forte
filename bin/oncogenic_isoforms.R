#!/usr/local/bin/Rscript
# __author__      = "Caryn Hale"
# __email__       = "halec@mskcc.org"
# __maintainer__ = "Alexandria Dymun (pintoa1@mskcc.org)"
# __version__     = "0.0.1"


suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(tidyr)
})

usage <- function() {
    message("Usage:")
    message(
        "oncogenic_isoforms.R --portcullis <*.pass.junctions.tab> --junctions <TARGET_reportable_junctions.tsv> --sample <sample_name>"
    )
}

args = commandArgs(TRUE)

if (is.null(args) | length(args) < 1) {
    usage()
    quit()
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x)
        scan(
            text = x,
            what = 'character',
            quiet = TRUE
        ))
    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z) {
        length(z) <- 2
        z
    })

    parsed_args <- structure(lapply(args_vals, function(x)
        x[2]),
        names = lapply(args_vals, function(x)
            gsub('-', '_', x[1])))
    parsed_args[!is.na(parsed_args)]
}

args_opt <- parse_args(paste(args, collapse = " "))

possible_args = c("portcullis", "junctions", "sample")
if (length(setdiff(names(args_opt), possible_args)) > 0) {
    message("Invalid options")
    usage()
    quit()
}

required_args <- c("portcullis", "junctions", "sample")
if (length(setdiff(required_args, names(args_opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

output_cols <- c(
    'TumorId',
    "event",
    "WT_junctions",
    "Onco_junctions",
    "percent_oncogenic",
    "score",
    "Gene1",
    "Transcript1",
    "Chr1",
    "Pos1",
    "Gene2",
    "Transcript2",
    "Chr2",
    "Pos2",
    "Annotation",
    "Note",
    "oncokb_sv_type",
    "Site1Description",
    "Site2Description",
    "Fusion",
    "TotalReadSupport",
    "CallMethod",
    "FrameCallMethod",
    "Position",
    "Significance",
    "action"
)


junctionsOI <- fread(args_opt$junctions, data.table = F)
tab <- fread(args_opt$portcullis, data.table = F)

if (nrow(tab) == 0) {
    report <- data.frame(matrix(ncol = length(output_cols), nrow = 0))
    colnames(report) <- output_cols
    drop <- data.frame(matrix(ncol = length(output_cols), nrow = 0))
    colnames(drop) <- output_cols
    write.table(
        report,
        file = paste0(args_opt$sample, "_oncogenic_isoforms.txt"),
        quote = F,
        row.names = F,
        sep = "\t"
    )
    write.table(
        drop,
        file = paste0(args_opt$sample, "_oncogenic_isoforms_dropped.txt"),
        quote = F,
        row.names = F,
        sep = "\t"
    )
}

pre_Report <- left_join(junctionsOI,
                        tab,
                        by = c(
                            "chr" = "refname",
                            "Pos1" = "start",
                            "Pos2" = "end"
                        )) %>% mutate(score = rel2raw)

toReport <- pre_Report %>%
    select(event, junction_type, nb_raw_aln)  %>%
    mutate(nb_raw_aln = if_else(is.na(nb_raw_aln), 0, nb_raw_aln))

toReport <- toReport %>% pivot_wider(names_from = junction_type,
                                     values_from = nb_raw_aln,
                                     values_fill = 0)
toReport <- toReport %>% mutate(
    percent_oncogenic = OI / (WT + OI),
    action = ifelse(percent_oncogenic > 0 &
                        WT != 0 , "REPORT", "drop")
) %>%
    rename(WT_junctions = WT, Onco_junctions = OI)

toReport <- left_join(toReport, junctionsOI[junctionsOI$junction_type == "OI", ], by = "event")

toReportF <- left_join(toReport, pre_Report[pre_Report$junction_type == "OI", c("event", "score")], by = c("event"))
toReportF$TumorId <- args_opt$sample
toReportF <- toReportF[, output_cols]

report <- toReportF %>% filter(action == "REPORT") %>% select(-action)
drop <- toReportF %>% filter(action != "REPORT") %>% select(-action)

write.table(
    report,
    file = paste0(args_opt$sample, "_oncogenic_isoforms.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
)
write.table(
    drop,
    file = paste0(args_opt$sample, "_oncogenic_isoforms_dropped.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
)