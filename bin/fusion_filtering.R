#!/usr/local/bin/Rscript
# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __contributor__ = "Caryn Hale (halec@mskcc.org)"
# __version__     = "0.0.2"


suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(purrr)
    library(tidyr)
})

usage <- function() {
    message("Usage:")
    message(
        "fusion_filtering.R --cff <*.final.cff> --starfusion <*.starfusion.abridged.coding_effect.tsv> --fusioncatcher <*.fusioncatcher.fusion-genes.txt> --arriba <*.fusions.tsv> --clinical_genes <clinical_genes.txt> --out_prefix <prefix>"
    )
}

args = commandArgs(TRUE)

fc_flags = c(
    "banned",
    "bodymap2",
    "hpa",
    "1000genomes",
    "cacg",
    "gtex",
    "conjoing",
    "paralogs",
    "distance1000bp",
    "ensembl_fully_overlapping",
    "ensembl_same_strand_overlapping",
    "mt",
    "pair_pseudo_genes",
    "refseq_fully_overlapping",
    "refseq_same_strand_overlapping",
    "rrna",
    "similar_symbols",
    "ucsc_fully_overlapping",
    "ucsc_same_strand_overlapping"
)
sf_flags = c(
    "GTEx_recurrent_StarF2019",
    "BodyMap",
    "DGD_PARALOGS",
    "Greger_Normal",
    "Babiceanu_Normal",
    "ConjoinG"
)
cis_sage_allow = c("FGFR3::TACC3",
                   "TACC3::FGFR3",
                   "FGFR2::TACC2",
                   "TACC2::FGFR2")

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

possible_args = c("cff",
                  "starfusion",
                  "fusioncatcher",
                  "arriba",
                  "clinical_genes",
                  "out_prefix")
if (length(setdiff(names(args_opt), possible_args)) > 0) {
    message("Invalid options")
    usage()
    quit()
}

required_args <- c("cff",
                   "starfusion",
                   "fusioncatcher",
                   "arriba",
                   "clinical_genes",
                   "out_prefix")
if (length(setdiff(required_args, names(args_opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

#format fusioncatcher false positive flags
check_fc_flags <- function(fusion_description) {
    out <- sapply(fc_flags, function(flag) {
        ifelse(grepl(flag, fusion_description), flag, NA)
    })
    if (any(!is.na(out))) {
        return(paste0("fc:", paste0(out[!is.na(out)], collapse = ",")))
    } else{
        return(NA)
    }
}

#format starfusion false positive flags
check_sf_flags <- function(annots) {
    out <- sapply(sf_flags, function(flag) {
        ifelse(grepl(flag, annots), flag, NA)
    })
    if (any(!is.na(out))) {
        return(paste0("sf:", paste0(out[!is.na(out)], collapse = ",")))
    } else{
        return(NA)
    }
}

# return tool codes as A S F and G
## return alphabetically
get_tool_code <- function(tools) {
    out <- case_when(
        tools == "arriba" ~ "A",
        tools == "starfusion" ~ "S",
        tools == "fusioncatcher" ~ "F",
        tools == "agfusion" ~ "G"
    )
    out <- paste0(sort(out), collapse = "")
    return(out)
}

# internal scoring system for agfusion frame statuses.
# up for debate
fusion_effect_values <- function(Fusion_effect) {
    case_when(
        Fusion_effect == "in-frame" ~ 10,
        Fusion_effect == "in-frame (with mutation)" ~ 9,
        Fusion_effect == "out-of-frame" ~ 8,
        Fusion_effect == "exon-exon" ~ 7,
        Fusion_effect == "exon-intron" ~ 6,
        Fusion_effect == "intron-exon" ~ 6,
        Fusion_effect == "intron-intron" ~ 5.5,
        Fusion_effect == "exon-3UTR" ~ 5,
        Fusion_effect == "exon-5UTR" ~ 5,
        Fusion_effect == "3UTR-exon" ~ 5,
        Fusion_effect == "5UTR-exon" ~ 5,
        Fusion_effect == "exon-5UTR (end)" ~ 5,
        Fusion_effect == "exon-3UTR (start)" ~ 5,
        Fusion_effect == "5UTR (end)-exon" ~ 5,
        Fusion_effect == "3UTR (start)-exon" ~ 5,
        Fusion_effect == "intron-3UTR" ~ 4,
        Fusion_effect == "intron-5UTR" ~ 4,
        Fusion_effect == "5UTR-intron" ~ 4,
        Fusion_effect == "3UTR-intron" ~ 4,
        Fusion_effect == "intron-3UTR (start)" ~ 4,
        Fusion_effect == "intron-5UTR (end)" ~ 4,
        Fusion_effect == "5UTR (end)-intron" ~ 4,
        Fusion_effect == "3UTR (start)-intron" ~ 4,
        Fusion_effect == "5UTR-3UTR" ~ 3,
        Fusion_effect == "5UTR-3UTR (start)" ~ 3,
        Fusion_effect == "5UTR (end)-3UTR" ~ 3,
        Fusion_effect == "5UTR (end)-3UTR (start)" ~ 3,
        Fusion_effect == "3UTR-5UTR" ~ 2,
        Fusion_effect == "3UTR (start)-5UTR" ~ 2,
        Fusion_effect == "3UTR-5UTR (end)" ~ 2,
        Fusion_effect == "3UTR (start)-5UTR (end)" ~ 2,
        Fusion_effect == "3UTR-3UTR" ~ 1,
        Fusion_effect == "3UTR (start)-3UTR" ~ 1,
        Fusion_effect == "3UTR-3UTR (start)" ~ 1,
        Fusion_effect == "3UTR (start)-3UTR (start)" ~ 1,
        Fusion_effect == "5UTR-5UTR" ~ 1,
        Fusion_effect == "5UTR (end)-5UTR" ~ 1,
        Fusion_effect == "5UTR-5UTR (end)" ~ 1,
        Fusion_effect == "5UTR (end)-5UTR (end)" ~ 1,
        is.na(Fusion_effect) ~ 0,
        .default = 0
    )
}

### select a single breakpoint per cluster
select_breakpoint <- function(cluster_df) {
    unique_brs <- length(unique(cluster_df$breakpoint_key))
    br <- cluster_df %>% group_by(breakpoint_key) %>% mutate(n_tools = n_distinct(tool)) %>% ungroup()
    br <- br %>% mutate(
        gene5_ensg = grepl("^ENSG", reann_gene5_symbol) |
            is.na(reann_gene5_symbol),
        gene3_ensg = grepl("^ENSG", reann_gene3_symbol) |
            is.na(reann_gene3_symbol),
        both_symbol = !gene5_ensg & !gene3_ensg,
        one_symbol = (!gene5_ensg &
                          gene3_ensg) | (gene5_ensg & !gene3_ensg)
    )

    ## if any of the annotations have at least one symbol, select that annotation
    ### select both symbol breakpoint first
    if (any(br$both_symbol) & unique_brs > 1) {
        br <- br %>% filter(both_symbol)
        unique_brs <-  length(unique(br$breakpoint_key))
    }
    ### select at least one symbol breakpoint
    if (any(br$one_symbol) & unique_brs > 1) {
        br <- br %>% filter(one_symbol)
        unique_brs <-  length(unique(br$breakpoint_key))
    }
    ### if any tool inframe metafusion use the inframe
    if (unique_brs > 1 &
        any(grepl("in-frame", br$Fusion_effect) |
            (br$Arriba_Inframe & !is.na(br$Arriba_Inframe)))) {
        br <- br %>% filter(grepl("in-frame", br$Fusion_effect) |
                                (Arriba_Inframe &
                                     !is.na(Arriba_Inframe)))
        unique_brs <- length(unique(br$breakpoint_key))
    }
    ### select breakpoint with most tools calling
    if (unique_brs > 1 & any(br$n_tools > 1)) {
        br <- br %>% filter(n_tools == max(n_tools))
        unique_brs <-  length(unique(br$breakpoint_key))
    }
    ### select breakpint with max read support
    if (unique_brs > 1 & length(unique(br$total_support)) > 1) {
        br <- br %>% filter(total_support == max(total_support))
        unique_brs <-  length(unique(br$breakpoint_key))
    }
    ### if there are still multiple breakpoints, pick the lowest genomic position by 5' then 3'
    if (unique_brs > 1) {
        br <- br %>% filter(gene5_breakpoint == min(gene5_breakpoint)) %>%
            filter(gene3_breakpoint == min(gene3_breakpoint))  %>%
            mutate(fusion_effect_val = fusion_effect_values(Fusion_effect)) %>%
            filter(fusion_effect_val == max(fusion_effect_val))
        unique_brs <-  length(unique(br$breakpoint_key))
    }


    ## if still multiple select annotation with lowest transcript ID, sometimes no transcript is selected
    br <- br %>% mutate(trans5_num = as.numeric(gsub("ENST", "", gene5_transcript_id)),
                        trans3_num = as.numeric(gsub("ENST", "", gene3_transcript_id))) %>%
        filter(case_when(is.na(trans5_num) ~ T, trans5_num == min(trans5_num) ~ T)) %>%
        filter(case_when(is.na(trans3_num) ~ T, trans3_num == min(trans3_num) ~ T))
    ### need to return all rows that match the breakpoint key for frame status evaluation
    return(unique(br$breakpoint_key))
}

# select a single annotation per cluster
get_fusion_info <- function(br) {
    ### situations where there is only one breakpoint, but multiple annotations
    ### Attempt to select annotation where one or both genes are clinical genes
    ### if no clinical genes found, follow logic for finding breakpoint
    ## num tools will always be the same as they are grouped by breakpoint
    ## check in frame status then check total_support
    if (any(c(br$reann_gene5_symbol, br$reann_gene3_symbol) %in% clinical_genes$V1)) {
        br <- br %>% filter(
            reann_gene5_symbol %in% clinical_genes$V1 |
                reann_gene3_symbol %in% clinical_genes$V1
        )
    }
    if (any(grepl("in-frame", br$Fusion_effect) |
            (br$Arriba_Inframe & !is.na(br$Arriba_Inframe)))) {
        br <- br %>% filter(grepl("in-frame", Fusion_effect) |
                                (Arriba_Inframe &
                                     !is.na(Arriba_Inframe)))
    }
    ### select annotation with max read support
    ### if still too many annotations, select via fusion_effect value
    br <- br %>% filter(total_support == max(total_support)) %>%
        mutate(fusion_effect_val = fusion_effect_values(Fusion_effect)) %>%
        filter(fusion_effect_val == max(fusion_effect_val))
    ## if still have multiple select annotation where at least one symbol is a real symbol
    br <- br %>% mutate(
        gene5_ensg = grepl("^ENSG", reann_gene5_symbol) |
            is.na(reann_gene5_symbol),
        gene3_ensg = grepl("^ENSG", reann_gene3_symbol) |
            is.na(reann_gene3_symbol),
        both_symbol = !gene5_ensg & !gene3_ensg,
        one_symbol = (!gene5_ensg &
                          gene3_ensg) | (gene5_ensg & !gene3_ensg)
    )
    ## if any of the remaining annotations have at least one symbol, select that annotation
    ## select both symbols first
    if (any(br$both_symbol)) {
        br <- br %>% filter(both_symbol)
    }
    ## select with at least one symbol
    if (any(br$one_symbol)) {
        br <- br %>% filter(one_symbol)
    }
    ## if still multiple select annotation with lowest transcript ID, sometimes no transcript is selected
    br <- br %>% mutate(trans5_num = as.numeric(gsub("ENST", "", gene5_transcript_id)),
                        trans3_num = as.numeric(gsub("ENST", "", gene3_transcript_id))) %>%
        filter(case_when(is.na(trans5_num) ~ T, trans5_num == min(trans5_num) ~ T)) %>%
        filter(case_when(is.na(trans3_num) ~ T, trans3_num == min(trans3_num) ~ T))

    five_prime_gene <- ifelse(is.na(unique(br$reann_gene5_symbol)), ".", unique(br$reann_gene5_symbol))
    three_prime_gene <- ifelse(is.na(unique(br$reann_gene3_symbol)), ".", unique(br$reann_gene3_symbol))
    Fusion_effect <- unique(br$Fusion_effect)
    tx5 <- unique(br$gene5_transcript_id)
    tx3 <- unique(br$gene3_transcript_id)
    reciprocal_cluster_id <- unique(br$reciprocal_cluster_id)
    reciprocal_cluster <- unique(br$reciprocal_cluster)

    out <- c(
        "Fusion" = paste0(five_prime_gene, "::", three_prime_gene),
        "tx5" = tx5,
        "tx3" = tx3,
        "Fusion_effect" = Fusion_effect,
        "reciprocal_cluster_id" = reciprocal_cluster_id,
        "reciprocal_cluster" = reciprocal_cluster
    )

    return(out)
}

# assess cluster, determine if dropping or reporting
get_cluster_action <- function(cluster_df) {
    ## All breakpoints in cluster have less than 5 reads supports is LOW_SUPPORT
    read_support <- ifelse(all(cluster_df$total_support < 5), "LOW_SUPPORT", NA)
    # Cluster only has a single caller is ONE_CALLER
    caller_count <- ifelse(length(unique(cluster_df$tool)) < 2, "ONE_CALLER", NA)
    # Neither gene in the fusion belongs to the clinical_gene list is NO_SIG_GENE
    clinical_gene <- ifelse(any(unique(
        c(
            cluster_df$reann_gene5_symbol,
            cluster_df$reann_gene3_symbol
        )
    ) %in% clinical_genes$V1), NA, "NO_SIG_GENE")
    # Any breakpoint in the cluster is identified as a false positive by a tool is FP
    fp <- ifelse(all(is.na(cluster_df$FP_flag)), NA, "FP")
    ## only cis_sage if not in the allow list of cis-sage
    cis_sage <- ifelse(any(grepl("cis", cluster_df$cluster)) &
                           !any(cluster_df$symbol_id %in% cis_sage_allow),
                       "CIS_SAGE",
                       NA)


    reason <- paste(na.omit(c(
        clinical_gene, read_support, caller_count, fp, cis_sage
    )), collapse = ",")
    action  <- case_when(
        reason == "" ~ "REPORT",
        reason == "NO_SIG_GENE" ~ "NOVEL",
        reason == "CIS_SAGE" ~ "READ_THROUGH",
        .default = "drop"
    )

    return(c('action' = action, 'reason' = reason))
}

format_final_out <- function(cluster_df) {
    cluster_action <- get_cluster_action(cluster_df)
    ### return all tools that belong to the cluster
    tool = get_tool_code(unique(cluster_df$tool))
    ### get the tools within the cluster that are in frame
    frame_status_cl_inframe <- unique(c(cluster_df$tool[cluster_df$Tool_Inframe], na.omit(ifelse(
        grepl("in-frame", cluster_df$Fusion_effect),
        "agfusion",
        NA
    ))))
    frame_status_cl <- get_tool_code(frame_status_cl_inframe)
    fp_flags <- ifelse(any(!is.na(cluster_df$FP_flag)), paste(na.omit(unique(
        cluster_df$FP_flag
    )), collapse = ";"), NA)
    fp_tools <- ifelse(any(!is.na(cluster_df$FP_flag)), paste(na.omit(unique(
        gsub("\\:.*", "", cluster_df$FP_flag)
    )), collapse = ","), NA)

    ### ensure getting high quality breakpoints for non-dropped fusions
    if (cluster_action['action'] != "drop") {
        cluster_df <- cluster_df[cluster_df$total_support >= 5, ]
    }
    br_key <- select_breakpoint(cluster_df)
    br <- cluster_df[cluster_df$breakpoint_key %in% br_key, ]
    ###if there are still multiple breakpoints after selection, change strand to be ?
    if (length(br_key) > 1) {
        br$Diff_strand5 <- any(br$gene5_strand != br$gene5_strand[1])
        br$Diff_strand3 <- any(br$gene3_strand != br$gene3_strand[1])
        if (any(br$Diff_strand3) & any(br$Diff_strand5)) {
            br_key <- unique(gsub("(\\-|\\+)", "?", br_key))
        }
        if (any(br$Diff_strand3)) {
            br_key <- unique(gsub("(\\-$|\\+$)", "?", br_key))
        }
        if (any(br$Diff_strand5)) {
            br_key <- unique(gsub("(\\-\\||\\+\\|)", "?|", br_key))
        }
    }

    ### get the fusion name related to the chosen breakpoint
    fusion <- get_fusion_info(br)

    ### get the tools with the particular breakpoint that are inframe
    frame_status_br_inframe <- unique(c(br$tool[br$Tool_Inframe], na.omit(ifelse(
        any(grepl("in-frame", br$Fusion_effect)), "agfusion", NA
    ))))
    frame_status_br <- get_tool_code(frame_status_br_inframe)
    ### need to select max total support if the same breakpoint exists in multiple tools
    total_support <- max(br$total_support)

    out <- data.frame(
        "sample" = unique(cluster_df$sample),
        "cluster" = unique(cluster_df$cluster),
        "tool" = tool,
        "fusion" = fusion[["Fusion"]],
        "total_support" = total_support,
        "breakpoint" = br_key,
        "action" = cluster_action[["action"]],
        "reason" = cluster_action[["reason"]],
        "fp_tools" = fp_tools,
        "fp_flags" = fp_flags,
        "frame_status_br" = frame_status_br,
        "frame_status_cl" = frame_status_cl,
        "tx5" = fusion[["tx5"]],
        "tx3" = fusion[["tx3"]],
        "Fusion_effect" = fusion[["Fusion_effect"]],
        "reciprocal_cluster" = fusion[["reciprocal_cluster"]],
        "reciprocal_cluster_id" = fusion[["reciprocal_cluster_id"]]
    )
    return(out)

}

cff = fread(args_opt$cff, data.table = F)
output_headers <- c(
    "sample",
    "cluster",
    "tool",
    "fusion",
    "total_support",
    "breakpoint",
    "action",
    "reason",
    "fp_tools",
    "fp_flags",
    "frame_status_br",
    "frame_status_cl",
    "tx5",
    "tx3",
    "Fusion_effect"
)
cvr_output_headers <- c(
    "TumorId",
    "Chr1",
    "Pos1",
    "Str1",
    "Chr2",
    "Pos2",
    "Str2",
    "Gene1",
    "Transcript1",
    "Site1Description",
    "Gene2",
    "Transcript2",
    "Site2Description",
    "Fusion",
    "TotalReadSupport",
    "CallMethod",
    "FrameCallMethod",
    "Note",
    "Annotation",
    'Position',
    'oncokb_sv_type',
    "Significance"
)

if (nrow(cff) == 0) {
    final_outputfile <- data.frame(matrix(nrow = 0, ncol = length(output_headers)))
    colnames(final_outputfile) <- output_headers
    final_outputfile_cvr <- data.frame(matrix(nrow = 0, ncol = length(cvr_output_headers)))
    colnames(final_outputfile_cvr) <- cvr_output_headers
    cis_sage_output <- data.frame(matrix(nrow = 0, ncol = length(output_headers)))
    colnames(cis_sage_output) <- output_headers

    write.table(
        final_outputfile_cvr,
        file = paste0(args_opt$out_prefix, "_filtered_fusions_cvr.tsv"),
        quote = F,
        row.names = F,
        sep = "\t"
    )

    write.table(
        final_outputfile,
        file = paste0(args_opt$out_prefix, "_filtered_fusions.tsv"),
        quote = F,
        row.names = F,
        sep = "\t"
    )
    write.table(
        cis_sage_output,
        file = paste0(args_opt$out_prefix, "_cis_sage_fusions.tsv"),
        quote = F,
        row.names = F,
        sep = "\t"
    )
    q()
}


starfusion = fread(args_opt$starfusion)
fusioncatcher = fread(args_opt$fusioncatcher)
arriba = fread(args_opt$arriba)
clinical_genes = fread(args_opt$clinical_genes, header = F)
### cff FIDs are matching order of arriba/starfusion/fusioncatcher row input
## this is important for when there are duplicated fusions within a tool (frequently happens with arriba)
### make a key to match up information between files
cff <- cff %>% mutate(FID_num = as.numeric(gsub("F", "", FID))) %>% arrange(FID_num) %>%
    mutate(
        key = paste0(gene5_chr, gene5_breakpoint, gene3_chr, gene3_breakpoint, tool),
        key = ifelse(duplicated(key), paste0(key, rowid(key)), key)
    )

if (nrow(arriba) != 0) {
    arriba <- arriba %>% mutate(
        key = paste0(
            "chr",
            gsub(":", "", breakpoint1),
            "chr",
            gsub(":", "", breakpoint2),
            "arriba"
        ),
        key = ifelse(duplicated(key), paste0(key, rowid(key)), key),
        total_support = split_reads1 + split_reads2 + discordant_mates,
        Tool_Inframe = grepl("in-frame", reading_frame),
        Arriba_Inframe = Tool_Inframe,
        FP_flag = NA
    )  %>% select(key, total_support, Tool_Inframe, Arriba_Inframe, FP_flag)
} else{
    arriba <- setNames(
        data.frame(matrix(ncol = 5, nrow = 0)),
        c(
            "key",
            "total_support",
            "Tool_Inframe",
            "Arriba_Inframe",
            "FP_flag"
        )
    )
}


if (nrow(fusioncatcher) != 0) {
    fusioncatcher <- fusioncatcher %>% mutate(
        key = paste0(
            "chr",
            gsub(
                "(:-|:\\+|:)",
                "",
                fusioncatcher$`Fusion_point_for_gene_1(5end_fusion_partner)`
            ),
            "chr",
            gsub(
                "(:-|:\\+|:)",
                "",
                fusioncatcher$`Fusion_point_for_gene_2(3end_fusion_partner)`
            ),
            "fusioncatcher"
        ),
        key = ifelse(duplicated(key), paste0(key, rowid(key)), key),
        total_support = Spanning_pairs,
        Arriba_Inframe = NA,
        Tool_Inframe = grepl("in-frame", Predicted_effect)
    ) %>%
        rowwise() %>% mutate(FP_flag = check_fc_flags(Fusion_description)) %>% ungroup()  %>% select(key, total_support, Tool_Inframe, Arriba_Inframe, FP_flag)
} else{
    fusioncatcher <- setNames(
        data.frame(matrix(ncol = 5, nrow = 0)),
        c(
            "key",
            "total_support",
            "Tool_Inframe",
            "Arriba_Inframe",
            "FP_flag"
        )
    )
}

if (nrow(starfusion) != 0) {
    starfusion <- starfusion %>% mutate(
        key = paste0(
            gsub("(:-|:\\+|:)", "", starfusion$LeftBreakpoint),
            gsub("(:-|:\\+|:)", "", starfusion$RightBreakpoint),
            "starfusion"
        ),
        key = ifelse(duplicated(key), paste0(key, rowid(key)), key),
        total_support = JunctionReadCount + SpanningFragCount,
        Arriba_Inframe = NA,
        Tool_Inframe = grepl("INFRAME", PROT_FUSION_TYPE)
    ) %>%
        rowwise() %>% mutate(FP_flag = check_sf_flags(annots)) %>% ungroup() %>% select(key, total_support, Tool_Inframe, Arriba_Inframe, FP_flag)
} else{
    starfusion <- setNames(
        data.frame(matrix(ncol = 5, nrow = 0)),
        c(
            "key",
            "total_support",
            "Tool_Inframe",
            "Arriba_Inframe",
            "FP_flag"
        )
    )
}


tool_info <- rbind(arriba, fusioncatcher, starfusion)
cff <- merge(cff, tool_info, by = "key")
cff <- cff  %>% mutate(
    breakpoint_key = paste0(
        gene5_chr,
        ":",
        gene5_breakpoint,
        ":",
        gene5_strand,
        "|",
        gene3_chr,
        ":",
        gene3_breakpoint,
        ":",
        gene3_strand
    ),
    symbol_id = paste(reann_gene5_symbol, reann_gene3_symbol, sep =
                          "::"),
    reciprocal_id = paste(reann_gene3_symbol, reann_gene5_symbol, sep = "::")
) %>% arrange(cluster) %>%
    group_by(cluster, symbol_id)  %>%  mutate(reciprocal_cluster_id = cur_group_id())


cff$reciprocal_cluster <- sapply(cff$key, function(key) {
    row <- cff[cff$key == key, ]
    cluster_df <- cff[cff$cluster == row$cluster, ]
    if (any(cluster_df$reciprocal_id == row$symbol_id)) {
        return(min(
            c(
                cluster_df$reciprocal_cluster_id[cluster_df$reciprocal_id == row$symbol_id],
                row$reciprocal_cluster_id
            )
        ))
    } else{
        return(NA)
    }
})
## after grouping turn from tibble back to data frame (mandatory for filtering)
cff <- as.data.frame(cff)

final_outputfile <- map_dfr(unique(cff$cluster), function(cluster_info) {
    ## reciprocal clusters need to be handled separately and carefully (do not wnat to return reciprocalcluster when a better breakpoint annotation exists)
    cluster_df <- cff[cff$cluster == cluster_info, ]
    out <- format_final_out(cluster_df)

    reciprocal_cluster_exists <- cluster_df[cluster_df$reciprocal_cluster == out[["reciprocal_cluster"]] &
                                                !is.na(cluster_df$reciprocal_cluster), ]
    if (length(unique(reciprocal_cluster_exists$reciprocal_cluster_id)) > 1) {
        cluster_df_2 <- reciprocal_cluster_exists[reciprocal_cluster_exists$reciprocal_cluster_id != out[["reciprocal_cluster_id"]], ]
        out2 <- format_final_out(cluster_df_2)

        out <- rbind(out, out2)
    }
    out <- out[, output_headers]
    return(out)
})

## separate cis sage clusters, if in cis sage allowlist, send to final output keep cis_sage cluster if clinical gen ****************
cis_sage_output <- final_outputfile %>% filter(
    grepl("cis_sage", cluster) &
        !fusion %in% cis_sage_allow & grepl("NO_SIG_GENE", reason)
)
final_outputfile <- final_outputfile %>% filter(!grepl("cis_sage", cluster) |
                                                    fusion %in% cis_sage_allow |
                                                    (grepl("cis_sage", cluster) &
                                                         !grepl("NO_SIG_GENE", reason)))

final_outputfile_cvr <- final_outputfile %>% filter(action == "REPORT") %>% select(sample,
                                                                                   tool,
                                                                                   fusion,
                                                                                   total_support,
                                                                                   breakpoint,
                                                                                   frame_status_cl,
                                                                                   tx5,
                                                                                   tx3)
final_outputfile_cvr <- final_outputfile_cvr %>% mutate(breakpoint = gsub("chr","",breakpoint)) %>% separate_wider_delim(fusion, "::", names = c("Gene1", "Gene2")) %>%
    separate_wider_delim(breakpoint, "|", names = c("bp1", "bp2")) %>%
    separate_wider_delim(bp1, ":", names = c("Chr1", "Pos1", "Str1")) %>%
    separate_wider_delim(bp2, ":", names = c("Chr2", "Pos2", "Str2"))
setnames(
    final_outputfile_cvr,
    c(
        'sample',
        "tx5",
        "tx3",
        "frame_status_cl",
        "total_support",
        "tool"
    ),
    c(
        "TumorId",
        "Transcript1",
        "Transcript2",
        "FrameCallMethod",
        "TotalReadSupport",
        "CallMethod"
    )
)
add_these <- setdiff(cvr_output_headers, colnames(final_outputfile_cvr))
final_outputfile_cvr[, add_these] <- NA
final_outputfile_cvr <- final_outputfile_cvr[, cvr_output_headers]

write.table(
    final_outputfile,
    file = paste0(args_opt$out_prefix, "_filtered_fusions.tsv"),
    quote = F,
    row.names = F,
    sep = "\t"
)

write.table(
    final_outputfile_cvr,
    file = paste0(args_opt$out_prefix, "_filtered_fusions_cvr.tsv"),
    quote = F,
    row.names = F,
    sep = "\t"
)

write.table(
    cis_sage_output,
    file = paste0(args_opt$out_prefix, "_cis_sage_fusions.tsv"),
    quote = F,
    row.names = F,
    sep = "\t"
)