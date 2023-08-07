#!/usr/local/bin/Rscript
# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __contributor__ = "Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "0.0.1"
# __status__      = "Dev"



library(dplyr)
library(data.table)
    args <- commandArgs(TRUE)
    if (length(args) != 5) {
        stop(
            "5 arguments are required as input in the following order: filtered_cff cluster_file cis_sage_file problematic_chromosomes_file sample_name"
        )
    }

    filtered_cff <- fread(args[1],data.table = F)
    header_filtered <-
        c(
            "gene5_chr",
            "gene5_breakpoint",
            "gene5_strand",
            "gene3_chr",
            "gene3_breakpoint",
            "gene3_strand",
            "library",
            "sample",
            "T_N",
            "disease",
            "tool",
            "max_split_cnt",
            "max_span_cnt",
            "gene5_renamed_symbol",
            "gene5_tool_annotation",
            "gene3_renamed_symbol",
            "gene3_tool_annotation",
            "FusionType",
            "reann_gene5_symbol",
            "reann_gene5_region",
            "reann_gene3_symbol",
            "reann_gene3_region",
            "reann_gene5_on_bndry",
            "reann_gene5_close_to_bndry",
            "reann_gene3_on_bndry",
            "reann_gene3_close_to_bndry",
            "score",
            "coding_id_distance",
            "gene_interval_distance",
            "dnasupp",
            "FID",
            "gene5_seq",
            "gene3_seq",
            "is_inframe",
            "closest_exon5",
            "closest_exon3",
            "captured_reads",
            "gene5_transcript_id",
            "gene3_transcript_id"
        )
    colnames(filtered_cff) <- header_filtered
    cluster <- fread(args[2],data.table = F)
    header_cluster <-
    c(
        "gene5_renamed_symbol",
        "gene3_renamed_symbol",
        "gene5_chr",
        "gene5_breakpoint",
        "gene3_chr",
        "gene3_breakpoint",
        "max_split_cnt",
        "max_span_cnt",
        "T_N",
        "disease",
        "tool",
        "FusionType",
        "sample",
        "cancer_db_hits",
        "FID"
    )
    colnames(cluster) <- header_cluster

    sample_name <- args[5]

    cluster_fids <- strsplit(cluster$FID, ",")
    df_cluster <- data.frame(FID = unlist(cluster_fids),
                            cluster = rep(seq_along(cluster_fids), lapply(cluster_fids, length)))

    cis_sage <- tryCatch({fread(args[3],data.table = F)},warning = function(cond){return( NULL)})
    if(!is.null(cis_sage)){
        header_cis <-
        c(
            "test",
            "gene5_renamed_symbol",
            "gene3_renamed_symbol",
            "max_split_cnt",
            "max_span_cnt",
            "T_N",
            "disease",
            "tool",
            "FusionType",
            "reann_gene5_on_bndry",
            "reann_gene5_close_to_bndry",
            "reann_gene3_on_bndry",
            "reann_gene3_close_to_bndry",
            "is_inframe",
            "sample",
            "gene5_chr",
            "gene5_breakpoint",
            "gene3_chr",
            "gene3_breakpoint",
            "closest_exons5",
            "closest_exons3",
            "FID"
        )

        colnames(cis_sage) <- header_cis
        cis_fids <- strsplit(cis_sage$FID, ",")
        df_cis <- data.frame(FID = unlist(cis_fids),
                        cluster = paste0("cis_sage_", rep(
                        seq_along(cis_fids), lapply(cis_fids, length)
                        )))
        df_cluster <- rbind(df_cluster, df_cis)
        filtered_cff <- merge(filtered_cff, df_cluster)
    }

    filtered_cff$Metafusion_flag <- apply(filtered_cff, 1, function(row) {
    if (is.na(row["reann_gene5_symbol"]) ||
        is.na(row["reann_gene3_symbol"])) {
        return("Gene_or_loc_not_in_bed")
    } else if (row["gene3_renamed_symbol"] != row["reann_gene3_symbol"] ||
                row["gene5_renamed_symbol"] != row["reann_gene5_symbol"]) {
        return("Annot_gene_different_from_renamed")
    } else{
        return(NA)
    }
    })

    weird_chromosomes <- tryCatch({fread(args[4],data.table = F)},warning = function(cond){return( NULL)})
    if(!is.null(weird_chromosomes)){
        colnames(weird_chromosomes) <-
        c(
            "gene5_chr",
            "gene5_breakpoint",
            "gene5_strand",
            "gene3_chr",
            "gene3_breakpoint",
            "gene3_strand",
            "library",
            "sample",
            "T_N",
            "disease",
            "tool",
            "max_split_cnt",
            "max_span_cnt",
            "gene5_renamed_symbol",
            "gene5_tool_annotation",
            "gene3_renamed_symbol",
            "gene3_tool_annotation"
        )
        weird_chromosomes[, colnames(filtered_cff)[!colnames(filtered_cff) %in% colnames(weird_chromosomes)]] <- NA
        weird_chromosomes$Metafusion_flag <- "Chromosome_not_in_bed"
        filtered_cff <- rbind(filtered_cff,weird_chromosomes)
    }
    write.table(
    filtered_cff[,c(header_filtered,"Metafusion_flag","cluster")],
    paste0(sample_name, "_metafusion_cluster.cff"),
    row.names = F,
    append = F,
    quote = F,
    sep = "\t"
    )
