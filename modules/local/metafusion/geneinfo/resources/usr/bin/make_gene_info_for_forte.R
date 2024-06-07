#!/usr/local/bin/Rscript
# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __contributor__ = "Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "0.0.1"
# __status__      = "Dev"


suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
})

usage <- function() {
    message("Usage:")
    message("make_gene_info_for_forte.R --primary_gtf <file.gtf> --fc_custom_bed_gene_names <custom_genes.bed> --star_fusion_ref <ref_annot.gtf> --fusioncatcher_ref <organism.gtf> --outputDir <outputpath>")
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

opt <- parse_args(paste(args,collapse=" "))

required_args <- c(
    "primary_gtf",
    "fc_custom_bed_gene_names",
    "star_fusion_ref",
    "fusioncatcher_ref",
    "outputDir"
)
if (length(setdiff(required_args,names(opt))) > 0) {
    message("Missing required arguments")
    usage()
    quit()
}

### primary gtf is v75, also used in arriba
primary_gtf <-  as.data.frame(rtracklayer::import(opt$primary_gtf))

### Fusion catcher has custom gene names/gene_ids....
## https://github.com/ndaniel/fusioncatcher/blob/ebc46fd1a8046fc909a56e09944a2ec2d69cc808/bin/add_custom_gene.py#L704-L715
fc_custom_bed_gene_names <- read.table(opt$fc_custom_bed_gene_names)
fc_custom_bed_gene_names$gene_name <- str_split_fixed(fc_custom_bed_gene_names$V4,"-",n=2)[,1]
fc_custom_bed_gene_names$gene_id <- str_split_fixed(fc_custom_bed_gene_names$V4,"-",n=3)[,2]
star_fusion_ref <- as.data.frame(rtracklayer::import(opt$star_fusion_ref))
fusioncatcher_ref <- as.data.frame(rtracklayer::import(opt$fusioncatcher_ref))


all_my_gene_ids_and_names <- list(primary_gtf,fc_custom_bed_gene_names,star_fusion_ref,fusioncatcher_ref)
### whichever gtf you label as primary should also be the same version of the gene_bed file generated for metafusion
names(all_my_gene_ids_and_names) <- c("primary","one","two","three")
unique_id_to_names <- lapply( all_my_gene_ids_and_names,function(gtf) {
    ### If gene id has versions, strip them off
    if(all(grepl("\\.",gtf$gene_id))){
        gtf$gene_id_with_version <- gtf$gene_id
        gtf$gene_id <- str_split_fixed(gtf$gene_id_with_version,"\\.",n=2)[,1]
        return(unique(gtf[,c("gene_name","gene_id","gene_id_with_version")]))
    } else{
        return(unique(gtf[,c("gene_name","gene_id")]))
    }
})


gene_info <- unique_id_to_names$primary
### tack on missing gene_ids from other references to gene info
versioned_gtf <-unlist(sapply(names(unique_id_to_names)[names(unique_id_to_names) != "primary"],function(name){
    if(any(colnames(unique_id_to_names[[name]]) == "gene_id_with_version")){
    return(name)
    }
}))


add_these_excess_gene_ids <- do.call(rbind,lapply(names(unique_id_to_names)[names(unique_id_to_names) != "primary"],function(name){
    add_symbols_and_ids <- unique_id_to_names[[name]]
    add_symbols_and_ids <- add_symbols_and_ids[!add_symbols_and_ids$gene_id %in% gene_info$gene_id,]
    if(name %in% versioned_gtf){
    add_symbols_and_ids <-add_symbols_and_ids[,c("gene_name","gene_id_with_version")]
    colnames(add_symbols_and_ids) <- c("gene_name","gene_id")
    }
    return(add_symbols_and_ids)

}))
# Excess genes being added (genes will be flagged as gene not in v75)
gene_info <- rbind(gene_info,add_these_excess_gene_ids)

gene_info <- merge(gene_info,do.call(rbind,unique_id_to_names[versioned_gtf])[,c("gene_id","gene_id_with_version")],by = "gene_id",all.x = T, all.y = F)

gene_info$Synonyms <- ifelse(is.na(gene_info$gene_id_with_version),gene_info$gene_id,paste0(gene_info$gene_id,"|",gene_info$gene_id_with_version))
gene_info$Symbol <- gene_info$gene_name

gene_info <- gene_info[,c("Symbol","Synonyms")]

write.table(
    gene_info,
    paste0(opt$outputDir,"/gene.info"),
    sep ="\t",
    quote = F,
    row.names = F

)
write.table(
    add_these_excess_gene_ids,
    paste0(opt$outputDir,"/excess_gene_ids.txt"),
    sep ="\t",
    quote = F,
    row.names = F
)
