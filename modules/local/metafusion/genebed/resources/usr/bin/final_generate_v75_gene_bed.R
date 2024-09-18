#!/usr/local/bin/Rscript

# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __contributor__ = "Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "0.0.2"
# __status__      = "Dev"


suppressPackageStartupMessages({
    library(plyr)
    library(dplyr)
    library(data.table)
    library(stringr)
    options(scipen = 999)
})

usage <- function() {
    message("Usage:")
    message("final_generate_v75_gene_bed.R <in.gff> <out.bed>")
}

args = commandArgs(TRUE)

if (length(args)!=2) {
    usage()
    quit()
}

# Utilized gtf from igenomes for FORTE This corresponds to GRCh37 ensembl 75
# Add introns to gtf, convert to gff3
# bsub -R "rusage[mem=64]" -o add_introns_agat_%J.out singularity exec -B /juno/ \\
# -B /tmp -B /scratch/ docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0  \\
# /bin/bash -c "agat_sp_add_introns.pl -g /juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf\\
# -o genes.INTRONS.gff3"

gtf <- rtracklayer::import(args[1])
gtf_df <- as.data.frame(gtf)
#remove incomplete transcripts mRNA_end_NF and mRNA_start_NF (not finished)
gtf_df <- gtf_df[!grepl("NF",gtf_df$tag),]

file.to_write <- args[2]

### ensure start is 0 based
gtf_df <- gtf_df %>%
    rename(
        chr = seqnames
    ) %>%
    select(c(chr, start, end, transcript_id, type, strand, gene_name, gene_id)) %>%
    filter(type %in% c("exon","intron","UTR","CDS","cds","utr")) %>% mutate(start = start-1)


#START CLOCK
ptm <- proc.time()
print(ptm)

# Index each transcript feature, incrementing when an intron is passed
## metafusion expects exon count 0 to (N(exons)-1)
## Forward strand: Exon 0 == Exon 1
### Reverse strand: Exon 0 == LAST EXON IN TRANSCRIPT

print(dim(gtf_df))
print(length(unique(gtf_df$transcript_id)))

modify_transcript <- function(transcript){

    # Remove exons if coding gene, since "exon" and "CDS" are duplicates of one another
    if ("CDS" %in% transcript$type){
        transcript <- transcript[!transcript$type == "exon",]
        }
    # Order features by increasing bp
    transcript <- transcript[order(transcript$start, decreasing = FALSE),]
    # Index features
    idx <- 0
    for (i in 1:nrow(transcript)){
        transcript$idx[i]<- idx
        if (transcript$type[i] == "intron"){
            idx <- idx + 1
        }
    }
    # REFORMAT TRANSCRIPT
    #Change strand info (+ --> f, - --> r)
    if (unique(transcript$strand) == "+"){
        transcript$strand <- 'f'
    } else if  (unique(transcript$strand) == "-"){
        transcript$strand <- 'r'
    } else {
        errorCondition("Strand info for this transcript is inconsistent")
    }
    #Add "chr" prefix to chromosomes
    transcript$chr <- sapply("chr", paste0,  transcript$chr)
    #Change CDS --> cds ### IF A TRANSCRIPT LACKS "CDS" THIS LINE WILL DO NOTHING, Changing exon values to UTRs later
    transcript <- transcript %>% mutate(type = as.character(type))
    transcript <- transcript %>% mutate(type=ifelse(type == "CDS","cds",type))
    ## DETERMING UTR3 and UTR5
    ### INSTEAD OF START AND STOP, USE CDS LOCATIONS AND STRAND INFORMATION.....
    if ("UTR" %in% unique(transcript$type)){
        if( unique(transcript$strand) == "f"){
            #Forward strand
            start_coding <- min(transcript[transcript$type == "cds","start"])
            stop_coding <-  max(transcript[transcript$type == "cds","end"])
            transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr5"
            transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr3"
        }else {
            #Reverse strand
            start_coding <- max(transcript[transcript$type == "cds","end"])
            stop_coding <- min(transcript[transcript$type == "cds","start"])
            transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr3"
            transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr5"
        }
    }
    #### Any exon that remains after teh cds change, is likely and untranslated region. change below

    # Basically, subfeatures which are "exon" need to be changed (i.e. exon --> utr3/utr5)
    #Forward strand
    transcript$type[transcript$strand == "f" &  transcript$type == "exon" ] <- "utr5"
    #Reverse strand
    transcript$type[transcript$strand == "r" &  transcript$type == "exon"]<- "utr3"
    #transcript <- transcript[,c("chr", "start", "end", "transcript_id", "type", "idx", "strand", "gene_name", "gene_id" )]
    expected_types <- c("cds","intron","utr3","utr5")
    transcript <- transcript[transcript$type %in% c(expected_types),]
    return(transcript)
}

if(file.exists(file.to_write) ) {file.remove(file.to_write)}

gtf_df_modified <- gtf_df %>%
    group_by(transcript_id,.drop = FALSE) %>%
    group_modify(~ modify_transcript(.x)) %>%
    select(c(chr, start, end, transcript_id, type, idx, strand, gene_name, gene_id )) %>%
    arrange(chr,start,end)

time <- proc.time() - ptm
print(time)

write.table(
    gtf_df_modified,
    file.to_write,
    sep="\t",
    quote=F,
    row.names=F,
    col.names=F
)
