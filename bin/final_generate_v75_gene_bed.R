
# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __contributor__ = "Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "0.0.1"
# __status__      = "Dev"



library(dplyr)
library(data.table)
library(stringr)
gtf <- rtracklayer::import('/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf')
gtf_df <- as.data.frame(gtf)


# Utilized gtf from igenomes for FORTE This corresponds to GRCh37 ensembl 75
# Add introns to gtf, convert to gff3
# bsub -R "rusage[mem=64]" -o add_introns_agat_%J.out singularity exec -B /juno/ \\
# -B /tmp -B /scratch/ docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0  \\
# /bin/bash -c "agat_sp_add_introns.pl -g /juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf\\
# -o genes.INTRONS.gff3"
# gff2bed < genes.INTRONS.gff3 > genes.INTRONS.agat.bed


total.introns.bed <- fread(file="/work/ccs/pintoa1/metafusion_refs/meta_fusion_bed_generation/genes.INTRONS.agat.bed", header = FALSE, stringsAsFactors = F, sep="\t", na.strings = "",data.table = F)
colnames(total.introns.bed) <- c("chr","start","end","gene_id","tmp","strand","gene_biotype","type","V9","description")
total.introns.bed$transcript_id <- gsub("\\;.*","",str_split_fixed(total.introns.bed$description,"transcript_id=",n=2)[,2])
total.introns.bed$gene_name <-gsub("\\;.*","",str_split_fixed(total.introns.bed$description,"gene_name=",n=2)[,2])

transcript_ids <- unique(total.introns.bed$transcript_id)
file.to_write <- "/work/ccs/pintoa1/metafusion_refs/meta_fusion_bed_generation/cleaned_metafusion_v75_gene.bed"

if(file.exists(file.to_write) ) {file.remove(file.to_write)}

#START CLOCK: THE INDEXING TAKES A LONG TIME, LIKE 5 HOURS
ptm <- proc.time()

# Index each transcript feature, incrementing when an intron is passed
## metafusion expects exon count 0 to (N(exons)-1)
## Forward strand: Exon 0 == Exon 1
### Reverse strand: Exon 0 == LAST EXON IN TRANSCRIPT
for (id in transcript_ids){
    transcript <- total.introns.bed[total.introns.bed$transcript_id == id,]
    # Remove exons if coding gene, since "exon" and "CDS" are duplicates of one another
    if ("CDS" %in% transcript$type){
        transcript <- transcript[!transcript$type == "exon",]
        }
    # Order features by increasing bp
    transcript <- transcript[order(transcript$start, decreasing = FALSE),]
    # Index features
    idx <- 0
    for (i in 1:nrow(transcript)){
        transcript$idx [i]<- idx
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
    if ("CDS" %in% unique(transcript$type)){transcript[transcript$type == "CDS",]$type <- "cds"}
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
            start_coding <- max(transcript[transcript$type == "cds","end"])
            stop_coding <- min(transcript[transcript$type == "cds","start"])
            transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr3"
            transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr5"
        }
    }
    transcript <- transcript[,c("chr", "start", "end", "transcript_id", "type", "idx", "strand", "gene_name", "gene_id" )]
    write.table(transcript, file.to_write, append=TRUE, sep="\t", quote=F,  row.names=F, col.names=F)
}

time <- proc.time() - ptm
time
#
# user    system   elapsed
# 16657.116    32.227 16741.382


new.bed <- fread(file.to_write,data.table = F)
colnames(new.bed) <- c("chr","start","end","transcript_id","type","idx","strand","gene_name","gene_id")

#### Any exon that remains after teh cds change, is likely and untranslated region. change below

# Basically, subfeatures which are "exon" need to be changed (i.e. exon --> utr3/utr5)
#Forward strand
new.bed$type[new.bed$strand == "f" &  new.bed$type == "exon" ] <- "utr5"
#Reverse strand
new.bed$type[new.bed$strand == "r" &  new.bed$type == "exon"]<- "utr3"

expected_types <- c("cds","intron","utr3","utr5")
new.bed.ready <- new.bed[new.bed$type %in% c(expected_types),]

write.table(new.bed.ready, "/work/ccs/pintoa1/metafusion_refs/meta_fusion_bed_generation/v75_gene.bed",  sep="\t", quote=F,  row.names=F, col.names=F)
