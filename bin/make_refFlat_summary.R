#!/home/amoahk1/miniconda3/bin/Rscript
# __author__      = "Kofi Amoah"
# __email__       = "amoahk1@mskcc.org"

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

reflat <- read.delim(args[1], header=F, sep="\t")

data<-reflat %>%
	transform(V7 = as.numeric(gsub(pattern = "exon", replacement = "", V7))) %>%
		group_by(V5,V6) %>%
		mutate(first_exon = min(V7), last_exon = max(V7),
						pos1 = min(V2,V3), pos2 = max(V2,V3)) %>%
		select(V1, V4, V5, V6, first_exon, last_exon, pos1, pos2) %>%
		distinct()

colnames(data) <- c("Chrom", "Strand", "Gene", "Transcript", "first_exon", "last_exon", "pos1", 'pos2')

write.table(data, file=args[2],
	col.names = T, sep = "\t", row.names = F, quote = F)

