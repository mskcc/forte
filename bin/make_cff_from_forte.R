#!/usr/local/bin/Rscript

library(dplyr)
library(data.table)
library(stringr)
cff_format <- c("chr1","pos1","strand1","chr2","pos2","strand2","library","sample_name",
                "sample_type","disease","tool",'split_cnt',"span_cnt","t_gene1","t_area1",
                "t_gene2","t_area2")


opt <- commandArgs(TRUE)

cff_format_df <- setNames(data.frame(matrix(ncol = length(cff_format), nrow = 0)), cff_format)

make_arriba <- function(sample_file){
df <- as.data.frame(matrix(ncol = 0,nrow=nrow(sample_file)))

df$t_gene1 <- sample_file$gene_id1
df$t_gene2 <- sample_file$gene_id2
df$chr1 <- str_split_fixed(sample_file$breakpoint1,":",2)[,1]
df$pos1 <- str_split_fixed(sample_file$breakpoint1,":",2)[,2]
df$chr2 <-  str_split_fixed(sample_file$breakpoint2,":",2)[,1]
df$pos2 <- str_split_fixed(sample_file$breakpoint2,":",2)[,2]
df$strand1 <- str_split_fixed(sample_file$`strand1(gene/fusion)`,"/",2)[,1]
df$strand2 <- str_split_fixed(sample_file$`strand2(gene/fusion)`,"/",2)[,1]
df$tool <- "arriba"
df$split_cnt <- ifelse(!is.na(sample_file$split_reads1), sample_file$split_reads1, -1)
df$span_cnt <- ifelse(!is.na(sample_file$discordant_mates), sample_file$discordant_mates, -1)
df$t_area1 <- sample_file$site1
df$t_area2 <- sample_file$site2

return(df)
}

make_fusioncatcher <- function(sample_file){
df <- as.data.frame(matrix(ncol = 0,nrow=nrow(sample_file)))

df$t_gene1 <- sample_file[,"Gene_1_id(5end_fusion_partner)"]
df$t_gene2 <- sample_file[,"Gene_2_id(3end_fusion_partner)"]
df$chr1 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"],":",3)[,1]
df$pos1 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"],":",3)[,2]
df$chr2 <-  str_split_fixed(sample_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"],":",3)[,1]
df$pos2 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"],":",3)[,2]
df$strand1 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"],":",3)[,3]
df$strand2 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"],":",3)[,3]
df$tool <- "fusioncatcher"
df$split_cnt <- ifelse(!is.na(sample_file$Spanning_unique_reads), sample_file$Spanning_unique_reads, -1)
df$span_cnt <- ifelse(!is.na(sample_file$Spanning_pairs), sample_file$Spanning_pairs, -1)
df$t_area1 <- sample_file$Predicted_effect
df$t_area2 <- sample_file$Predicted_effect


return(df)

}

make_starfusion <- function(sample_file){

df <- as.data.frame(matrix(ncol = 0,nrow=nrow(sample_file)))
df$t_gene1 <- str_split_fixed(sample_file$LeftGene,"\\^",2)[,2]
df$t_gene2 <- str_split_fixed(sample_file$RightGene,"\\^",2)[,2]
df$chr1 <- str_replace(str_split_fixed(sample_file$LeftBreakpoint,":",3)[,1],"chr","")
df$pos1 <- str_split_fixed(sample_file$LeftBreakpoint,":",3)[,2]
df$chr2 <-    str_replace(str_split_fixed(sample_file$RightBreakpoint,":",3)[,1],"chr","")
df$pos2 <-  str_split_fixed(sample_file$RightBreakpoint,":",3)[,2]
df$strand1 <- str_split_fixed(sample_file$LeftBreakpoint,":",3)[,3]
df$strand2 <- str_split_fixed(sample_file$RightBreakpoint,":",3)[,3]
df$tool <- "starfusion"
df$split_cnt <- ifelse(!is.na(sample_file$JunctionReadCount), sample_file$JunctionReadCount, -1)
df$span_cnt <- ifelse(!is.na(sample_file$SpanningFragCount), sample_file$SpanningFragCount, -1)
df$t_area1 <- sample_file$SpliceType
df$t_area2 <- sample_file$SpliceType

return(df)


}

sample_file <- fread(opt[2],data.table = F)

tool_cff <-  setNames(data.frame(matrix(ncol = length(cff_format), nrow = 0)), cff_format)
if(opt[1] == "arriba"){

if(nrow(sample_file) > 0){
    tool_cff <- make_arriba(sample_file)
    }
} else if( opt[1]  == "fusioncatcher") {

if(nrow(sample_file) > 0){
    tool_cff <-make_fusioncatcher(sample_file)
    }
} else if (opt[1]  == "starfusion"){

if(nrow(sample_file) > 0){
    tool_cff <-make_starfusion(sample_file)
    }
}
if(nrow(tool_cff) > 0 ){
tool_cff$sample_name <- opt[3]
tool_cff$library <- "RNA"
tool_cff$sample_type <- "Tumor"
tool_cff$disease <- NA
}

tool_cff$strand1[tool_cff$strand1 == "."] <- NA
tool_cff$strand2[tool_cff$strand2 == "."] <- NA
tool_cff <- tool_cff[,cff_format]
write.table(tool_cff,opt[4],sep = "\t",row.names = F,quote =  F,col.names =F)
