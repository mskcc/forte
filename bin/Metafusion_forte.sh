#!/bin/bash
#STEPS
rename=1
annotate=1
annotate_exons=1
merge=1
output_ANC_RT_SG=1
RT_call_filter=1
blck_filter=1
ANC_filter=1
rank=1
fusionannotator=1
benchmark=1

# Loop through arguments and process them
while test $# -gt 0;do
    case $1 in
        -n=*|--num_tools=*)
        num_tools="${1#*=}"
        shift
        ;;
        --outdir)
        outdir="$2"
        shift 2
        ;;
        --cff)
        cff="$2"
        shift 2
        ;;
        --gene_bed)
        gene_bed="$2"
        shift 2
        ;;
        --gene_info)
        gene_info="$2"
        shift 2
        ;;
        --genome_fasta)
        genome_fasta="$2"
        shift 2
        ;;
        --truth_set)
        truth_set="$2"
        shift 2
        ;;
        --recurrent_bedpe)
        recurrent_bedpe="$2"
        shift 2
        ;;
        --scripts)
        fusiontools="$2"
        shift 2
        ;;
        --fusion_annotator)
        FA=1
        shift
        ;;
        --annotate_exons)
        exons=1
        shift
        ;;
        *)
        #OTHER_ARGUMENTS+=("$1")
        shift # Remove generic argument from processing
        ;;
    esac
done
#echo num_tools $num_tools
#echo outdir $outdir
#echo cff  $cff
#echo gene_bed $gene_bed
#echo gene_info $gene_info
#echo genome_fasta $genome_fasta
#echo truth_set $truth_set

#METAFUSION WORKFLOW
mkdir $outdir

#Check CFF file format:
#Remove entries with nonconformming chromosome name
#Remove "." from strand field and replace with "NA"
cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'BEGIN{FS=OFS="\t"} $3 !~ /^[-+]$/{$3="NA"} 1' | awk 'BEGIN{FS=OFS="\t"} $6 !~ /^[-+]$/{$6="NA"} 1'   > $outdir/$(basename $cff).reformat
#cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'FS=OFS="\t"{if ($3==".") {$3="NA"}; print}'|  awk 'FS=OFS="\t"{if ($6==".") {$6="NA"}; print}'  > $outdir/$(basename $cff).reformat
#cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' | awk 'FS=OFS="\t"{ if ($3 !="+" || $3 != "-" ); $3="NA"; print }' | awk 'FS=OFS="\t"{ if ($6 !="+" || $6 != "-" ); $6="NA"; print }' > $outdir/$(basename $cff).reformat
cff=$outdir/$(basename $cff).reformat
#NEED TO INSERT +/-/NA for strand, make NA if other


#Rename cff
if [ $rename -eq 1 ]; then
echo Rename cff
echo rename_cff_file_genes.MetaFusion.py $cff $gene_info \> $outdir/$(basename $cff).renamed
rename_cff_file_genes.MetaFusion.py $cff $gene_info > $outdir/$(basename $cff).renamed
cff=$outdir/$(basename $cff).renamed
fi

#Annotate cff
if [ $annotate -eq 1 ]; then
if [ $genome_fasta ]; then
    echo Annotate cff, extract sequence surrounding breakpoint
    reann_cff_fusion.py --cff $cff --gene_bed $gene_bed --ref_fa $genome_fasta > $outdir/$(basename $cff).reann.WITH_SEQ
else
    echo Annotate cff, no extraction of sequence surrounding breakpoint
    reann_cff_fusion.py --cff $cff --gene_bed $gene_bed > $outdir/$(basename $cff).reann.NO_SEQ
fi
fi
# Assign .cff based on SEQ or NOSEQ
if [ $genome_fasta ]; then
cff=$outdir/$(basename $cff).reann.WITH_SEQ
echo cff $cff
else
cff=$outdir/$(basename $cff).reann.NO_SEQ
echo cff $cff
fi

if [ $annotate_exons -eq 1 ] && [ $exons -eq 1 ]; then
echo Add adjacent exons to cff
extract_closest_exons.py $cff $gene_bed $genome_fasta  > $outdir/$(basename $cff).exons
fi
# assign cff as ".exons" if --annotate_exons flag was specified
if [ $exons -eq 1 ]; then
cff=$outdir/$(basename $cff).exons
fi

#Merge
cluster=$outdir/$(basename $cff).cluster
if [ $exons -eq 1 ] && [ $merge -eq 1 ]; then
echo Merge cff by exons
RUN_cluster_exons.sh $cff $outdir > $cluster
elif [ $merge -eq 1 ]; then
echo Merge cff by genes and breakpoints
RUN_cluster_genes_breakpoints.sh $cff $outdir > $cluster
fi

#output ANC_RT_SG file
if [ $output_ANC_RT_SG -eq 1 ]; then
echo output cis-sage.cluster file
output_ANC_RT_SG.py $cluster > $outdir/cis-sage.cluster
fi

#ReadThrough Callerfilter
if [ $RT_call_filter -eq 1 ]; then
echo ReadThrough, callerfilter $num_tools
cat $cluster | grep ReadThrough > $outdir/$(basename $cluster).ReadThrough
echo python callerfilter_num.py --cluster $cluster  --num_tools $num_tools \| grep -v ReadThrough \> $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
callerfilter_num.py --cluster $cluster  --num_tools $num_tools | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
fi
cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools

# Blocklist Filter
if [ $blck_filter -eq 1 ]; then
echo blocklist filter
#echo bash blocklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe
blocklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
fi
cluster=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter


# Adjacent Noncoding filter
if [ $ANC_filter -eq 1 ]; then
echo ANC adjacent noncoding filter
filter_adjacent_noncoding.py $cluster > $outdir/$(basename $cluster).ANC_filter
fi
cluster=$outdir/$(basename $cluster).ANC_filter

#Rank and generate final.cluster
if [ $rank -eq 1 ]; then
echo Rank and generate final.cluster
rank_cluster_file.py $cluster > $outdir/final.n$num_tools.cluster
fi
cluster=$outdir/final.n$num_tools.cluster
### Generate filtered FID file
out=`awk -F '\t' '{print $15}' $cluster  | tail -n +2`
for this in echo ${out//,/ }; do grep $this $cff; done >> $outdir/$(basename $cff).filtered.cff

### TURNINGO FF FUSION ANNOATOR AND BENCHMARK< DO THIS LATER
#fusionannotator
#if [ $fusionannotator -eq 1 ] && [ $FA -eq 1 ]; then
#  echo Running FusionAnnotator
#  bash RUN_FusionAnnotator.sh $outdir $cluster
#  echo Adding FusionAnnotator database hits to final.cluster.CANCER_FUSIONS file
#  FA_db_file=$outdir/cluster.preds.collected.gencode_mapped.wAnnot.CANCER_FUSIONS
#  python add_db_hits_to_cluster.py $cluster $FA_db_file > $outdir/$(basename $cluster).CANCER_FUSIONS
#  cluster=$outdir/$(basename $cluster).CANCER_FUSIONS
#fi

#Benchmark
#if [ $benchmark -eq 1 ] && [ $truth_set ]; then
# echo Running benchmarking
#  bash benchmarking_cluster.MetaFusion.sh $outdir $truth_set $cluster
#else
#	echo no benchmarking performed
#fi


