#!/bin/bash
#STEPS

# __author__      = "Alexandria Dymun"
# __email__       = "pintoa1@mskcc.org"
# __contributor__ = "Anne Marie Noronha (noronhaa@mskcc.org)"
# __version__     = "0.0.1"
# __status__      = "Dev"


output_ANC_RT_SG=1
RT_call_filter=1
blck_filter=1
ANC_filter=1
usage() {
    echo "Usage: Metafusion_forte.sh [--num_tools=<minNumToolsCalled> --genome_fasta <FASTA adds SEQ to fusion>  --recurrent_bedpe <blacklistFusions> ] --outdir <outputDirectory> --cff <cffFile> --gene_bed <geneBedFile> --gene_info <geneInfoFile>" 1>&2;
    exit 1;
}

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
        --recurrent_bedpe)
        recurrent_bedpe="$2"
        shift 2
        ;;
        *)
        #OTHER_ARGUMENTS+=("$1")
        shift # Remove generic argument from processing
        ;;
    esac
done

if [[ ! $cff || ! $gene_info  || ! $gene_bed ]]; then
    echo "Missing required argument"
    usage
fi


mkdir $outdir

#Check CFF file format:
#Remove entries with nonconformming chromosome name
#Remove "." from strand field and replace with "NA"
cat $cff | awk '$1 ~ /[0-9XY]/ && $4 ~ /[0-9XY]/ ' |  awk 'BEGIN{FS=OFS="\t"} $3 !~ /^[-+]$/{$3="NA"} 1' | awk 'BEGIN{FS=OFS="\t"} $6 !~ /^[-+]$/{$6="NA"} 1'   > $outdir/$(basename $cff).reformat
cff=$outdir/$(basename $cff).reformat


#Rename cff
echo Rename cff
rename_cff_file_genes.MetaFusion.py $cff $gene_info > $outdir/$(basename $cff).renamed
cff=$outdir/$(basename $cff).renamed

### remove any chromosomes that do not exist in genebed
all_gene_bed_chrs=`awk -F '\t' '{print $1}' $gene_bed | sort | uniq `
awk -F " " -v arr="${all_gene_bed_chrs[*]}" 'BEGIN{OFS = "\t"; split(arr,arr1); for(i in arr1) dict[arr1[i]]=""} $1 in dict && $4 in dict' $cff >  $outdir/$(basename $cff).cleaned_chr
grep -v -f $outdir/$(basename $cff).cleaned_chr $cff > problematic_chromosomes.cff
cff=$outdir/$(basename $cff).cleaned_chr

#Annotate cff
if [ $genome_fasta ]; then
    echo Annotate cff, extract sequence surrounding breakpoint
    reann_cff_fusion.py --cff $cff --gene_bed $gene_bed --ref_fa $genome_fasta > $outdir/$(basename $cff).reann.WITH_SEQ
else
    echo Annotate cff, no extraction of sequence surrounding breakpoint
    reann_cff_fusion.py --cff $cff --gene_bed $gene_bed > $outdir/$(basename $cff).reann.NO_SEQ
fi

# Assign .cff based on SEQ or NOSEQ
if [ $genome_fasta ]; then
    cff=$outdir/$(basename $cff).reann.WITH_SEQ
    echo cff $cff
else
    cff=$outdir/$(basename $cff).reann.NO_SEQ
    echo cff $cff
fi

echo Add adjacent exons to cff
extract_closest_exons.py $cff $gene_bed $genome_fasta  > $outdir/$(basename $cff).exons

# assign cff as ".exons" if --annotate_exons flag was specified

cff=$outdir/$(basename $cff).exons


#Merge
cluster=$outdir/$(basename $cff).cluster
echo Merge cff by genes and breakpoints
RUN_cluster_genes_breakpoints.sh $cff $outdir > $cluster

#output ANC_RT_SG file
if [ $output_ANC_RT_SG -eq 1 ]; then
    echo output cis-sage.cluster file
    output_ANC_RT_SG.py $cluster > $outdir/cis-sage.cluster
fi

#ReadThrough Callerfilter
if [ $RT_call_filter -eq 1 ]; then
    echo ReadThrough, callerfilter $num_tools
    cat $cluster | grep ReadThrough > $outdir/$(basename $cluster).ReadThrough
    callerfilter_num.py --cluster $cluster  --num_tools $num_tools | grep -v ReadThrough  > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
    cluster_RT_call=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools
fi
# Blocklist Filter
if [ $recurrent_bedpe ]; then
    echo blocklist filter
    blocklist_filter_recurrent_breakpoints.sh $cff $cluster_RT_call $outdir $recurrent_bedpe > $outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter

    cluster=$outdir/$(basename $cluster).RT_filter.callerfilter.$num_tools.blck_filter
fi
# Adjacent Noncoding filter
if [ $ANC_filter -eq 1 ]; then
    echo ANC adjacent noncoding filter
    filter_adjacent_noncoding.py $cluster > $outdir/$(basename $cluster).ANC_filter

    cluster=$outdir/$(basename $cluster).ANC_filter
fi
#Rank and generate final.cluster
echo Rank and generate final.cluster
rank_cluster_file.py $cluster > $outdir/final.n$num_tools.cluster
cluster=$outdir/final.n$num_tools.cluster
### Generate filtered FID file
out=`awk -F '\t' '{print $15}' $cluster  | tail -n +2`
out2=`awk -F '\t' '{print $22}' $outdir/cis-sage.cluster | tail -n +2`
out3=`echo $out $out2`

for this in echo ${out3//,/ }; do grep $this $cff; done >> $outdir/$(basename $cff).filtered.cff