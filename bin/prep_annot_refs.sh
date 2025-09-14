#!/bin/bash

set -eo pipefail

# __author__      = "Kofi Amoah"
# __email__       = "amoahk1@mskcc.org"

usage() {
	echo "Usage: prep_annot_refs.sh \
			--target_genes <xgen panel genes> \
			--gtf <Ensemblv111 gtf gzipped> \
			--ref_txs <Reference transcripts ID file> \
			--svtable <sv_table for iAnnotateSV> \
			--refFlat <refFLAT> \
			--ref_summary <refFLAT summary> \
			--kinase_domains <kinase domains>" 1>&2;
	exit 1;
}

# Process args
while [[ "$#" -gt 0 ]];
do
	case "$1" in
		-t|--target_genes)
			target_genes="$2"
			shift 2
			;;
		-g|--gtf)
			gtf="$2"
			shift 2
			;;
		-p|--pfam)
			pfam="$2"
			shift 2
			;;
		-r|--ref_txs)
			ref_txs="$2"
			shift 2
			;;
		-s|--svtable)
			svtable="$2"
			shift 2
			;;
		-f|--refFlat)
			refFlat="$2"
			shift 2
			;;
		-u|--ref_summary)
			refFlat_summary="$2"
			shift 2
			;;
		-k|--kinase_domain)
			kinase_domain="$2"
			shift 2
			;;
	esac
done

if [[ ! $target_genes || ! $gtf ]]; then
	echo "Missing required arguments"
fi


### Get Ensembl canonical tx ids
get_ensembl_txids.py -t $target_genes -g $gtf -o $ref_txs
echo "Making Ref_TxID file....DONE!"

### Make sv_table for iAnnotate_sv
echo -e "#name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\texonFrames" > $svtable
gtfToGenePred -genePredExt -geneNameAsName2 $gtf tmp
less tmp | awk -F'\t' '{print $1, "chr"$2, $3, $4, $5, $6, $7, $8, $9, $10, $12, $NF}' OFS="\t"  >>  $svtable
rm tmp
echo "Making SV_table for iAnnotateSV....DONE!"

### Make Ensembl version of refFlat and summary files
make_refFlat.py -i $svtable -t $ref_txs -o $refFlat
make_refFlat_summary.R $refFlat $refFlat_summary
echo "Making refFlat and refFlat_summary files....DONE!"

### Reformat kinase domain data to include gene and transcript ids
less $pfam | grep -P "Pkinase\t" | grep -v _ | cut -f2-7 > raw_pkinase.txt
reformat_kinase_domains.py -i raw_pkinase.txt -r $refFlat_summary -o $kinase_domain
rm $pfam raw_pkinase.txt
echo "Making kinase domain table....DONE"
