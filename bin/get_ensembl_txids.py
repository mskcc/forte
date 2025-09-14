#!/usr/bin/env python

import argparse
import sys
import time
from gtfparse import read_gtf
from time import strftime
from collections import defaultdict as dd

def main(argv):
	print('job starts', strftime('%a, %d %b %Y %I:%M:%S'))
	start_time = time.time()


	gene2_txids = dd(list)
	gtf = read_gtf(opts.g, result_type="pandas")

	gtf_tx = gtf[gtf["feature"] == "transcript"]

	# 	print(gene_name, tx_id, tx_version)

	target_genes = set()
	with open(opts.t) as f:
		for line in f:
			gene = line.strip().split('\t')[3].split('_')[0]
			target_genes.add(gene)


	with open(opts.o, 'w') as out:
		out.write('IsPanel\tDescription\tGene\tLookup_Transcript\tReported_Transcript\tAssay_Type\n')

		is_panel = '1'
		assay_type = "TARGET"

		for gene in target_genes:

			tx_df = gtf_tx[gtf_tx['gene_name'] == gene]
			for index, row in tx_df.iterrows():
				reported_tx_id = row['transcript_id']
				tx_version = row['transcript_version']

				lookup_tx_id = f"{reported_tx_id}.{tx_version}"

				#reported_tx_id = lookup_tx_id.split('.')[0]
				out.write(f'{is_panel}\tRv1\t{gene}\t{lookup_tx_id}\t{reported_tx_id}\t{assay_type}\n')



	print("--- %s seconds ---" % (time.time() - start_time))
	print('DONE!', strftime('%a, %d %b %Y %I:%M:%S'))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-t', metavar='targetfile', required=True, help='RNA-seq targets')
	parser.add_argument('-g', metavar='v11_gtf', required=True, help='Ensemblev111 GTF')
	parser.add_argument('-o', metavar='outf', required=True, help='Reference transcript IDs for Emsembl v111 genes')

	opts = parser.parse_args()

	main(sys.argv[1:])