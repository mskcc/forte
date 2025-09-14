#!/usr/bin/env python

import argparse
import sys
import time
import numpy as np
from time import strftime
from collections import defaultdict as dd

def main(argv):
	print('job starts', strftime('%a, %d %b %Y %I:%M:%S'))
	start_time = time.time()

	ref_txids = dd(str)
	with open(opts.t) as f:
		next(f)
		for line in f:
			ispanel, desc, gene, lookup_tx, reported_tx, assay_type = line.strip().split('\t')
			ref_txids[(gene, reported_tx)] = lookup_tx

	refFlat_data = dd(list)
	with open(opts.i) as f:
		next(f)
		for line in f:
			tx, chrm, strd, tx_start, tx_end, cds_start, cds_end, exon_count, exon_starts, exon_ends, gene, exon_frames = line.strip().split('\t')

			if not (gene, tx) in ref_txids: continue

			tx = ref_txids[(gene, tx)] #use the tx with dot version
			chrm = chrm[3:]
			tx_start = int(tx_start)
			tx_end = int(tx_end)
			exon_count = int(exon_count)
			exon_starts = np.array(exon_starts.split(',')[:-1]).astype(np.int32)
			exon_ends = np.array(exon_ends.split(',')[:-1]).astype(np.int32)
			exon_num = 1
			if strd == '-':
				exon_num = exon_count

			i = 0
			while i < exon_count:
				e_start, e_end = exon_starts[i], exon_ends[i]
				refFlat_data[(gene, tx)].append((chrm, e_start, e_end, strd, f"exon{exon_num}"))
				if strd == '+':
					exon_num += 1
				else:
					exon_num -= 1

				i += 1 


	with open(opts.o, 'w') as out:
		for (gene, tx) in refFlat_data:
			for (chrm, e_start, e_end, strd, exon_num) in refFlat_data[(gene, tx)]:
				out.write(f"{chrm}\t{e_start}\t{e_end}\t{strd}\t{gene}\t{tx}\t{exon_num}\n")


	print("--- %s seconds ---" % (time.time() - start_time))
	print('DONE!', strftime('%a, %d %b %Y %I:%M:%S'))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', metavar='sv-table', required=True, help='sv-table')
	parser.add_argument('-t', metavar='tx_refs', required=True, help='Reference txids')
	parser.add_argument('-o', metavar='outf', required=True, help='Custom refFlat file')

	opts = parser.parse_args()

	main(sys.argv[1:])