#!/usr/bin/env python

import sys
import argparse
import time
from time import strftime
from collections import defaultdict as dd


def main(argv):
	print('job starts', strftime('%a, %d %b %Y %I:%M:%S'))
	start_time = time.time()

	raw_pkinase = dd(set)
	with open(opts.i) as f:
		for line in f:
			chrm, kstart, kend, typ, zero, strd = line.strip().split('\t')
			chrm = chrm[3:]
			kstart, kend = int(kstart), int(kend)
			raw_pkinase[(chrm, strd)].add((kstart, kend))

	refFlat = dd(set)
	with open(opts.r) as f:
		next(f)
		for line in f:
			chrm, strd, gene, tx, first_exon, last_exon, cstart, cend = line.strip().split('\t')
			cstart, cend = int(cstart), int(cend)
			refFlat[(chrm, strd)].add((cstart, cend, gene, tx))


	with open(opts.o, 'w') as out:
		out.write('Chr\tStart\tEnd\tEnsembl_transcriptID\tGene\n')
		for (chrm, strd) in raw_pkinase:
			for (kstart, kend) in raw_pkinase[(chrm, strd)]:
				for (cstart, cend, gene, tx) in refFlat[(chrm, strd)]:
					if kstart >= cstart and kstart <= cend and kend >= cstart and kend <= cend:
						out.write(f"{chrm}\t{kstart}\t{kend}\t{tx}\t{gene}\n")


	print("--- %s seconds ---" % (time.time() - start_time))
	print('DONE!', strftime('%a, %d %b %Y %I:%M:%S'))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', metavar='inf', required=True, help='Raw Pkinase from UCSC')
	parser.add_argument('-r', metavar='ref', required=True, help='Transcript reference IDs')
	parser.add_argument('-o', metavar='outf', required=True, help='Kinase domains')

	opts = parser.parse_args()
	print("Inf:", opts.i)
	print("Ref:", opts.r)
	print("Outf:", opts.o)

	main(sys.argv[1:])