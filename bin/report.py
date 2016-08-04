#!/usr/bin/env python
# Author: Xie Chao
import json
from glob import glob
from datetime import datetime
import argparse
import re

version = '1.2'

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-in', dest = 'input', type = file, required = True, help = 'Input file')
	parser.add_argument('-out', required = True, help = 'Output JSON')
	parser.add_argument('-subject', required = True, help = 'Subject ID')
	parser.add_argument('-sample', required = True, help = 'Sample ID')
	args = parser.parse_args()

	creation_time = datetime.utcnow().isoformat()
	creation_time = creation_time[0:creation_time.index('.')] + 'Z'

#	suffix = re.compile(r'\D+$')

	hla = list()
	count = {}
	back = {}
	for line in args.input:
		data = line.split('\t')
		if data[0] == '0':
#			allele = re.sub(suffix, '', data[1])
			allele = data[1]
			gene = allele[0:allele.index('*')]
			back[gene] = allele
			if gene in count:
				count[gene] += 1
			else:
				count[gene] = 1
			hla.append(allele)

	for gene in count:
		if count[gene] == 1:
			hla.append(back[gene])

	hla.sort()

	traits = {
		'subject_id': args.subject, 
		'sample_id': args.sample, 
		'report_type': 'hla_typing',
		'report_version': version,
		'creation_time': creation_time, 
		'hla': {
			'alleles': hla
		}
		
	}
	json.dump(traits, open(args.out, 'w'), indent = True)


if __name__ == '__main__':
	main()
