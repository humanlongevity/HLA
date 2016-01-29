#!/usr/bin/env python
import json
from glob import glob
from datetime import datetime
import argparse
import re

version = '1.0'

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-in', dest = 'input', type = file, required = True, help = 'Input file')
	parser.add_argument('-out', required = True, help = 'Output JSON')
	parser.add_argument('-subject', required = True, help = 'Subject ID')
	parser.add_argument('-sample', required = True, help = 'Sample ID')
	args = parser.parse_args()

	creation_time = datetime.utcnow().isoformat()
	creation_time = creation_time[0:creation_time.index('.')] + 'Z'

	suffix = re.compile(r'\D+$')

	hla = set()
	for line in args.input:
		data = line.split('\t')
		if data[0] == '1':
			allele = re.sub(suffix, '', data[1])
			if 'DPB' not in allele:
				hla.add(allele)

	traits = {
		'subject_id': args.subject, 
		'sample_id': args.sample, 
		'report_type': 'hla',
		'report_version': version,
		'creation_time': creation_time, 
		'hla': {
			'alleles': list(hla)
		}
		
	}
	json.dump(traits, open(args.out, 'w'), indent = True)


if __name__ == '__main__':
	main()
