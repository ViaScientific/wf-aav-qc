#!/usr/bin/env python3

import argparse, HTSeq, sys, gzip

def main(args):

	upstream_itr_end, downstream_itr_start = process_annotation(args)
	target_reference = args.target

	read_key = {}
	
	with gzip.open(args.per_read, 'rt') as infile:
		header = infile.readline()
		for line in infile:
			cur = line.rstrip().split('\t')
			read_key[cur[0]] = (cur[3], cur[4], cur[6])
	
	print("reference\tstart\tend\talignment\talignment_quality\tstrand\ttype\tsubtype\tlabel\tcontinuous\tnon-continuous\timpurity")

	with HTSeq.BAM_Reader(args.tagged_bam) as bam:
		
		for alignment in bam:
			
			if alignment.aligned:
	
				if not alignment.supplementary:
					
					reference = alignment.iv.chrom
					start = alignment.iv.start
					end = alignment.iv.end
					strand = alignment.iv.strand
						
					alignment_length = end - start
	
					continuous = 0
					noncontinuous = 0
					impurity = 0
	
					if read_key[alignment.read.name][0] == 'ssAAV':
						
						if start < upstream_itr_end and end > downstream_itr_start:
							continuous = downstream_itr_start - upstream_itr_end
						elif start < upstream_itr_end and end > upstream_itr_end:
							noncontinuous = end - upstream_itr_end
						elif end > downstream_itr_start and start < downstream_itr_start:
							noncontinuous = downstream_itr_start - start
						else:
							noncontinuous = alignment_length
					
					elif read_key[alignment.read.name][0] == 'backbone':
						impurity = alignment_length
	
			
					elif read_key[alignment.read.name][0] == 'NA':
						
						if read_key[alignment.read.name][2] == 'repcap':
							impurity = alignment_length
	
						elif read_key[alignment.read.name][2] == 'host':
							impurity = alignment_length
						elif read_key[alignment.read.name][2] == 'helper':
							impurity = alignment_length
						elif read_key[alignment.read.name][2] == 'chimeric-vector' or read_key[alignment.read.name][2] == 'chimeric-nonvector':
							result = process_alignment(target_reference, reference, alignment.cigar, alignment.optional_fields)
							noncontinuous = result[0]
							impurity = result[1]
	
			
					elif read_key[alignment.read.name][0] == 'other-vector':
						
						## Same rules for complex, snapback, tandem, unclassified, unresolved-dimer
						## If main alignment section spans the ITRs then call continous + noncontinous for auxilary sections.
						## If main alignment section doesn't span ITRs then all are noncontinous.
	
						if start < upstream_itr_end and end > downstream_itr_start:
							continuous = downstream_itr_start - upstream_itr_end
							result = process_alignment(target_reference, reference, None, alignment.optional_fields)
						else:
							result = process_alignment(target_reference, reference, alignment.cigar, alignment.optional_fields)
						noncontinuous = result[0]
						impurity = result[1]
		
					print('%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d' % (reference, start, end, alignment.read.name, alignment.aQual, strand, read_key[alignment.read.name][0], read_key[alignment.read.name][1], read_key[alignment.read.name][2], continuous, noncontinuous, impurity))

def process_annotation(args):

	if args.annotation != None and args.upstream_label != None and args.downstream_label != None:

		upstream_candidates = []
		downstream_candidates = []

		with open(args.annotation) as infile:
			for line in infile:
				cur = line.rstrip().split('\t')
				if cur[3] == args.upstream_label:
					upstream_candidates.append(int(cur[2]))
				if cur[3] == args.downstream_label:
					downstream_candidates.append(int(cur[1]))

		upstream_candidates.sort()
		downstream_candidates.sort()

		upstream_itr_end = upstream_candidates[0]
		downstream_itr_start = downstream_candidates[-1]

	elif args.downstream_itr_start != None and args.upstream_itr_end != None:
		upstream_itr_end = args.upstream_itr_end
		downstream_itr_start = args.downstream_itr_start
	else:
		sys.exit("Error: Please provide either (--annotation and --upstream-label and --downstream-label) or (--upstream-ITR-end and --downstream-ITR-start)")

	return upstream_itr_end, downstream_itr_start

def cigar_length(cigar):
	bases = 0
	for i in cigar:
		if i.type == '=' or i.type == 'I' or i.type == 'M' or i.type == 'X':
			bases += i.size
	return(bases)

def extract_SA_field(input_list):
	for field in input_list:
		if field[0] == 'SA':
			return(field[1])

def prepare_SA_field(SA):
	data = []
	for section in SA.split(';')[:-1]: # Since SA string ends in ; there will be an empty element at end of list. This is removed by [:-1]
		reference, position, strand, cigar, mapQ, nm = section.rsplit(',', maxsplit=5)
		data.append((reference, cigar_length(HTSeq.parse_cigar(cigar, int(position), reference, strand))))
	return(data)

def process_alignment(target_reference, reference, cigar, tags):
	noncontinuous = 0
	impurity = 0
	if cigar == None:
		sections = []
	else:
		sections = [(reference, cigar_length(cigar))]
	sections.extend(prepare_SA_field(extract_SA_field(tags)))
	for section in sections:
		if section[0] == target_reference:
			noncontinuous += section[1]
		else:
			impurity += section[1]

	return noncontinuous, impurity

def parseArguments():
	parser = argparse.ArgumentParser(prog="process_alignments", description='', usage='%(prog)s [options]')
	
	args = parser.add_argument_group('Input')
	args.add_argument('-a', '--annotation', help='Position of downstream ITR start. Provide this annotation file and labels for the ITRs OR both the upstream ITR start an downstream ITR end.', metavar='', dest='annotation')
	args.add_argument('-l', '--upstream-label', help='Upstream ITR label in the annotation file.', metavar='', dest='upstream_label')
	args.add_argument('-m', '--downstream-label', help='Downstream ITR label in the annotation file.', metavar='', dest='downstream_label')
	args.add_argument('-e', '--upstream-ITR-end', type=int, help='Position of upstream ITR end', metavar='', dest='upstream_itr_end')
	args.add_argument('-s', '--downstream-ITR-start', type=int, help='Position of downstream ITR start', metavar='', dest='downstream_itr_start')
	args.add_argument('-t', '--target', required=True, help='Name of target', metavar='', dest='target')
	args.add_argument('-p', '--per-read', required=True, help='Per read output from LAAVA', metavar='', dest='per_read')
	args.add_argument('-b', '--tagged-bam', required=True, help='Tagged bam output from LAAVA', metavar='', dest='tagged_bam')
	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)