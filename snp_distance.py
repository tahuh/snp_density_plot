#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, required=True, help="Infile SNP")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file")

args = parser.parse_args()
outfile = open(args.outfile, "w")
outfile.write("distance\tcount\n")
d = {}
with open(args.infile) as infile:
	# Do not care memory usage
	lines = infile.readlines()
	while 1:
		line = lines.pop(0)
		if line[0] == '#':
			continue
		else:
			break # chrom line
	N = len(lines) - 1
	for i in range(N):
		line1 = lines[i].rstrip()
		line2 = lines[i+1].rstrip()
		data1 = line1.split("\t")
		data2 = line2.split("\t")
		is_exon1 = (data1[5] == "YES")
		is_exon2 = (data2[5] == "YES")
		if data1[2] == '.' or data2[2] == '.':
			continue
		if is_exon1 or is_exon2:
			continue
		dist = abs(int(data1[1]) - int(data2[1]))
		try:
			d[dist] += 1
		except KeyError:
			d[dist] = 1

for k, v in d.items():
	outfile.write(str(k) + '\t' + str(v) + '\n')
