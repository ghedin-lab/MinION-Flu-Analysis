#!/usr/bin/python2.6

import sys
import re
from Bio import AlignIO

if (len(sys.argv) != 3):
  sys.exit("specify: unaligned_sequence_file unaligned_sequence_file_format (output will be fasta-format to stdout)")

alignment = AlignIO.read(open(sys.argv[1]), sys.argv[2])

count = 0
mismatch = 0
deletion = 0

print "RefName\tRefLen\t#Variants\t#Deletions"

for record in alignment :

	if count == 0:
		ref = record.seq
		refName = record.id
		count = 1
	else:
		for i in range(0,len(record.seq)):
			if record.seq[i] == "-":
				deletion += 1
				print refName,i,"del"
			elif record.seq[i] != ref[i]:
				mismatch += 1
				print refName,i,"snv"
			else:
				print refName,i,"NA"

		print refName,"\t",len(ref),"\t",mismatch,"\t",deletion
