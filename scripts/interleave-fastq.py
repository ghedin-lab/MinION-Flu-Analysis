#!/usr/bin/python

import sys
import gzip


def interleave(r1,r2):
	while True:
		line = r1.readline()
		if line.strip() == "":
			break
			
		print line.strip()
		
		for i in xrange(3):
			print r1.readline().strip()
		
		for i in xrange(4):
			print r2.readline().strip()

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print "Please provide two fastq files"
		sys.exit(1)
	
	f1 = sys.argv[1]
	f2 = sys.argv[2]
	
	if "gz" in f1:
		r1 = gzip.open(f1, "rb")
		r2 = gzip.open(f2, "rb")
		interleave(r1, r2)
		r1.close()
		r2.close()
	else:
		with open(f1, "rb") as r1:
			with open(f2, "rb") as r2:
				interleave(r1, r2)