#!/usr/bin/python

import sys
from Bio.Blast import NCBIXML as xml
import re

hits = dict()
# segments = dict()
#
# for i in ["H1N1","H3N2","FluB"]:
# 	segments[i] = dict()
# 	for x in ["PB1","PB2","PA","NA","HA","NP","NS","MP"]:
# 		segments[i][x] = 0

with open(sys.argv[1],"rb") as f:
	for record in xml.parse(f):
		# for result in record.hit_id:
		# 	print result
		for alignment in record.alignments:
			hitDef = alignment.hit_def
			seg = re.search('^.+Segment:(\d).+',hitDef)
			seg = int(seg.group(1))
			#print seg
			segment = ""
			if seg == 1:
				segment = "PB2"
			elif seg == 2:
				segment = "PB1"
			elif seg == 3:
				segment = "PA"
			elif seg == 4:
				segment = "HA"
			elif seg == 5:
				segment = "NP"
			elif seg == 6:
				segment = "NA"
			elif seg == 7:
				segment = "MP"
			elif seg == 8:
				segment = "NS"

			# use this for CDS only
#			seg = re.search('^.+Name:((.{2,3})(-|\s)|nonstructural).+',hitDef)
#			segment = seg.group(1)
# 			if segment in ["BM2","M1","M2","M42"]:
# 				segment = "MP"
# 			elif segment in ["NS1","NS2"] or segment == "nonstructural":
# 				segment = "NS"
# 			elif segment == "NB":
# 				segment = "NA"
			
			gb = re.search('^gb:(.{8}).+', hitDef)
			acc =  gb.group(1)
			sub = re.search('^.+Subtype:(.{4}).+', hitDef)
			subtype = sub.group(1)
			if subtype == "null":
				subtype = "FluB"
			# print "sub",subtype
		# 	print "seg",segment
		# 	print "gb",acc
		# 	print "def",hitDef
			if acc not in hits:
				hits[acc] = dict()
				hits[acc]["seg"] = segment
				hits[acc]["id"] = alignment.hit_id
				hits[acc]["def"] = hitDef
				hits[acc]["sub"] = subtype
				hits[acc]["count"] = 1
			else:
				hits[acc]["count"] += 1
				
for acc in hits:
	print acc+"\t"+hits[acc]["sub"]+"\t"+hits[acc]["seg"]+"\t"+str(hits[acc]["count"])+"\t"+hits[acc]["def"]
