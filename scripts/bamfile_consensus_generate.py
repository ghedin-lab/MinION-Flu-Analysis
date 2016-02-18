import pysam
import sys
if (len(sys.argv) > 1):
	namedfile = sys.argv[1]
	segment = sys.argv[2]
# THIS IS WAHT IT SHOULD LOOK LIKE U NOOB
#python bamfile_consensus_generate.py file.bam H1N1-HA
samfile = pysam.AlignmentFile(namedfile, "rb" )

print ">"+segment
conseq = []
for pileupcolumn in samfile.pileup(segment):

	#count = 0	

	ntdict = {}
	for pileupread in pileupcolumn.pileups:

	#	count += 1

		# if pileupread.is_refskip == 1:
		# 	print 'WAHTAHTHA'
		if not pileupread.is_del:# and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
			soment = pileupread.alignment.query_sequence[pileupread.query_position]
			if ntdict.has_key(soment):
				ntdict[soment] = ntdict[soment] + 1
			else:
				ntdict[soment] = 1
		else:
			soment = '-'
			if ntdict.has_key(soment):
				ntdict[soment] = ntdict[soment] + 1
			else:
				ntdict[soment] = 1
			# print pileupread.is_refskip
	#print ntdict
	#print count
	sorted_dict = sorted(ntdict.iteritems(), key=lambda item: -item[1])
#	print sorted_dict
	conseq.append(sorted_dict[0][0])
print ''.join(conseq)
samfile.close()
