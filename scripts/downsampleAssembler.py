#!/usr/bin/python

import argparse
import sys
import subprocess
#from subprocess import Popen, PIPE
import pysam
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument("-r","--reference",help="Reference file",type=str,required=True)
group = parser.add_mutually_exclusive_group()
group.add_argument("-f","--fastq",help="Prefix of paired fastq filenames used to map",type=str, default=None)
group.add_argument("-b","--bam",help="Bam file used to partition reads",type=str, default=None)
# parser.add_argument("-f","--fastq",help="Prefix of paired fastq filenames",type=str,required=False)
# parser.add_argument("-b","--bam",help="Bam file used to partition reads",type=str,required=False)
parser.add_argument("-c","--coverage",help="Target average coverage to downsample to", type=int,default=100)
args = parser.parse_args()

if args.fastq is None and args.bam is None:
	parser.error("argument -b/--bam or -f/--fastq must be specified")

spadesPath = "/usr/local/Cellar/spades/3.7.0/bin/spades.py"


def getRefNames(bam):
	f = pysam.AlignmentFile(bam,"rb")
	refdict = {}
	for i in f.header['SQ']:
		# remove when finished
		if 'FluB' in i['SN']:
			refdict[i['SN']] = i['LN'] 
	return refdict


def calcAvgCov(bam, ref, refLen):
	f = pysam.AlignmentFile(bam,"rb")
	counts = f.count_coverage(ref,0,refLen)
	avgCov = 0
	for i in range(0,len(counts)):
		for x in range(0,len(counts[i])):
			avgCov += counts[i][x]
	return float(avgCov/refLen)
	
	
def indexBam(bam):
	subprocess.call(["samtools","index",bam])
	
	
def downsampleBam(bam, ref, refLen, avgCov):
	outBam = ref+"."+str(args.coverage)+"x.bam"
	s = args.coverage/avgCov
	
	ps = subprocess.Popen(("samtools","view","-hs",str(s),bam),stdout=subprocess.PIPE)
	output = subprocess.check_output(("samtools","view","-sB","-","-o"+outBam),stdin=ps.stdout)
	ps.wait()
	indexBam(outBam)
	
	#print ref,calcAvgCov(outBam,ref,refLen)
	return outBam


def partitionBam(bam, ref):
	outBam = ref+".bam"
	ps = subprocess.Popen(("samtools","view","-h",bam, ref), stdout=subprocess.PIPE)
	output = subprocess.check_output(("samtools", "view", "-Sb", "-", "-o"+outBam),stdin=ps.stdout)
	ps.wait()
	indexBam(outBam)
	return outBam

def bamToFastq(bam):
	m = re.match("(^.+).bam",bam)
	outFastq = m.group(1)
	outFastq = outFastq+".fastq"
	with open(outFastq,"w") as f:
		subprocess.call(["bamtools","convert","-in",bam,"-format","fastq"],stdout=f)
	return outFastq

def assemble(fastq):
	m = re.match("(^.+).fastq",fastq)
	spadesOut = m.group(1)
	spadesOut = spadesOut+"-spades-out"
	subprocess.call(["spades.py","-s",fastq,"-o",spadesOut],stdout=open(os.devnull, 'wb'))
	return spadesOut

def main():
	refNames = getRefNames(args.bam)
	
	avgCov = {}
	
	for ref in refNames:
		partBam = partitionBam(args.bam,ref)
		avgCov[ref] = calcAvgCov(partBam,ref,refNames[ref])
		downBam = downsampleBam(partBam,ref,refNames[ref],avgCov[ref])
		downFastq = bamToFastq(downBam)
		spadesOut = assemble(downFastq)


if __name__ == "__main__":
	main()