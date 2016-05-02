#!/bin/bash 

seg='PB2'
segNum='1'
segSize='2.3k'

cd "$(dirname "$0")"

grep -P -B6 'Influenza B.+Segment:'$segNum flu-11-9.2d.fludb.xml | grep 'Iteration_query-def' | \
perl -pe 's/^.+>(.+) FAS.+/$1/g' | grep --no-group-separator -A3 -F -f - flu-11-9.2d.fastq  > flu-11-9.2d.flub-$seg-reads.fastq

for i in .5 .3 .1 .035 .01
do

	echo $i

	canu -p flu -d flu-11-9-flub-$seg-only-low-cov-canu-out genomeSize=$segSize corMhapSensitivity=high corMinCoverage=2 \
	errorRate=0.035 minOverlapLength=499 corMaxEvidenceErate=$i -nanopore-raw flu-11-9.2d.flub-$seg-reads.fastq
	
	grep -c '^>' flu-11-9-flub-$seg-only-low-cov-canu-out/flu.unassembled.fasta

done

:<<'END'

cd flu-11-9-flub-$seg-only-low-cov-canu-out/

# align to each other to see how different they are
nucmer --maxmatch flu.unassembled.fasta flu.unassembled.fasta
show-coords -lcT <(delta-filter -q -l 500 out.delta) | perl -pe 's/\t/ | /g' > flu.coords

ln -s ../flu-11-9.2d.fastq 

bwa index flu.unassembled.fasta
bwa mem -t 30 -x ont2d flu.unassembled.fasta flu-11-9.2d.fastq > flu-11-9.sam

samtools view -b -o flu-11-9.bam flu-11-9.sam
samtools sort -o flu-11-9.$seg.sort.bam flu-11-9.bam
samtools index flu-11-9.$seg.sort.bam

rm flu-11-9.sam flu-11-9.bam

samtools idxstats flu-11-9.$seg.sort.bam 