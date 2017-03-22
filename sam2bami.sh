#!/bin/bash

if [ $1 == "-t" ] ; then 
	shift
	threads=$1;
	shift
fi
 
threads=${threads:-1}
samfile=$1
bamfile=${samfile%.sam}.bam
if [[ $samfile =~ ".bam" ]] ; then 
	echo "BAM input" >&2
	bamfile=$samfile
	samtools sort -@ $threads ${bamfile} > ${bamfile}.temp2
else
	samtools view -S -b $samfile > ${bamfile}.temp
	samtools sort -@ $threads ${bamfile}.temp > ${bamfile}.temp2
fi
mv ${bamfile}.temp2 $bamfile
rm ${bamfile}.temp
samtools index $bamfile
