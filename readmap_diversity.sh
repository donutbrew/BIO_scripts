#!/usr/bin/env bash
#
# Script to find read variation across a genome
#
#
# Author: Clint Paden <cpaden@cdc.gov>
# Date: 13 July 2015
#
# Input: Paired read data, reference
# Output: Consensus sequences, BAM for reference, BAM for self reference (Future: VCF?)
#
# Usage: $0 -f <fasta-reference> -1 <read1> -2 <read2> -o <output prefix>

clean=1
trimends=1
keepclean=0
threads=6
avg=1

#tempdir=$(mktemp -d -p /scratch)
tempdir=$(mktemp -d)
#mkdir mytmp
#tempdir=mytmp
module load samtools bowtie2

while getopts "1:2:o:f:cmkt:a" opt; do
	case $opt in
		1)	read1=$OPTARG ;;
		2)	read2=$OPTARG ;;
		o)	root=$OPTARG ;;
		f)	reference=$OPTARG ;;
		c)	clean=0 ;;
		m)	trimends=0 ;;
		k)	keepclean=1 ;;
		t)	threads=$OPTARG ;;
		a)	avg=0 ;; #no average
		*)	echo "Usage: $0 -f <fasta-reference> -1 <read1> -2 <read2> -o <output prefix> [-c dont clean]  [-m dont trim from ends] [-t threads]"
			exit 1
			;;
	 esac
done
shift $((OPTIND-1)) 


[ -n "$read1" ] && [ -n "$read2" ] && [ -n "$reference" ] && [ -n "$root" ] ||  exit 1;

echo "Finding read diversity using: $read1 and $read2"
echo "Hostname is $(hostname)"
if [ $clean -eq 1 ]; then
	echo "Quality trimming . . . "
	pdd_trim.sh -1 $read1 -2 $read2 -o ${root} -p
	read1=${root}_Clean_1.fastq
	read2=${root}_Clean_2.fastq
 
else echo "Skipping adapter and quality trimming"
fi

if [ $trimends -eq 1 ] ; then 
	echo "Trimming ends . . ."
	seqtk trimfq -b 26 -e 26 $read1 >  $root.trimtmp1.fastq &
	seqtk trimfq -b 26 -e 26 $read2 > $root.trimtmp2.fastq	
	wait
	[ $clean -eq 1 ] && read1=${root}_Cleantrim_1.fastq || read1=${root}_trim_1.fastq 
	[ $clean -eq 1 ] && read2=${root}_Cleantrim_2.fastq || read2=${root}_trim_2.fastq
	mv $root.trimtmp1.fastq $read1
	mv $root.trimtmp2.fastq $read2
else echo "Skipping end trimming"
fi
 
echo "Mapping reads to reference . . ."
bowtie2-build $reference $tempdir/ref  &>> $root.bowtie2.log 
samfile1=$root.$(basename $reference); samfile1=${samfile1%.*}.sam
bamfile1=${samfile1%.*}.bam
bowtie2 -x $tempdir/ref -p $threads --fast-local -X 1000 -1 $read1 -2 $read2 -S $samfile1 &>> $root.bowtie2.log
samtools view -@ $threads -S -b $samfile1 > $bamfile1
samtools sort -@ $threads $bamfile1 -o ${bamfile1%.bam}.sorttmp.bam
mv ${bamfile1%.bam}.sorttmp.bam $bamfile1
samtools index $bamfile1
echo "Generating consensus sequence . . ."
samtools mpileup -uf $reference $bamfile1 | bcftools call -c - | vcfutils.pl vcf2fq -d 5 -D 1000000 > $root.cons.fastq

seqtk seq -A $root.cons.fastq > $root.cons.fasta

echo "Mapping reads to consensus . . ."
samfile2=$root.self.sam
bamfile2=$root.self.bam
bowtie2-build $root.cons.fasta $tempdir/self &>> $root.bowtie2.log 
bowtie2 -x $tempdir/self -p $threads --fast-local -X 1000 -1 $read1 -2 $read2 -S $samfile2  &>> $root.bowtie2.log
#samtools view -@ $threads -S -b $samfile2 > $bamfile2
#grep -v "^@"  $samfile2 | awk '$3 != "*" {print$4"\t"$17}' > $root.NM.tsv # collect mapping location and edit distance (col 17 for bt2 output)
echo "Getting edit distance information . . ."
if [ "$avg" -eq 0 ]; then 
	samtools view -S -F4 $samfile2 | awk '{print$4"\t"substr($17,6)}' | sort -k1 -n > $root.NM.tsv    # collect mapping location and edit distance (col 17 for bt2 output)
else
	samtools view -S -F4 $samfile2 | awk '{print$4"\t"substr($17,6)}' | sort -k1 -n > $root.NM.tsv
	averagebycol1.pl $root.NM.tsv |tail -n +2|sort -k1 -n > $root.NMavg.tsv    # collect mapping location and edit distance (col 17 for bt2 output)
fi

if [ $quality -eq 1 ]; then
	avgqscore.pl -S $samfile2 > $root.Q.tsv
fi

#if [ $clean -eq 1 ] || [ $trimends -eq 1 ] && [ $keepclean -ne 1 ] ; then rm $read1 $read2 $samfile1 $bamfile1; fi
echo  "Tool is complete. Please see $root.NM.tsv"

#rm -rf $tempdir

