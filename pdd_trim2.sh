#!/bin/bash
#
# PDD Standardized quality trimming
# Default Trim B2 adapter, Illumina adapter, Quality trim to Q20, Dust, Min Length 50
# 
# Requires cutadapt and prinseq-lite
#
# Usage: pdd_trim.sh -1 <R1> [-2 <R2>] -o <output-file prefix>
#
# Output is one or two FASTQ files with the name XXX_Clean.fastq / XXX_1_Clean.fastq XXX_2_Clean.fastq
#   If no -o is specified, XXX is everything before the first underscore
#
#
# Author: Clint Paden     Modified: 08 Apr 2015
# 

min_q=25
min_len=50
threads=8

sge=0;read1=""; read2=""; root="";zipout=0;filet='*.gz';batchmode=0;prinseqrun=1

error () {
	echo $1
	exit 1
}
while getopts "1:2:o:zpd:st:" opt; do
	case $opt in
		1)	read1=$OPTARG ;;
		2)	read2=$OPTARG ;;
		o)	root=$OPTARG ;;
		z)	zipout=1 ;;
		d)	batchdir=$OPTARG ; batchmode=1 ;;
		s)	sge=1 ;;
		t)	threads=$OPTARG ;;
		p)	prinseqrun=0 ;;
		*)	echo "Usage: pdd_trim.sh -1 <R1.fastq/.gz> [-2 <R2.fastq/.gz>] -o <output-file prefix> [-p (skip dusting)]"
			exit 1
			;;
	 esac
done
shift $((OPTIND-1)) 

if [[ ! -f "$read1" && $batchmode -ne 1 ]] ; then 
	echo "Error read file 1 not found!" ; exit 1
elif [ ! -f "$read2" ] ; then
	echo "Single end mode"
	semode=1
else
	echo "Paired end mode"
	semode=0
fi

if [ $batchmode -eq 1 ]; then 
	ls $filet >/dev/null || filet='*.fastq'; ls $filet >/dev/null || filet='*.fq'; ls $filet >/dev/null || error "CANNOT DETECT FILETYPE" ;
	for i in $batchdir/$filet; do 
		if [[ $i =~ "_R1" && $semode -eq 0 ]]; then read1=$i
		elif [[ $i =~ "R2" && $semode -eq 0 ]]; then read2=$i
			if [ $sge -eq 1 ]; then
				each./pdd_trim2.sh -d . -o ./trim/root=${read1##*/}; eachroot=${eachroot%.fastq}; eachroot=${eachroot%%_L00*};
				mkdir -p $root
				echo "Command: qsub -V -cwd -pe smp $threads -b y -N ptr_$eachroot pdd_trim2.sh -1 $read1 -2 $read2 -t $threads -o $root "
				qsub -V -cwd -pe smp $threads -b y -N ptr_${read1##*/} pdd_trim2.sh -1 $read1 -2 $read2 -t $threads -o $root
			else 
				$0 -1 $read1 -2 $read2 -o $root &
			fi
		elif [[ $semode -eq 1 ]]; then
			read1=$i
			if [ $sge -eq 1 ]; then
				eachroot=${read1##*/}; eachroot=${eachroot%.fastq}; eachroot=${eachroot%%_L00*};
				mkdir -p $root
				echo "Command: qsub -V -cwd -pe smp $threads -b y -N ptr_$eachroot pdd_trim2.sh -1 $read1 -2 $read2 -t $threads -o $root "
				qsub -V -cwd -pe smp $threads -b y -N ptr_${read1##*/} pdd_trim2.sh -1 $read1 -t $threads -o $root
			else 
				$0 -1 $read1 -o $root &
			fi
		else
			echo "ERROR WITH FILES" >&2
			exit 1
		fi
	done
	echo "Jobs submitted"
	if [ $sge -ne 1 ];  then echo "waiting on jobs to finish..." ; wait; fi
	exit 0
fi




if [ -z "$root" ] ; then
	echo "-o not set, using default"
	root=${read1##*/}
	root=${root%%_*}
	if [ -e "$root.1.trim.fastq" ] ; then 
		echo "Output exist, use -o option"; exit 1
	fi 
fi

# The case where -o is a directory (e.g. for batch processing)
if [ -d "$root" ] ; then
	echo "$root is a directory, will deposit in there using the R1 filename as a prefix"
	tmpstem=${read1##*/}
	tmpstem=${tmpstem%%_R*}
	root="$root/$tmpstem"
fi
# State outputs
if [ $semode -eq 1 ]; then
#	echo "Sorry, single end not supported" ; exit 1;
	echo "Final output file will be  $root""_Clean.fastq" >&2
else
	echo "Final output files will be  $root""_1_Clean.fastq and $root""_2_Clean.fastq" >&2
fi
# Use cutadapt to remove adapters 
#    Currently, this is for Primer B2 and Illumina adapters. May make user configurable later.
if [ $semode -eq 1 ]; then
echo "SE NOT SUPPORTED for parallel, using regular";
	cutadapt -g GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGC -q $min_q -m $min_len -n 5 -o $root.1.trim.fastq $read1  
else		    
	echo "Command: cutadapt-parallel.sh --threads $threads -g GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -G GTTTCCCAGTCACGATA -A TATCGTGACTGGGAAAC -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q $min_q -m $min_len -n 5 -o $root.1.trim.fastq -p $root.2.trim.fastq $read1 $read2"
	cutadapt-parallel.sh --threads $threads -g GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -G GTTTCCCAGTCACGATA -A TATCGTGACTGGGAAAC -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q $min_q -m $min_len -n 5 -o $root.1.trim.fastq -p $root.2.trim.fastq $read1 $read2
fi

# Use Prinseq-lite to dust junk and Quality trim - currently no bad stuff gets out
if [ $prinseqrun -eq 1 ] ; then 
	if [ $semode -eq 1 ]; then
		prinseq-lite.pl -fastq $root.1.trim.fastq -out_good $root"_Clean" -lc_method dust -lc_threshold 7 -min_len $min_len # -out_bad $root.trim.dust_bad 
	else
		prinseq-lite.pl -fastq $root.1.trim.fastq -fastq2 $root.2.trim.fastq -out_good $root"_Clean" -lc_method dust -lc_threshold 7 -min_len $min_len 
	fi
else
	mv  $root.1.trim.fastq ${root}_Clean_1.fastq ; mv $root.2.trim.fastq ${root}_Clean_2.fastq
fi

# Option to compress output
if [ "$zipout" -eq 1 ]; then 
	for i in "$root*_Clean.fastq"; do 
		gzip $i &
		waitpid=$$
	done
	while ps axc|grep $waitpid >/dev/null; do sleep 8s; done # Wait on the last gzip proc. Last proc may not be the best indicator, but whatev.
fi
		
# Cleanup

#rm "${root}"*.trim.fastq
#rm "${root}"*singletons.fastq
#rm "${root}"*prinseq_bad*.fastq
