#!/bin/bash

tempdir=$(mktemp -d)
threads=1
while getopts "gt:" opt; do
	case $opt in
                    #    f)      fasta=$OPTARG ;;
			g)	get_gi=1 ;;
			t)	threads=$OPTARG ;;
                        *)      echo "Usage: $0 [-g get by gi] <gi or fasta file> <fastq file>" ; exit 1        ;;
                 esac
done
shift $((OPTIND-1))

if [ "$#" -lt 2 ] ; then exit 1; fi

if  [ "$get_gi" -eq 1 ]; then 
	get_genbankfasta.pl $1 > $tempdir/temp.fa
else 	
	cp $1 $tempdir/temp.fa
fi

bowtie2-build -f $tempdir/temp.fa $tempdir/tempidx

bowtie2 --fast-local -p $threads -x $tempdir/tempidx -U $2 -S ${3:-BOWTIE2OUT.sam}
 rm -rf $tempdir
#echo $tempdir






