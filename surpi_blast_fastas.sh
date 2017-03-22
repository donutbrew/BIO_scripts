#!/bin/bash
# 
# BASH script to go through a SURPI OUTPUT_* directory and BLAST all the reads under a certain * per genus
# in order to "verify" them. A followup parser will need to be generated.

surpioutput=$PWD
maxreads=20
outdir=""
blastdb="/blast/db/nt"
blastjobs=50
exclude="$HOME/surpi_tmp/genusfilterlist.txt"

module load ncbi-blast+ || exit 1

while getopts "d:m:o:h" opt; do 
	case $opt in 
		d) surpioutput=$OPTARG ;;
		m) maxreads=$OPTARG ;;
		o) outdir=$OPTARG ;;
		h) echo -e "Usage: $0 -d <SURPI OUTPUT Directory> -m <Only BLAST genuses with fewer than X reads> -o <blast output dir>"; exit 1 ;;
	esac
done

mkdir $outdir || exit 1


for gendir in $surpioutput/genus.bar.*fasta; do
	for fastafile in $gendir/SnRa.* ; do
		if [ $(grep -c ">" $fastafile) -lt $maxreads ]; then
			newhead=${fastafile##*/}; newhead=$(cut -f2,4,5 -d'.' <<< $newhead)
#			sed "s/>/>$newhead|/g"  $fastafile >> $outdir/query_file.fa
			sed "s/>/>$newhead|/g"  $fastafile >> $outdir/$newhead.query.fa
		fi
	done
done

mkdir $outdir/junk
find $outdir |grep -f $exclude|xargs -I {} mv {} $outdir/junk/
echo DONE STAGING; sleep 10s

mkdir -p $outdir/log
x=1
numthreads=1
for qfile in $outdir/*.query.fa; do 
	numseq=$(grep -c ">" $qfile)
	if [ $numseq -lt 5 ] ; then numthreads=1
	elif [ $numseq -lt 20 ] ; then numthreads=5
	else numthreads=10
	fi
	
#	qsub -V -cwd -b y -N blast_$x -pe smp $numthreads -o $outdir/log/ -e $outdir/log/ blastn -num_threads $numthreads -db $blastdb -max_target_seqs 10 -task blastn -outfmt \'6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle\' -query $qfile -out $qfile.blast
        qsub -V -cwd -b y -N blast_$x -pe smp $numthreads -o $outdir/log/ -e $outdir/log/ blastn -num_threads $numthreads -db $blastdb -max_target_seqs 10 -max_hsps 1 -task blastn -outfmt \'6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle\' -query $qfile -out $qfile.blast  # Added  -max_hsps 1 to make sure we only get one hit per subject 1-11-17

	((x++))
done
