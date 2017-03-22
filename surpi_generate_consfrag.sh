#!/bin/bash
#
# $0 index surpioout genus output-root
# index file has to be named the same as .fa

# TO do:

# get fragment 
# awk '{print substr($0, start, go)}'
# make index

usage () {
	echo "$0 [ (-i <gi>) / (-g) choose-by-genus / (<index-base>) ] <SURPI-OUTPUT_DIR> <genus> <outputdir> [mincoverage]"; exit
}

if [ $1 == "-h" ]; then usage; fi

module load samtools bwa bcftools

### Pull reference and make a bwa index
if [ $1 == "-i" ] ; then 
	echo "Initialize and get reference and make index $2.$4.fa (genus.gi.fa) (start, len)..." >&2
	get_genbankfasta.pl $2 > $4.$2.fasta
### Fix this code sometime to pull best common fragment out
#	if [[ $4 ]] && [[ $5 ]]; then
#		fastascrub.sh -i $2.$4.fasta
#		endc=$(($2 + $3))
#		awk -v st="$4" -v len="$5" '$1 !~ "^>" {print substr($0, st, len)}; $1 ~ "^>" {print$0}' >$2.$3.${4}-${endc}.fa
#		bwa index $2.$3.${4}-${endc}.fa
#	else
		bwa index $4.$2.fasta
		index=$4.$2.fasta
		shift;  #only since -i takes an argument
#	fi
#	exit 0
fi



surpiout=$2
genus=$3
output=$4
mincov=${5:-1}
mkdir -p $output

### Use genus to pick a gi and index
if [ $1 == "-g" ] ; then
	mygi=$(ls $surpiout/genus.bar.*/*$genus*|grep -oE "SnRa.*"|cut -f6 -d'.'|sort |uniq -c|sort -nr|head -1|awk '{print $2}')
	echo "Best gi for $genus is $mygi" >&2
	get_genbankfasta.pl $mygi > $genus.$mygi.fasta
	bwa index $genus.$mygi.fasta
	index=$genus.$mygi.fasta
fi

index=${index:-$1}





for i in $surpiout/genus.bar*.fasta; do 
	if ls $i/*$genus* &>/dev/null ; then
		fastafile=$(ls $i/*$genus*.fa)
##		if [[ ! $fastafile =~ "48752"|"48681"|"S13" ]]; then continue; fi
		cp $fastafile $output
		fastafile=$output/${fastafile##*/}
		samfile=${fastafile%%.Blastn*}.sam; samfile=$output/${samfile##*/}
##		bowtie2 --very-sensitive-local -f -p 8 -x $index -U $fastafile -S $samfile
		bwa mem -k 8 $index $fastafile > $samfile
		bamfile=${samfile%.sam}.bam
		samtools view -S -b $samfile > $bamfile
		rm $samfile
		if [ $(samtools view -F4 $bamfile|head|wc -l) -gt 0 ]; then
			samtools sort $bamfile -o $bamfile.s
			mv $bamfile.s $bamfile
			samtools mpileup -d 100000 -uf $index $bamfile |bcftools call -c|vcfutils.pl vcf2fq -d $mincov > ${bamfile%.bam}.cons.fq
			sed -i "s/^@.*$/@$Consensus_${fastafile##*/}_${genus}/g" ${bamfile%.bam}.cons.fq
			fastqconvert.pl ${bamfile%.bam}.cons.fq |sed  's/^n*//g' > ${bamfile%.bam}.cons.fa # sed is to remove starting n's
			tail -n +2 ${bamfile%.bam}.cons.fa | tr '[:lower:]' 'n' |sed  's/^n*//g' >${bamfile%.bam}.cons.fa.temp2
			head -1 ${bamfile%.bam}.cons.fa > ${bamfile%.bam}.cons.fa.temp1
			cat ${bamfile%.bam}.cons.fa.temp1 ${bamfile%.bam}.cons.fa.temp2 > ${bamfile%.bam}.cons.fa 


			rm ${bamfile%.bam}.cons.fa.temp1 ${bamfile%.bam}.cons.fa.temp2
			rm ${bamfile%.bam}.cons.fq
		else
			rm $bamfile
		fi
	fi
done

rm $output/*.ex.fa.fai  $output/*.ex.fa
