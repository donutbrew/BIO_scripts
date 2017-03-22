#!/bin/bash

if [ $1 == "-h" ]; then echo "$0 <countfile> [<outdir>]" ; fi
countfile=$1
outdir=${2:-surpibycol}

cutoff=1
wholetax=1
combine=0

mkdir -p $outdir

if    [[ $countfile =~ family ]]; then startcol=2
 elif [[ $countfile =~ genus ]]; then startcol=3
 elif [[ $countfile =~ species ]]; then startcol=4
 elif [[ $countfile =~ "gi.counttable" ]]; then startcol=5
fi

cleancountfile=$(mktemp)
cat $countfile | sed "s:\t\t:\t-\t:g" |sed "s:\t \t:\t-\t:g"| tr ' ' '_'|sed 's/_\t/\t/g' > $cleancountfile


numsamples=$(head -1 $countfile |tr '\t' '\n'|grep -c "bar#")
numcols=$((numsamples + $startcol - 1))


declare -a samplenames
header=1
while IFS= read -r line; do
	if [ $header -eq 1 ]; then outfile="$outdir/${line##*#}"
		samplenames=(${line//bar#/})
		header=0
	else	echo $(echo $line|cut -f1 -d$' ')
		declare -a reads=(${line/  / - })
		for ((coln=$((startcol-1)); coln<=$((numcols-1)); coln++)); do
			outfile=$outdir/${samplenames[$coln]}.tmp
			if [ ${reads[$coln]} -ge $cutoff ]; then
				if [ $wholetax -eq 1 ]; then
					for ((j=$((startcol - 2)); j >= 0; j--)); do
						echo -en "${reads[$j]}"
						if [ $j -ne 0 ]; then echo -en ":"; else echo -en "\t"; fi
					done 
				else
					echo -en "${reads[0]}\t"
				fi
				echo ${reads[$coln]}
			fi >> $outfile	
			
		done 
	fi
	sed -i '/^-\t-$/d' $outfile
done < $cleancountfile

for i in $outdir/*tmp; do
	name=${i%.tmp}; name=${name##*/}
	( echo ${name} && sort -nrk2 $i && rm $i )  > ${i%.tmp}.txt &
done

wait
#NEEDS TOP5 FN

if [ $combine -eq 1 ]; then
	linecount=$(mktemp)
	wc -l $outdir/*| sort -nr|tail -n +2 > $linecount
	maxline=$(head -1 $linecount|awk '{print$1}')
	while read line; do
		nline=$(echo $line|awk '{print$1}')
		fileloc=$(echo $line|awk '{print$2}')
		addlines=$((maxline - nline))
		for ((i=0;i<$addlines;i++)); do
			echo -en "\n-\t-" >> $fileloc
		done
	done < $linecount
#	for i in $outdir/*.txt; do sed -i 's/\t/,/g' $i;  done
	paste $outdir/* > $outdir/COMBINED.txt
fi
rm $cleancountfile
