#!/bin/bash
#
# $0 <SURPI-output> <genus> <output-root>

surpiout=$1
genus=$2
root=$3

mapdir=$(mktemp -d -p .)
tempdir=$(mktemp -d -p .)
topgi=$(ls $surpiout/OUTPUT_*/genus.bar.*/*$genus*|grep -oE "SnRa.*"|cut -f6 -d'.'|sort |uniq -c|sort -nr|head -1|awk '{print $2}')
for i in $surpiout/DATASETS_*/genus.bar.*; do cp $i/*$topgi*map $mapdir ; done 
z=1 ; 
for i in $mapdir/*.map; do 
	echo "coor" > $tempdir/temp0 ; cut -f1 $i >>  $tempdir/temp0; 
	echo ${i##*/}|cut -f4 -d '.' >  $tempdir/temp$z; 
	cut -f2 $i >>  $tempdir/temp$z ; 
	((z++)); 
done 
paste  $tempdir/temp* > ${root}_${genus}.${topgi}_comb.tsv; 
rm   -rf $tempdir

# bash consolidate_maps.sh < Respiro_cons.map > Respiro_cons_count.map

z=1
while read line; do
	if [ $z -eq 1 ]; then 
		((z++)); 
		continue; 	
	fi
	coor=$(cut -f1 -d' ' <<<$line)
	data=$(cut -f2- -d' '<<< $line)
	sum=0;
	#oper=${data// /+}
	declare -a mydata=($data)
	for i in ${data[@]}; do
		if [ $i -ge 1 ]; then
			((sum++))
		fi
	done
	echo -e "$coor\t$sum"
#	echo -e "$coor\t$((oper))"

done < ${root}_${genus}.${topgi}_comb.tsv > ${root}_${genus}_${topgi}_count.tsv

peak_map2.sh   ${root}_${genus}_${topgi}_count.tsv >  ${root}_${genus}_${topgi}_peak.tsv
