#!/bin/bash

if [[ $1 == "-h" ]]; then echo -e "Usage: $0 [-f Fix Headers] readcounts.<XXX>.BarcodeR1R2.log\n"; exit 0; fi
if [[ $1 == "-f" ]]; then fix=1; shift; fi

defrc=$(ls readcounts.*.BarcodeR1R2.log)
readcountfile=${1:-$defrc}
typelist=$(mktemp)
samples=$(mktemp)
#declare -a mytempf
mytempf=$(mktemp)
idx=0


fixnames () {
	name=$(echo $1|cut -f2- -d.)
	echo -n $name|sed 's/preprocessed.fastq/processed/;s/preprocessed.s20.h250n25d12xfu.human.snap.unmatched.fastq/Human/;s/NT.snap.matched.fulllength.all.annotated.sorted/NT-annotated/;s/NT.snap.matched.fl.Viruses.annotated/NT-Virus/;s/NT.snap.matched.fl.Bacteria.annotated/NT-Bacteria/;s/NT.snap.matched.fl.nonChordatEuk.annotated/NT-otherEuk/;s/NT.snap.unmatched.sam/NT-unmatched/;s/Contigs.and.NTunmatched.Viral.RAPSearch.e1.NR.e1.Viruses.annotated/NR-Virus/g;' $line
}


cut -f1 $readcountfile|uniq > $typelist
typelista=($(cut -f1 $readcountfile|uniq|tr '\n' '\t'))
(echo "Sample"; cut -f2 -d'#' $readcountfile |cut -f1 -d'/'|sort|uniq) > $samples
while read -r samplename; do
	if [ $samplename == "Sample" ]; then continue; fi
#	while read -r type; do
	for type in ${typelista[@]}; do
		if [[ $idx -eq 0 ]]; then
			for i in ${typelista[@]}; do fixnames $i; echo -en "\t";done >> $mytempf
			idx=1
			echo -en "\n" >> $mytempf
		fi
#echo $type >&2
 		readn=$(grep -P "^$type\b" $readcountfile|grep "$samplename" |awk '{print$3}') # >> ${mytempf[$idx]}
		if [[ "$type" =~ "preprocessed.fastq" ]]; then prepro=$readn; fi
		if [[ "$type" =~ "human.snap" ]]; then  readn=$((prepro - readn)) ;fi
		if [[ "$readn" == " " ||  "$readn" == "" ]]; then echo -en "00\t" >> ${mytempf}
		else echo -en "$readn\t" >> $mytempf
		fi

	done #< $typelist
	echo -en "\n" >> $mytempf
	idx=1
done < $samples
paste $samples ${mytempf} 
rm $samples $mytempf  $typelist

