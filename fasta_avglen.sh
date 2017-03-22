#!/bin/bash
#
# Gets average length of sequences in fasta file


fastafile=$1
counter=0
entrycount=0

while read line; do
	if [[ $line =~ ">" ]]; then 
		
		((entrycount++))
	else
		counter=$((counter + $(echo $line | tr -d  '[[:space:]]' | wc -c)  ))  
	fi
done < $fastafile

echo "$counter / $entrycount" |bc -l
echo -en "\n"

