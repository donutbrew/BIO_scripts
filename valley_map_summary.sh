#!/bin/bash
stuff=${1:-'*gap*.tsv'}
echo -e "Name\tTotal_gap\tGaps..."
for i in $stuff; do 
	gaps=""; totalgap=0
	while read -r line; do
		a=$(echo $line|cut -f1 -d' ')
		b=$(echo $line|cut -f2 -d' ')
		c=$(echo $line|cut -f3 -d' ')
		if [[ $a -eq $b ]]; then continue; fi
		gaps="$gaps\t$a..$b"
		totalgap=$((totalgap + c))
	done < <(tail -n +2 $i)
	echo -e "${i%%_L001*}\t$totalgap$gaps"
done
