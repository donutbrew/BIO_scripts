#!/bin/bash

if [ "$1" == "-h" ] ; then echo "Usage: $0 <2-col \"map\" file> [coverage level] > tsv file with [c1 c2 lin-size max-coor max-val";  exit ; fi

newcoor=0
firstc=0
sdata=0
#firstline=1
declare -a mapfile
declare -a segcount
mincoverage=${2:-1}

while IFS=$'\t' read -r -a mapline; do mapfile["${mapline[0]}"]=${mapline[1]} ; done < $1

echo -e "START_COOR\tEND_COOR\tLINEAR_SIZE\tMAX_COOR\tMAX_COV\tAVG_COV"

for i in ${!mapfile[@]}; do
	oldcoor=$newcoor
	newcoor=$i
	coverage=${mapfile[$i]}
	
	printf -v coverage "%.f" "$coverage" # rounds scientific notation to nearest integer...not the most precise.
	
	if   [[ $coverage -lt $mincoverage && ${mapfile[$oldcoor]} -lt $mincoverage ]]; then continue;
	elif [[ $coverage -lt $mincoverage && ${mapfile[$oldcoor]} -ge $mincoverage ]]; then
	######## Print data fro segment and reset segment
		max=0; maxi=0; init=0;counter=1; average=0
		for n in ${!segcount[@]}; do
			if [[ $init -eq 0 ]]; then firstc=$n; init=1; fi 
			if (( ${segcount[$n]} > $max )); then max=${segcount[$n]}; maxi=$n; fi
			average=$(($average + ${segcount[$n]})) # just adding here
			((counter++))
		done
		average=$(echo "scale=2;$average / $counter" |bc)
#		if [[ $firstline -eq 1 ]]; then echo -e "START_COOR\tEND_COOR\tLINEAR_SIZE\tMAX_COOR\tMAX_COV"; firstline=0;fi
		echo -e "${firstc}\t${n}\t$((n-firstc))\t$maxi\t$max\t$average" 
		unset segcount
		declare -a segcount
		sdata=0 #indicated that stored data has already been printed
	else
	######## Add to Segment
		segcount[$i]=$coverage
		sdata=1 # indicates there is data stored in segcount	
	fi
done


#to deal with remaining unprinted data
if [ $sdata -eq 1 ]; then
                max=0; maxi=0; init=0;counter=1; average=0
                for n in ${!segcount[@]}; do
                        if [[ $init -eq 0 ]]; then firstc=$n; init=1; fi 
                        if (( ${segcount[$n]} > $max )); then max=${segcount[$n]}; maxi=$n; fi
                        average=$(($average + ${segcount[$n]})) # just adding here
			((counter++))
                done
                average=$(echo "scale=2;$average / $counter" |bc)
	
#               if [[ $firstline -eq 1 ]]; then echo -e "START_COOR\tEND_COOR\tLINEAR_SIZE\tMAX_COOR\tMAX_COV"; firstline=0;fi
                echo -e "${firstc}\t${n}\t$((n-firstc))\t$maxi\t$max\t$average" 




		
fi
