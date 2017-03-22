#!/bin/bash
# Takes in rename file in format "Newname	Oldname" and file to rename

if [ $1 == "-h" ]; then echo "Usage: $0 [ -i (inline) ] <namekey(new\told)> <filein>; resutls to stdout unless -i is specified"; exit 1;fi

if [ $1 == "-i" ]; then inplace=1; shift; else inplace=0; fi

namekey=$1 
targetfile=$2
sedcom=""
append=1

if [ -f $targetfile ]; then 
while read line; do
		old=$(echo $line|awk  '{ print $2 }') ;# echo "Sample ID: $sampleno"
	        new=$(echo $line|awk  '{ print $1 }') ;# echo "Sample name is $sampname"
		if [ -z $old ] || [ -z $new ] ; then continue; fi
	        if [ $append -eq 1 ]; then
			sedcom=$sedcom"s:$old:${new}_$old:g;"
		else
			sedcom=$sedcom"s:$old:$new:g;"
		fi
	done < $namekey
	if [ $inplace -eq 1 ] ; then sed -i "$sedcom" "$targetfile"
	else sed "$sedcom" "$targetfile"
	fi

elif [ -d $targetfile ]; then 
	read -p "You are about to rename all the files in the directory \'$targetfile\'. Is this OK (y/N)? " confirm
	if [ $confirm == "y" ]; then
		while read line; do
			old=$(echo $line|awk  '{ print $2 }') ;# echo "Sample ID: $sampleno"
                	new=$(echo $line|awk  '{ print $1 }') ;# echo "Sample name is $sampname"
			for i in $targetfile/*; do
				if [[ $(echo $i) =~ "$old" ]]; then
					echo "Moving $i ${i/$old/${new}_$old}" 
					mv $i ${i/$old/${new}_$old}
				fi
			done

			
		done <$namekey
	fi



fi
