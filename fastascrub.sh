#!/bin/bash

if [ $1 == "-i" ] ; then 

file=$2; 
temp=$(mktemp -p $PWD) 

awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' $file > $temp
mv $temp $file

else file=$1; 



awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' $file

fi
