#!/bin/bash
# krakpipe classify module

# Modified 2015-08-03  NOT TESTED


# Defaults
krakdb=/ramdisk/defaulth
filter=0
customtax=0
custlabel="Archea"
batchdir=""
rawtot="-"
krakenpath=$HOME/bin/kraken
PATH="$krakenpath:$PATH"

readtot=0;unclass=0;humantot=0;bacteriatot=0;archaeatot=0;virustot=0;phixtot=0;phixtot=0; 


# wait_proc accepts a variable or comma-delimited set of processes
wait_proc() {
	for i in ${1//,/ }; do 
		while ps axc|grep -q $i ; do sleep 10s; done 
	done
}

usage () { echo "Usage: $0 -1 <read1> -2 <read2> -o <output-dir> -r <root-name> -d <database> [-t threads -f filter -z (zip) -b (batch) -q (quick) -w <# raw reads>] [-c <taxid> (Custom taxid to replace Archea) -l <label>]" }
}

while getopts "1:2:o:zd:b:t:r:fqw:c:l:" opt; do
	case $opt in
		1)	read1=$OPTARG ;;
		2)	read2=$OPTARG ;;
		o)	outdir=$OPTARG ;;
		r)	root=$OPTARG ;;
		b)	batchdir=$OPTARG ;;
		d)	krakdb=$OPTARG ;;
		t)	threads=$OPTARG ;;
		w)	rawtot=$OPTARG ;;
		f)	filter=1;;
		z)	zipin=1;;
		c)	customtax=$OPTARG ;;
		l)	custlabel=$OPTARG ;;
		q)	echo "quick  mode not supported yet";;
		*)	usage	;;
	 esac
done
shift $((OPTIND-1)) 

####### Check whether R1 / R2 / single-end / batch

if [ -n "$batchdir" ] ; then
	echo "BATCH MODE in $batchdir"
	echo "Not yet supported in standalone. "; exit 1;
	#assign R1 / R2
fi



mkdir -p "$outdir"
rootper=$root.${read1##*/}; rootper=${rootper%%_L00*}

krakenout="$rootper.krakenout"
echo "Starting the kraken..." >&2
				#Add logic for gzip conpressed $zipin and quick mode and paired
kraken --db $krakdb --threads $threads --fastq-input --output $outdir/$krakenout --paired $read1 $read2
if [ "$?" -ne 0 ]; then echo "Problem with Kraken on $read1" >&2; fi

if [ $filter -ne 1 ] ; then
	krakenreport="$rootper.krakenreport"
	kraken-report --db $krakdb "$outdir/$krakenout" > "$outdir/$krakenreport"
	awk -v OFS='\t' '{ print $2,$3 }' "$outdir/$krakenout" > "$outdir/$krakenout.krona"
elif [ $filter -eq 1] ; then
	krakenreport="$rootper.filter.krakenreport"
	kraken-filter --db $krakdb --threshold 0.2 "$outdir/$krakenout" > "$krakenout.filter"
	kraken-report --db $krakdb "$outdir/$krakenout.filter" > "$outdir/$krakenreport"
	krakenout="$krakenout.filter"
	awk -v OFS='\t' '{ print $2,$3,substr($5,3) }' "$outdir/$krakenout" > "$outdir/$krakenout.krona"
fi


#ktImportTaxonomy -o "$outdir/$krakenout.krona.html" "$outdir/$krakenout.krona" &>/dev/null &

samplename=${read1##*/}; samplename=${samplename%%_L0*}

readtot=$(wc -l $outdir/$krakenout|cut -f1 -d' ')
unclass=$(awk -F$'\t' '$5==0 { print $2 }' $outdir/$krakenreport); 
	unclassroot=$(awk -F$'\t' '$5==1 { print $3 }' $outdir/$krakenreport) ; unclass=$(($unclass + $unclassroot)); unclass=${unclass:-0}
humantot=$(awk -F$'\t' '$5==9606 { print $2 }' $outdir/$krakenreport); humantot=${humantot:-0}
bacteriatot=$(awk -F$'\t' '$5==2 { print $2 }' $outdir/$krakenreport); bacteriatot=${bacteriatot:-0}
if [ $customtax -eq 0 ]; then
	archaeatot=$(awk -F$'\t' '$5==2157 { print $2 }' $outdir/$krakenreport); archaeatot=${archaeatot:-0}
else
	archaeatot=$(awk -v mytax="$customtax" -F$'\t' '$5==mytax { print $2 }' $outdir/$krakenreport); archaeatot=${archaeatot:-0}
	#need to add archea to bacteria?
fi
virustot=$(awk -F$'\t' '$5==10239 { print $2 }' $outdir/$krakenreport); virustot=${virustot:-0}
phixtot=0		
phixtot=$(awk -F$'\t' '$5==374840 { print $2 }' $outdir/$krakenreport) ; phixtot=${phixtot:-0}
dataloc=$(pwd)
((realvirus=virustot-phixtot))   #Subtracts out PhiX from virus total
if [ $customtax -ne 0 ]; then ((realvirus=realvirus-archeatot)) ; fi # Subtracts out the custom tax from virus--this behavior may change, but I'm using it to get rid of AdV from 293 cells

if [ ! -f $outdir/$root.summary_report ]; then
	echo -e "Name	Raw_Reads	Filtered_Reads	Unclassified	Human	Bacteria	$custlabel	Virus	PhiX	File_read1	File_read2" >> $outdir/$root.summary_report
fi
echo -e "$samplename\t$rawtot\t$readtot\t$unclass\t$humantot\t$bacteriatot\t$archaeatot\t$realvirus\t$phixtot\t$read1\t$read2" >> $outdir/$root.summary_report

#echo -e "Files:\t$read1 $read2" >> $outdir/$root.summary_report
#echo -e "Read total:\t$readtot" >> $outdir/$root.summary_report
#echo -e "Unclassified:\t$unclass" >> $outdir/$root.summary_report
#echo -e "Human:\t$humantot" >> $outdir/$root.summary_report
#echo -e "Bacteria:\t$bacteriatot" >> $outdir/$root.summary_report
##echo -e "Viruses :\t$virustot" >> $outdir/$root.summary_report
#echo -e "Real virus:\t$realvirus" >> $outdir/$root.summary_report


#output as total of clean reads /percentage
# *.a.tmp is for the particular line, preprocessing--bc requires a CR. Adding tabs manually and remove CRs at the end and then concat the line into a comprehensive file--root.tmp
readtot2=$(echo "$unclass + $humantot + $bacteriatot + $archaeatot + $realvirus + $phixtot" | bc -l)
echo -en "$samplename\t"  >> $outdir/$rootper.a.tmp
echo  "$phixtot / $readtot2 * 100"     |bc -l >> $outdir/$rootper.a.tmp ; echo -en "\t" >> $outdir/$rootper.a.tmp
echo  "$realvirus /$readtot2 * 100"    |bc -l >> $outdir/$rootper.a.tmp ; echo -en "\t" >> $outdir/$rootper.a.tmp
echo  "$bacteriatot / $readtot2 * 100" |bc -l >> $outdir/$rootper.a.tmp ; echo -en "\t" >> $outdir/$rootper.a.tmp
echo  "$archaeatot / $readtot2 * 100"  |bc -l >> $outdir/$rootper.a.tmp ; echo -en "\t" >> $outdir/$rootper.a.tmp
echo  "$humantot / $readtot2 * 100"    |bc -l >> $outdir/$rootper.a.tmp ; echo -en "\t" >> $outdir/$rootper.a.tmp
echo  "$unclass / $readtot2 * 100"     |bc -l >> $outdir/$rootper.a.tmp ; echo -en "\t" >> $outdir/$rootper.a.tmp
echo  "0" >> $outdir/$rootper.a.tmp 

tr -d '\n'  < $outdir/$rootper.a.tmp >> $outdir/$root.tmp  # echo -n above causes error with bc, so we have to remove CRs here
echo -en "\n" >> $outdir/$root.tmp   # CR after this entry 
rm $outdir/$rootper.a.tmp

readtot3=$(echo "$readtot2 + ($rawtot-$readtot2)" | bc -l) 
echo -en "$samplename\t" >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "$phixtot / $readtot3 * 100"    |bc -l >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "$realvirus /$readtot3 * 100"   |bc -l >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "$bacteriatot / $readtot3 * 100"|bc -l >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "$archaeatot / $readtot3 * 100" |bc -l >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "$humantot / $readtot3 * 100"   |bc -l >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "$unclass / $readtot3 * 100"    |bc -l >> $outdir/$rootper.raw.a.tmp ; echo -en "\t" >> $outdir/$rootper.raw.a.tmp
echo  "($rawtot-$readtot2)/$readtot3 *100" |bc -l >> $outdir/$rootper.raw.a.tmp

tr -d '\n'  < $outdir/$rootper.raw.a.tmp >> $outdir/$root.raw.tmp
echo -en "\n" >> $outdir/$root.raw.tmp
rm $outdir/$rootper.raw.a.tmp





#sed -i 's/classified/classified as human/g' $filename
