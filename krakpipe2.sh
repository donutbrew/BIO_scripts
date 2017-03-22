#!/bin/bash
#
# PDD Standardized read sanitizer and classification pipeline using Kraken
# Default Trim B2 adapter, Illumina adapter, Quality trim to Q20, Dust, Min Length 50
# 
# Requires cutadapt and prinseq-lite
#
# Usage: krakpipe.sh -1 <R1> [-2 <R2>] -o <output-location>
#	 Use without options for the grahipcal interface (coming soon!)
#
# Output 
#
# requires xfig, gnuplot
#
# Author: Clint Paden     Modified: 08 Apr 2015
# 

#set -x
###### Enter defaults here

# Load modules for aspen
module load cutadapt prinseq kraken krona

krakendb=/scicomp/groups/OID/NCIRD/DVD/GRVLB/pdd/db/kraken/defaulth
#krakendb=/ramdisk/defaulth
threads=30
cleanram=0


min_q=20
min_len=50

rawreadn1=""; rawreadn2=""; creadn1=""; creadn2="";deletetrim=0

totalcpu=$(grep -c ^processor /proc/cpuinfo)
proc_count=$(mktemp); echo 1 > $proc_count # Process counter. Init to 1.

read1=""; read2=""; root="KRAKPIPE_$(date +%y%m%d)"
wpid="KLSJFGNDFKLGNJSLD" # Easy initialization. Just go with it.
noramdisk=0
force=0
adeno=0
skipstep=""


# Try and clean up the best I can on exit...
trap 'ec=$? ; if [ "$ec" -ne 0 ] && [ "$ec" -ne 2 ] ; then echo DEAD; elif [ "$ec" -eq 2 ] ; then echo "Goodbye.\nProgram is exiting" ; fi ; rm $proc_count; jobs -p ; kill $(jobs -p) ; ' EXIT


# wait_proc accepts a variable or comma-delimited set of processes
wait_proc() {
	for i in ${1//,/ }; do 
		if [ -n "$2" ]; then echo "Waiting on $2"; fi
		while ps axc|grep -qP "\b$i\b" ; do sleep 10s; done 
	done
}

usage () {
	echo "Usage: $0 -d input_dir -o output_dir [ -b (kdaken db dir) -z (zip output files)] [-k (delete fastqfiles)] [-t threads] [-r root_name] [-q min_quality] [-l min_length] [-m (do not load database in ramdisk)] [-n (clean up ramdisk after run)] [-x <trim> (skip step)]" 
	exit 1
}



################################################################################
if [ -n "$1" ] ; then  #check for CMD line opts, otherwise go to GUI 

	cmdline=1 # indicates that we are in cmdline mode
	while getopts "d:o:r:q:l:t:zkx:anmfvb:" opt; do
		case $opt in
			d)	batchdir=$(readlink -m $OPTARG) ;;
			o)	outdir=$(readlink -m $OPTARG) ;;
			r)	root=$OPTARG ;;
			q)	min_q=$OPTARG ;;
			l)	min_len=$OPTARG ;;
			z)	zipout=1 ;;
			k)	deletetrim=1 ;;
			t)	totalcpu=$OPTARG ;;
			m)	noramdisk=1 ;;
			n)	cleanram=1 ;;
			x)	skipstep=$OPTARG ;;
			f)	force=1 ;; # overwrite if folder is there.
			a)	adeno=1 ;; # Custom behavior to remove AdV reads...Temporary
			b)	krakendb=$OPTARG ;;
			v)	set -x ;;
			*)	usage	;;
		 esac
	done
	shift $((OPTIND-1)) 

outdir=${outdir:-"$batchdir/$root"}

if [ $force -ne 1 ]; then
	a=1; while [ -d "$outdir" ]; do echo "directory $outdir exists"; if [ $a -eq 1 ] ; then outdir=${outdir}_$a; else outdir=${outdir%_*}_$a; fi; ((a++)); done #simple to prevent overwrites  *****maybe generate a subdirectry
fi
mkdir -p $outdir
echo "Output directory will be $outdir"


### IMPLEMENT ZENITY INTERFACE HERE############################################# 
elif [ -z "$1" ] ; then

	zenity --info --text "GUI not yet finished. Sorry"
	exit 0;
fi
################################################################################



## Mount a ramdisk for kraken and load the database
if [ "$noramdisk" -ne 1 ] ; then
	# **later check to make sure there is enough ram	
	if ! mount |grep ramdisk >/dev/null ; then sudo mount -t tmpfs tmpfs /ramdisk/ -o size=110G; if [ $? -ne 0 ] ; then echo "ERROR MOUNTING RAMDISK" >&2; exit 1; fi ; fi
	krakenramdb=/ramdisk/${krakendb##*/}
	if [ ! -d "$krakenramdb" ] ; then
		echo "Copying database into RAM"
		touch /ramdisk/testfile # make sure you can write to it
		if [ "$?" -eq 0 ] ; then 
			(rsync -a $krakendb /ramdisk && touch /ramdisk/${krakendb##*/}/goodcopy ; chmod -R a+rwx /ramdisk/${krakendb##*/} ) & ramcopy="$!" # Copy over db and track pid to wait later.
			echo "Copying database into RAM"
			krakendb="$krakenramdb"
			rm /ramdisk/testfile
		else 
			echo "Unable to use RAMDISK. Using databas from disk ( $krakendb )."
		fi 
	else echo "Database already in RAM!"	
	fi
	krakendb=/ramdisk/${krakendb##*/} #set var to new location.
fi
## Copy over the kraken database





numfiles=$(( $(ls "$batchdir/"*.fastq* |wc -l|cut -f1 -d' ') / 2))
if [ $numfiles -lt $totalcpu ]; then threads=$((totalcpu / numfiles)); maxproc=$numfiles; # If the number of files is less than the number of cpus, then we may assisgn multiple threads to a process. maxproc downs't really matter here.
else threads=1; maxproc=$totalcpu ; # If we have a ton of files, only allow one thread per file (pair) and put all cpus to work at once.
fi

mkdir $outdir/fastq_raw #Directory to copy raw fastq
cd $batchdir                 			#NEED TO REMOVE THIS AND ADD IT INTO THE FOR LOOP. IT IS MESSING UP RELATIVE PATHS
for rfile in *.fastq* ; do 
	while [ "$(cat $proc_count)" -gt "$maxproc" ] ; do sleep 10s; done #2016-07-08 changed -ge to  -gt and if to while

	
	if [[ "$rfile" =~ "_R1" || "$rfile" =~ "Clean_1" ]] ; then    # Get Read 1
		read1=$rfile; echo Read 1 is $read1
		if [[ "$rfile" =~ "Clean_1" ]]; then skipstep="trim" ; fi # detect cleaned files automatically and skip additional cleaning
	
	elif [[ "$rfile" =~ "_R2" || "$rfile" =~ "Clean_2" ]] ; then  # Continue. Probably should put some error checking here to make sure the names are the same
	
		echo $(($(cat $proc_count) +1)) > $proc_count # Count the total processes using a file. Increment here.
						    	      # Only matters when there are more than total cpu anyway
	(	
		read2=$rfile; echo Read 2 is $read2


		# Start counting reads
			
		
		( echo Copying $read1 and $read2
		  cp $read1 $outdir/fastq_raw/ || echo "Problem with read1 file $read1" >&2 &
		  cp $read2 $outdir/fastq_raw/ || echo "Problem copying read2 file $read2" >&2 
		  echo -e "$(date)\tUnzipping $read1 and $read2"
		  wait # Must wait on both files to transfer
		  if [[ $read1 =~ ".gz" ]]; then
			gunzip $outdir/fastq_raw/$read1 || echo "Problem unzipping $read1" >&2 & 
		  	gunzip $outdir/fastq_raw/$read2 || echo "Problem unzipping $ $read2" >&2
		  fi
		  wait # wait for subshell procs to finish
		) # Subshell to make copies complete simultaneously before moving on
		
		
		read1=${read1%.gz}; read1=$outdir/fastq_raw/${read1##*/}  # Change read file locations to the fastq_raw
		read2=${read2%.gz}; read2=$outdir/fastq_raw/${read2##*/}
		
		echo "Raw reads" >> "$outdir"/$root.log
		rawreadn1=$(( $(cat $read1 | wc -l | cut -f1 -d' ') / 4))
		rawreadn2=$(( $(cat $read2 | wc -l | cut -f1 -d' ') / 4))
		echo -e "$read1\t$rawreadn1" >> $outdir/$root.log 
		echo -e "$read2\t$rawreadn2" >> $outdir/$root.log
		echo 


		readbasen=${read1%%_R1*}; readbasen=${readbasen##*/} #generate base name without extension or abs directory
		
		if [ "$skipstep" != "trim" ]; then
			pdd_trim.sh -p -1 $read1 -2 $read2 -o "$outdir/$readbasen" &>>$outdir/$root.trim.log
		
		
			cread1="$outdir/${readbasen}_Clean_1.fastq"
			cread2="$outdir/${readbasen}_Clean_2.fastq"
		else 
			cread1="${read1}"
			cread2="${read2}"
		fi
		# Count filtered reads
		(	
			echo "Trimmed, filtered, reads"
			creadn1=$(( $(wc -l $cread1 | cut -f1 -d' ') / 4))
			creadn2=$(( $(wc -l $cread2 | cut -f1 -d' ') / 4))
			echo -e "$cread1\t$creadn1" >> $outdir/$root.log 
			echo -e "$cread2\t$creadn2" >> $outdir/$root.log
			echo 


		) & wpid="$!" 


#		wait_proc $wpid  # this can be moved if we dont have to pass the raw value tp kp_c


		wait_proc $ramcopy "RAMdisk copywait"
		if [ ! -f $krakendb/goodcopy ] && [ $noramdisk -ne 1 ] ; 
			then echo "Database copy error. Try again, or with the -m flag to avoid the RAMdisk.">&2; 			
			rm -rf /ramdisk/$krakendb 
			kill $$
		fi
		echo -e "$(date)\tRAMDISK copy complete"
		
		# makes a bunch of files-- $root.tmp is CR-delim file
		if [ $adeno -eq 1 ] ; then
		echo -e "$(date)\tAdV subtraction on for $cread1"
		krakpipe_classify.sh -1 $cread1 -2 $cread2 -o $outdir -d $krakendb -t $threads -r $root -w $rawreadn1 -c 10508 -l Adenovirus &>> $outdir/$root.err
		else
		echo -e "$(date)\tDefault classification for $cread1"
		krakpipe_classify.sh -1 $cread1 -2 $cread2 -o $outdir -d $krakendb -t $threads -r $root -w $rawreadn1 &>> $outdir/$root.err
		fi

		#add kraken data into log
		#fi
		
		
		echo $(($(cat $proc_count) -1)) > $proc_count # Count the total processes by using a file. Decrement here.
	) & wpid="$!,$wpid" 
	fi
	
done




echo -e "$(date)\tFinishing kraken classification..."

wait




echo "Generating plots..."

ktImportTaxonomy -o "$outdir/$root.krona.html" "$outdir"/*.krona &>> $outdir/$root.err &

(head -n 1 $outdir/$root.raw.tmp && tail -n +2 $outdir/$root.raw.tmp|sort) > $outdir/$root.raw.percent
rm $outdir/*.raw.tmp
(head -n 1 $outdir/$root.tmp && tail -n +2 $outdir/$root.tmp|sort)> $outdir/$root.percent
rm $outdir/*.tmp
#cat $outdir/*raw.tmp > $outdir/$root.raw.percent
#rm $outdir/*.raw.tmp
#cat $outdir/*.tmp > $outdir/$root.percent
#rm $outdir/*.tmp
(head -n 1 $outdir/$root.summary_report && tail -n +2 $outdir/$root.summary_report|sort) > $outdir/$root.tmp.summary_report ; mv $outdir/$root.tmp.summary_report $outdir/$root.summary.tsv


### Bar graph section
if [ $(wc -l < $outdir/$root.raw.percent) -lt 25 ]; then fontsize=9; else fontsize=7; fi
#echo "=stackcluster;PhiX;Virus;Bacteria;Archaea;Human;Unclassified;Bad_reads
echo "=stacked;PhiX;Virus;Bacteria;Archaea;Human;Unclassified;Bad_reads
max=100
min=0
=nogridy
=noupperright
legendx=right
legendy=center
#=nolegoutline
yscale=1.6
xscale=1
yformat=%g%%
xlabel=Sample
ylabel=Percentage of total reads
rotateby=-45
colors=black,red,yellow,light_green,med_blue,grey5,grey7
#leading_space_mul=1
=nocommas
font=Helvetica
fontsz=$fontsize
legendfontsz=9
intra_space_mul=0.01
barwidth=1



=table"  > $outdir/$root.perf
#multimulti=Clean Only" > $outdir/$root.perf


#cat $outdir/$root.percent >> $outdir/$root.perf
#if [ $numfiles -eq 1 ] ; then echo -ne "\n__\t0\t0\t0\t0\t0\t0" >> $outdir/$root.perf ; fi
#echo -e "\nmultimulti=Total" >> $outdir/$root.perf
cat $outdir/$root.raw.percent >> $outdir/$root.perf
#if [ $numfiles -eq 1 ] ; then echo -e "\n__\t0\t0\t0\t0\t0\t0" >> $outdir/$root.perf ; fi


#make bar graphs
bargraph.pl -pdf $outdir/$root.perf > "$outdir/${root}_ReadDistribution.pdf"

wait

# Cleanup

rm *.krona
mkdir "$outdir"/kraken_files
mv "$outdir"/*.kraken* "$outdir"/kraken_files
rm $proc_count

rm -rf $$outdir/fastq_raw $outdir/*.summary_report

if [ "$deletetrim" -eq 1 ] ; then
	rm "$outdir"/*Clean*.fastq
fi

if [ $cleanram -eq 1 ]; then rm -rf /ramdisk/defaulth ; fi

echo -e "$(date)\tDone!"



