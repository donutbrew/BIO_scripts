#!/bin/bash
#DATASETS_B226-RAB2/genus.bar.SL148717.B226-RAB2.plotting/SnRa.Nairovirus.bar.SL148717.B226-RAB2.223019527.Nairovirus.1e-15.Blastn - identity
#                                                         bar.SL148717.Nairovirus.B226-RAB2.78191716.1e-15.Report - name #reads
# OUTPUT        genus.bar.SL148717.B226-RAB2.Blastn.fasta/SnRa.Nairovirus.bar.SL148717.B226-RAB2.78191716.Nairovirus.1e-15.Blastn.uniq.ex.fa - gi of winner

persum=0
#ls *.config || exit 1
root=$(ls *.fastq|cut -f1 -d'.')|| root=*
outfile=${1:-"${root}_surpi_blast_summary.tsv"}
echo $outfile
tempdir=$(mktemp -d)
filter=0
filterfile=$HOME/surpi_tmp/genusfilterlist.txt

#for sampledir in OUTPUT_$root/genus.bar*; do
for sampledir in $(find OUTPUT_* -maxdepth 1 -xtype d -name "genus.bar.*.fasta") ; do #OUTPUT_*/genus.bar*; do  
	if [ ! -d $(readlink -f $sampledir) ]; then continue; fi # From when we started generating files that are called genus*
	sampleid=$(echo $sampledir|sed "s/.*genus.bar.\(.*\).*.Blastn.fasta/\1/"); #echo $sampleid; continue 
	for aln in $sampledir/*.fa; do
		persum=0; reads=0;
 		aln2=${aln##*/}
		gi=$(echo $aln2|cut -f6 -d'.')
		genus=$(echo $aln2|cut -f7 -d'.')
#               blastfile=DATASETS_$root/genus.bar.$sampleid.$root.plotting/${aln2%.uniq.ex.fa}
#		blastfile=DATASETS_*/genus.bar.${sampleid}.plotting/${aln2%.uniq.ex.fa}
                blastfile=DATASETS_*/genus.bar.${sampleid}.plotting/${aln2%.ex.fa}
#               reportfile=DATASETS_$root/genus.bar.$sampleid.$root.plotting/bar.$sampleid.$genus.$root.$gi.1e-15.Report
		reportfile=DATASETS_*/genus.bar.${sampleid}.plotting/bar.${sampleid%.*}.$genus.*.$gi.1e-15.Report
		reads=$(wc -l < $blastfile) 
		giplus=$(head -n 3 $reportfile|tail -n 1)
		covfrac="$(grep "Coverage in bp" $reportfile|awk '{print$5}')/$(grep "Reference sequence length" $reportfile|awk '{print$5}')"
		covperc=$(grep \%Coverage $reportfile|awk '{print$3}')
		covdepth=$(grep "Average depth of coverage" $reportfile|awk '{print$6}')


		
		average=$(awk '{sum += $3} END {print sum/NR}' $blastfile)		
		echo -e "$sampleid\t$genus\t$reads\t$covfrac\t$covperc\t$covdepth\t$average\t$giplus" >> $tempdir/$sampleid
				
	done
	
	sort -k5 -nr $tempdir/$sampleid > $tempdir/$sampleid.sorted
	rm $tempdir/$sampleid
done

echo -e "SampleID\tGenus\tNum_reads\tCov_frac\tCov_perc\tCov_depth\tAvg_BLAST_ID\tGI" > $outfile
if [ $filter -eq 1 ] ; then
	terms=$(tr '\n' '|' < $filterfile|sed 's;||;;g'|sed 's;|$;;g')
	for i in $tempdir/*; do cat $i; echo -en "\n" ; done | grep -vE "$terms" >> $outfile
else
	cat $tempdir/* >> $outfile
fi

rm -rf $tempdir

