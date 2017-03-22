#!/bin/bash

echo -e "Sample_ID\tGenus\tGI\tNumber_reads\tCovered/Total\tCoverage_Percentage\tAvg_Depth"
for genusdir in genus.bar.*; do
	if [[ $(ls $genusdir/*.coverage.top|wc -l) -ne 1 ]]; then echo "*.coverage.top is > 1 for $genusdir" >&2 ; continue; fi
	
	while read toplist; do

		gi=$(echo $toplist|cut -f4 -d' ')
		sid=$(echo $toplist|cut -f2 -d' ')
		genusid=$(echo $toplist|cut -f3 -d' ')
		reportfile=$genusdir/*$gi*.Report
		covfrac="$(grep "Coverage in bp" $reportfile|awk '{print$5}')/$(grep "Reference sequence length" $reportfile|awk '{print$5}')"
		covperc=$(grep \%Coverage $reportfile|awk '{print$3}')
		codepth=$(grep "Average depth of coverage" $reportfile|awk '{print$6}')
		readn=$(grep "Number of reads" $reportfile|awk '{print$8}')
		echo -e "$sid\t$genusid\t$gi\t$readn\t$covfrac\t$covperc\t$codepth"
	done < $genusdir/*.report.coverage.top
	
done
