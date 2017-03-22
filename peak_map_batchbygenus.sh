#!/bin/bash

genus=${1:-Betacoronavirus}
outdir=${2:-peak_maps}

if [ $outdir != "-l" ]; then
	mkdir -p $outdir
fi

echo "Generating maps for $genus. " >&2
homedir=$PWD


for j in OUTPUT_*/genus.bar*.fasta; do 	
		fafile=$(ls $j/SnRa.$genus.* 2>/dev/null)
		if [ $? -ne 0 ]; then echo "No $genus for $j" >&2; continue; fi
		fafile=${fafile##*/}
		gi=$(echo $fafile|cut -f6 -d'.')
		barcode=$(echo $fafile|cut -f4 -d'.')
		if [ $outdir == "-l" ]; then 
			ls $homedir/DATASETS_*/genus.bar.$barcode*plotting/*$gi*.map
		else
			echo "found: $(ls $homedir/DATASETS_*/genus.bar.$barcode*plotting/*$gi*.map)" >&2
			peak_map2.sh $homedir/DATASETS_*/genus.bar.$barcode*plotting/*$gi*.map > $outdir/$genus.$barcode.$gi.tsv &
		fi
done

if [ $outdir != "-l" ]; then
	echo "All map files submitted. Processing . . ." >&2
	wait
	echo "Finished generating maps!" >&2
fi
