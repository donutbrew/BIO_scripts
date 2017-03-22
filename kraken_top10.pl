#!/usr/bin/env perl
#
# 
# Script to extrac the top 10 hits from a list of kraken report files
# 
# Usage:  $0 [family/genus/species] <list of files>  > top10 file

use Switch;

my $taxlevel = $ARGV[0];
shift;

my $limit = 10;

switch ($taxlevel){
	case /family/i	{$taxlevel="F"}
	case /genus/i	{$taxlevel="G"}
	case /species/i	{$taxlevel="S"}
	else			{die "Please specify Family, Genus, or Species report\n\t"}
}

my @filelist = @ARGV;
if ( $#filelist < 1) {die "Need to specify kraken report files\n\t"; }


print "Sample	%Reads	ROLLUP	TAXON	RANK	TAXID	NAME\n";
for my $file (@filelist) {
	open my $filehandle, $file;
	my $samplename = $file =~ s/(\w+)_L001.*$/$1/r; #Filename for typical Illumina files
	#my $samplename = $file =~ /(\w+)_L001.*$/$1/r; #Filename for typical Illumina files
	my @thisfile;
	while (<$filehandle>){
		my @line =  split "\t";
		@line = map { join(' ',split(' '))} @line; # Gets rid of the spaces
		if ("$line[3]" eq "$taxlevel") {
				push @thisfile, [ @line ];
		}
	}
	my @sorted = sort { $b->[2] <=> $a->[2]} @thisfile; # Sort array in reverse order by num reads
	my $print;
	if ($#sorted < $limit) {  $print = $#sorted + 1; } else {  $print = $limit;}
	for (my $i=0; $i < $print; $i++) { 
		print "$samplename\t";
		print join "\t", @{$sorted[$i]}; print "\n";
	}
	close $filehandle;
	print "\n";

}