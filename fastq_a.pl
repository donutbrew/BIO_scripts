#!/usr/bin/env perl
#
# Script to convert fastq files into fasta. 
# Usage: $0 [-f use only '@'' as header id] reads.fastq > reads.fasta
#
# Author: Clint Paden

my $force = 0;
if ($ARGV[0] eq "-f") { $force = 1; shift; } # use -f id the headers are not consistent. This can cause problems id q scores start with @
my $fastqfile = $ARGV[0];

open my $fastq, $fastqfile;
my $type = "first";
my $header = "@";

while (<$fastq>) {
	if ($type eq "first" && $force == 1) { $header = substr($_, 0, 4) ; } 
	if (/^$header/) {
		print "\n" unless $type eq "first";
		print ">";
		print substr($_,1);
		$type = "seq";		
	} elsif (/^\+/) {
		$type = "qual";
	} elsif ( $type eq "seq" ) {
		chomp;
		print $_;
	} 
	
	if (eof) { print "\n";}

}
