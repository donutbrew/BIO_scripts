#!/usr/bin/env perl
#

#use strict;
use Bio::SeqIO;
my $outfile1;
my $outfile2;
my $filter = 0;


if ($ARGV[0] == "-f") { 
	shift;
	$outfile1 = shift;
	$outfile2 = shift;
	$filter = 1;
}

my $seqio = Bio::SeqIO->new(-file => "$ARGV[0]", -format => "fasta") || die "could not access file $ARGV[0]!\n";

open(my $good, ">", $outfile1);
open(my $bad, ">", $outfile2);

my $line = 1 ;
while(my $seq = $seqio->next_seq) { 
#	if ($seq->id =~ />/ ) { 
#		print "OK";
#	}
        my $temp = $seq->seq;
	my $tempid = $seq->id;
	$temp =~ s/\r//g;

	if ( $tempid =~ /[,:">]/ or $tempid =~ /^\*/ ) { 
		print "Invalid header: $tempid\n";
		if ( $filter == 1 ) { print $bad ">$tempid\n$temp\n"; }
	}

	elsif  ($seq->seq =~ /[^ACTGactgNMRWSYKVHDBnmrwsykvhdb]/ ) {
		print "Sequence \"$tempid\" invalid: $temp\n";
		if ( $filter == 1 ) { print $bad ">$tempid\n$temp\n"; }
	}
	elsif ($temp =~ /^\s*$/) { 
		print "Empty sequence: $tempid\n"; 
		if ( $filter == 1 ) { print $bad ">$tempid\n$temp\n"; }
	}
	else { 
		 if ( $filter == 1 ) { print $good ">$tempid\n$temp\n"; }
	}
	$line++ ; 
		 
}
