#!/usr/bin/env perl
#
# Script to extract reads from a Kraken-dekined lineage
# Takes a taxid, builds descendent list, and extracts fastq for that lineage
#
# Usage $0 <taxid> <database> <kraken-out> <fastq-original> 
#
# Author: Clint Paden
# Modified date: 2017-02-22


my $findtaxid = $ARGV[0];
my $db = $ARGV[1];
my $koutfile = $ARGV[2];
my $fastqinfile = $ARGV[3];

my $nodesfile = "$db/taxonomy/nodes.dmp";
my $namesfile = "$db/taxonomy/names.dmp";

my $follow = 1; # whether or not to continue down the tax tree

my %deslist;
my %namelist;

open my $nodelist, $nodesfile || die "Cannot open nodes.dmp";

while (<$nodelist>) {
	my @line = split "\t";
	push(@{$deslist{$line[2]}}, $line[0]); 
}

close $nodelist;

my @tilist;
push @tilist, $findtaxid;

if ( $follow eq 1 ) {
	my @descend = getLeaves(\%deslist, $findtaxid);
	push @tilist, @descend;
}

print STDERR join " ", @tilist ; print STDERR "\n";


open my $kout, $koutfile || die "Cannot open kraken out file $koutfile";
my %readlist;
while (<$kout>) {
	my @line = split "\t", $_;
	if ( grep( /^$line[2]$/, @tilist) || $findtaxid == 0) {
		$readlist{"$line[1]"} = "$line[2]";
	}
		
}
close my $kout;

open my $names, $namesfile || die "Cannot open names.dmp";
while (<$names>) {
	my @line = split "\t", $_;
	$namelist{$line[0]} = "$line[2]";
}

open my $fastq, $fastqinfile;
while (<$fastq>) {
	chomp;
	my @line = split " ", $_;
	my $rname = $line[0];
	$rname =~ /.(.*)/; $rname = $1;
	if (exists $readlist{$rname}) {
		print "$_";
		print " | $readlist{$rname} | $namelist{$readlist{$rname}}\n";
		for ( my $i = 0; $i < 3; $i++) {
			my $newline = <$fastq>;
			print $newline;
		}
	}
}





# Take in ref to array tree and list of taxid
# return the leaves below
sub getLeaves {
	my ($list, @taxid) = @_;
	my @thislist;
	for my $i (@taxid) {
		my @lookuplist;
		if (defined $$list{$i} ) {
			for my $j (keys $$list{$i}) { 
				push @thislist, $$list{$i}->[$j];
				push @lookuplist, $$list{$i}->[$j]
			}
		}
		if (defined $$list{$i}->[0]) {
			my @nextlist =  getLeaves($list, @lookuplist);
			push @thislist, @nextlist;
		}
	}
	
	return @thislist;
}
