#!/usr/bin/env perl

use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;


my $m8file = $ARGV[0];
my $ctgfile = $ARGV[1];
my $searchtaxname = $ARGV[2];
my $searchtaxlevel = $ARGV[3];

print "Usage: $0 <m8file> <contig file> <Taxon name> <tax level>" if $ARGV[0] eq "-h";

open M8FILE, "<", $m8file || die "no M8 file\n";
open CTGFILE, "<", $ctgfile || die "no contig file\n";



my $tax = Bio::LITE::Taxonomy::NCBI->new (
                                               names=> "/scicomp/home/fep2/taxonomy/current/names.dmp",
                                               nodes=>"/scicomp/home/fep2/taxonomy/current/nodes.dmp"
					);
print STDERR "Done loading\n";
my $prevctg;
my %contigmatch;    

                                    
while (<M8FILE>) {
	chomp;
	next if /^#/ ;
	my @entry = split '\t', $_;
	next if $prevctg eq $entry[0]; # skip if this is the same contig and it's already been flagged
	my @gi = split '\|' , $entry[1];
	my $taxid = gi2tax($gi[1]) ;
	my $thistaxon = $tax->get_term_at_level($taxid,$searchtaxlevel);

	if ( $thistaxon =~ /$searchtaxname/i ) { 
		$contigmatch{"$entry[0]"} = [ "$gi[1]" ,  "$entry[2]"];
		#push @contigmatch, "$entry[0]" unless grep{$_ eq "$entry[0]"} @contigmatch; # Push contigname unless the value already exists
		$prevctg = $entry[0];
	}
	
}

my $firstline = 1;
my $printline = 0;

while (<CTGFILE>) {
	
	my $line=$_;
	
	if ( /^>/ ) {
#		$line =~ tr/>//d;
		my $key = $line =~ s/>(.*?)\s.*/$1/r;	# Some assemblies have additional info in header--cut before space or \n
		chomp $key;
		chomp $line; # must not chomp before regex match because I'm not spending any more time thinking about it.
		if ( defined $contigmatch{"$key"} ) {
			$printline = 1;
			$giname = `get_genbankfasta.pl $contigmatch{$key}->[0]|head -1|tr -d ">"|tr -d "\n"`;
			print "$line  | " . $contigmatch{$key}->[1] . " %ID | $giname\n" if $firstline == 1;
			print "\n$line  | " . $contigmatch{$key}->[1] . " %ID | $giname\n" if $firstline != 1;
			$firstline = 0;
		} else { $printline = 0; }
	} elsif ( $printline == 1 ) { chomp $line; print $line ; }
}

print "\n";


sub gi2tax($) {
	my $gi = $_[0];	
	my $start;
	my $stop;
	my $mid;
	my $word;
	my $fh;
	if ($#_ < 2) {
		$start = 0;
		$word = shift;
		open $fh, "/scicomp/home/fep2/taxonomy/161020/gi_taxid_prot.dmp";
		$stop = -s $fh;
		$mid = int(($start + $stop) / 2);
	} else {
		$start = shift;
		$stop = shift;
		$mid = int(($start + $stop) / 2);
		$word = shift;
		$fh = shift;
	}
	return 0 if $start == $mid or $stop == $mid;
	seek $fh, $mid, 0;
	<$fh>;
	my $line = <$fh>; chomp $line;
	($line, my $value) = split '\t', $line;
	return gi2tax($start, $mid, $word, $fh) if($word < $line);
	return gi2tax($mid, $stop, $word, $fh) if($word > $line);
	return $value;
}
