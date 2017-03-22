#!/usr/bin/env perl 



my %genusList;
my %readlist;
my %lineage;
my $minreads = $ARGV[1] // 10 ;
my $type = $ARGV[2] // "freq";

use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;

#if ($ARGV[0] =~ "-h" ) { print "Usage: $0 <list-of-name/taxid> <surpi-coverage-summary-file>\n" ; exit 1; }
if ($ARGV[0] =~ "-h" ) { print "Usage: $0 <list-of genus.counttbale> <min reads> <freq/reads\n" ; exit 1; }

open(FILELIST, "<", $ARGV[0]) || die "Cant open file list";

my $tax = Bio::LITE::Taxonomy::NCBI->new (
                                               names=> "/scicomp/home/fep2/taxonomy/current/names.dmp",
                                               nodes=>"/scicomp/home/fep2/taxonomy/current/nodes.dmp"
                                             );

#my @taxonomy = $tax->get_taxonomy_from_name('Betacoronavirus');
#	print join('\t',@taxonomy);

print "Count\tLineage\n";
while (<FILELIST>) {
	chomp;
	my $thisfile = $_;
	my @header;
	open(FILE, "<", $thisfile);
	while (<FILE>) {
		my @line = split '\t';
		if ( $line[0] =~ "Genus" ) { @header = @line; }
		else {
			for (my $i=2; $i <= $#line ; $i++ ) {
				if ( $line[$i] >= $minreads) {
					$line[0] =~ s/\s+$//; $line[0] =~ s/_/ /;
					$line[$i] =~ s/\s+$//;
#					if ( $line[0] =~ Flavivirus ) { print "$header[$i]   $line[$i]   $#line\n"; }
					$genuslist{"$line[0]"} = $genuslist{"$line[0]"} + 1;
					$readslist{"$line[0]"} += $line[$i];
				}		
			}
		}


	}
	close("FILE");

}

foreach my $key (keys(%genuslist)) {
#	my @taxonomy = $tax->get_taxonomy_from_name($key);
	my @taxonomy = getStdTaxonomy($key);
#	my $taxid = $tax->get_taxid_from_name($key);
#	my @taxonomy = $tax->get_taxonomy_with_levels($taxid);
	print "$genuslist{$key}\t" if ($type eq "freq") ;
	print "$readslist{$key}\t" if ($type eq "reads");
	print join("\t",@taxonomy);
	print "\n";
}




# Gets Kingdom-ss/-fam-genus-species
sub getStdTaxonomy {
	my $ingenus = $_[0];
	my $taxid = $tax->get_taxid_from_name($ingenus);
	my @temptax = $tax->get_taxonomy($taxid);
	my $super = "$temptax[0]";
	my $class; my $subclass;
	if ( $temptax[1] =~ "ssRNA" ) { 
		$subclass = "$temptax[2]"; $class = $temptax[1]; }
		else {$class = "$temptax[1]"; $subclass = '-'; }
	my $order = $tax->get_term_at_level($taxid,"order");
	my $family = $tax->get_term_at_level($taxid,"family");
	my $subfamily = $tax->get_term_at_level($taxid,"subfamily");
	my $genus = $tax->get_term_at_level($taxid,"genus");
 
	my @stdtax = ($super, $class, $subclass, $order, $family, $subfamily, $genus);
	for my $i (@stdtax) { $i =~ s/undef/-/; }
	return @stdtax;	

}

=pod

# Takes taxonomy array of arrays and gets the particular taxonomy
sub getTaxonByLevel  {
	my ($taxonlevel, @thistaxonomy) = @_;
	my $result = "undef";
	foreach my $i (@thistaxonomy) {
		my $z = $i->[1] ;
		if ("$z" eq "$taxonlevel") { 
			$result = $i->[0];
			last;
		}

	}
	return $result ;
}

# Gets the genuses
sub getSpecificLevel {
	my ($taxonlevel, @thistaxonomy) = @_;
	my $result = 0;
#	my $result = $thistaxonomy[$#thistaxonomy]->[0]; #default
	for my $i (@thistaxonomy) {
		if ($i->[1] eq $taxonlevel) {
			$result = "$i->[0]";
			last;
		}
	}
	return $result;


}

# Shows the level of taxonomy
sub getLevelofTaxonomy(@) {
	my @thistaxonomy = @_;
	return $thistaxonomy[$#thistaxonomy]->[1];
}
sub getLevelofTaxid($) {
	my $thistaxid = $_[0];
	my $level = "undef";
	my @taxonomy = $tax->get_taxonomy_with_levels($thistaxid);
	$level = getLevelofTaxonomy(@taxonomy);
	return $level ;
}

sub taxid2name($){

	my $taxid = $_[0];
	my @tax1 = $tax->get_taxonomy($taxid);
	return $tax1[$#tax1];



=cut


