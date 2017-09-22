#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# read in models and their orthologs for scoring then print out a table of scores for each representative

my %searchresults;
my %scores;

my @searchfiles = `ls */*.search_eggNOG`;
foreach my $search (@searchfiles) {
	chomp $search;
	&parse_hmmsearch_tbl($search);
}
#print Dumper (\%searchresults);

open IN, "orthologlist.tsv";
open OUT, "> eggNOGscores.tsv";
open CHECK, "> checkmodels";
my $count = 1;
while (<IN>) {
	my @split = split /\s+/;
	my $firstmodel;
    foreach my $pos (0..$#split) {      # get scores for each entry, and print out the name of the Typhimurium reference, if applicable
		my $entry = $split[$pos];
		if ($pos == 0 & $entry eq "NA") {
			print OUT "GENE$count-$firstmodel";
            $count++;
		}
		next if ($entry eq "NA");
		&get_domain_arch($entry);
		if ($pos == 0) {
			$firstmodel = $scores{$entry}[0];
			print OUT $entry, "_", $firstmodel;
		}
	}
	foreach my $entry (@split[0..$#split]) {
		if ($entry eq "NA") {
			print OUT "\tNA";
			print CHECK "\tNA";
		} else {
			print OUT "\t$scores{$entry}[1]";
			print CHECK "\t$scores{$entry}[0]";
		}
	}
	print OUT "\n";
	print CHECK "\n";
}

sub parse_hmmsearch_tbl {
	my ($tbl) = @_;
	open TBL, "<", $tbl;
	while(<TBL>){
		next if($_ =~ /^#/);
			chomp;
		my @splat = split /\s+/;
		#gene - domain, eval, score, start, end
		push @{$searchresults{$splat[0]}}, [$splat[3], $splat[7], $splat[6]];
	}
	close TBL;
}
sub get_domain_arch {		# takes the single top eggNOG hit for each gene
	my $identifier = shift @_;
	#	print Dumper (\@{$searchresults{$identifier}});
	if ($#{$searchresults{$identifier}}>=0) {
		my @sorted = sort {$a->[1] <=> $b->[1]} @{$searchresults{$identifier}}; #sort on bitscore
		$scores{$identifier} = @sorted[$#sorted];
	} else {
		$scores{$identifier} = ["NA","NA","NA"];
	}
}
