#!/usr/bin/perl
use strict;
use warnings;

# goes through CDS translations and trims anything after a stop codon

my $file = shift @ARGV;
open IN, $file;
open OUT, "> $file.trim";
while (<IN>) {
	if ($_ =~ /^(.*?)\*/) {
		print OUT $1, "\n";
	}
	else {
		print OUT $_;
	}
}
