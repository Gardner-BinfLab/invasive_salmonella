#!/usr/bin/perl

use warnings;
use strict;

my $input = shift @ARGV;
open IN, $input;
open OUT, "> $input.parse";
while (<IN>) {
    next if ($_=~ /^#/);
    $_ =~ s/^\#/@/g;            # a lot of the samples I'll work with use # as part of the ID
    my @split = split("\\s+");
    # print $#split, "\n";
    print OUT join("\t", @split[0..21]), "\t", join(" ", @split[22..$#split]), "\n";
}

close OUT;
close IN;
