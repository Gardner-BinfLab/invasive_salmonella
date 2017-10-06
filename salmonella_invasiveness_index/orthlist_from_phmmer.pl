#!/usr/bin/perl

# for random forest model testing - first you run your key indicator genes across your proteomes (fetch_seqs.sh), then you go through the phhmer output and generate an orthlist and fetch sequences for all of the top hits

use warnings;
use strict;
use Data::Dumper;

my %orthlist;
my $column=0;

my @files = `ls phmmer/*`;
my @strains;
foreach my $file (@files) {
    open IN, $file;
    my $strain;
    if ($file =~ /phmmer\/(\S+).phmmer/) {
        push @strains, $1;
    }
    my @seen;
    while (<IN>) {
        chomp;
        next if ($_=~ /^#/);
            my @split = split /\s+/;
        next if ($split[2] ~~ @seen);
        $orthlist{$split[2]}[$column] = $split[0];
        push @seen, $split[2];
    }
    close IN;
    $column++;
}

open OUT, "> orthlist.tsv";
print OUT "Gene\t", join ("\t", @strains), "\n";
foreach my $gene (keys(%orthlist)) {
    print OUT $gene, "\t", join ("\t", map { $_ // '' } @{$orthlist{$gene}}), "\n";
}
