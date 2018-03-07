#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# usage: ./parse_bitscores.pl <orthlist> <folder with .search files>

print "Pulling together bitscores for all orthologous genes. This could take some time...\n";

my $orthlist = shift @ARGV;
my %orthologs;

open ORTH, $orthlist;
while (<ORTH>) {
    chomp;
    my @split = split /\",\"/;
    for (@split) {s/\"//g};
    for (@split) {s/___[0-9]*//g};
    for (@split) {s/\r//g};
    next if ($split[0] eq "Gene" || $split[3] ==1);
    foreach my $gene (14..$#split) {
        @{$orthologs{$split[$gene]}} = ($split[0], $gene-14);
    }
}
# print Dumper(\%orthologs);
close ORTH;

my $folder = shift @ARGV;
my @files = `ls $folder/*.search`;
my %bitscores;
my @strains;

foreach my $file (@files) {
    if ($file =~ /\/+(\S+).search/) {
        push @strains, $1;
    }
    open IN, $file;
    while (<IN>) {
        chomp;
        next if ($_ =~ /^#/);
        my @split = split /\s+/;
        for (@split) {s/\"//g};
        if (defined($orthologs{$split[0]}) && $split[6]<0.001) {
            $bitscores{$orthologs{$split[0]}[0]}{$split[3]}[$orthologs{$split[0]}[1]] = $split[7];
        }
    }
}
# print Dumper (\%{$bitscores{bcfC}});
close IN;

open OUT, "> bitscores.tsv";
open MODELS, "> models_used.tsv";
print OUT "\t", join("\t", @strains), "\n";
foreach my $gene (keys(%bitscores)) {
    my $bestscore = 0;
    my $bestmodel;
    foreach my $model (keys(%{$bitscores{$gene}})) {
        my $sum;
        map { $sum += $_ } grep {defined} @{$bitscores{$gene}{$model}};
       if ($sum > $bestscore) {
            $bestmodel = $model;
           $bestscore = $sum;
        }
    }
    print OUT $gene, "\t", join("\t", map { $_ // '' } @{$bitscores{$gene}{$bestmodel}}), "\n";
    print MODELS $gene, "\t", $bestmodel, "\n";
}
