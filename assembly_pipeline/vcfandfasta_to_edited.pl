#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

# want to take a vcf file and the genbank file corresponding to the genome used for mapping, to produce a FASTA file of all of the protein coding genes from the reference genome, altered according to the snps in the vcf file

# usage: vcfandgff_to_fasta.pl snps.vcf genome.gff out_prefix

open MISS, "missinggenes";
my @miss;
while (<MISS>) {
    chomp;
    push @miss,$_;
}

my $vcf = shift @ARGV;
my %snps;

open VCF, $vcf or die "Can't open VCF file $vcf";
while (<VCF>) {
    next if ($_ =~ /^#/);
    my @split = split /\t/;
    next if ($split[4] =~ /,/);
    # snps{sequence}{position} = (reference, alternate)
    @{$snps{$split[0]}{$split[1]}} = ($split[3], $split[4]);
}

#print Dumper($snps);
close VCF;

my $fasta = shift @ARGV;
my %genome;
open FASTA, $fasta;
my $seq;
while (<FASTA>) {
    if ($_ =~ /^>(\S+)\s/) {
        $seq = $1;
        next;
    } else {
        chomp;
        my @split = split /\t/;
        push @{$genome{$seq}}, split("", $_);
    }
}

close FASTA;
#print Dumper(\@{$genome{STM0018}});

my $out = shift @ARGV;
open FNA, "> $out.fna";

my @genes = keys(%genome);

foreach my $seq (@genes) {
    next if ($seq ~~ @miss);
    my @fasta = @{$genome{$seq}};
    if (defined($snps{$seq})) {
        foreach my $pos (sort{ $b <=> $a } keys %{$snps{$seq}}) { # for each snp in the array
            # splice ARRAY,OFFSET,LENGTH,LIST
            # Removes the elements designated by OFFSET and LENGTH from an array, and replaces them with the elements of LIST
            splice @fasta, $pos-1, length($snps{$seq}{$pos}[0]), split("",$snps{$seq}{$pos}[1]) or die; # sub the ref for the alt
        }
    }
    print FNA ">$seq\n", join("", @fasta), "\n";
}
close FNA;

