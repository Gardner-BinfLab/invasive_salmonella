#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::Tools::CodonTable;

# want to take a vcf file and the genbank file corresponding to the genome used for mapping, to produce a FASTA file of all of the protein coding genes from the reference genome, altered according to the snps in the vcf file

# usage: vcfandgff_to_fasta.pl snps.vcf genome.gff outprefix

my $vcf = shift @ARGV;
my %snps;

open VCF, $vcf or die "Can't open VCF file $vcf";
while (<VCF>) {
    next if ($_ =~ /^#/);
        my @split = split /\t/;
    next if ($split[4] =~ /,/);
    # snps{chromosome}{position} = (reference, alternate)
    @{$snps{$split[0]}{$split[1]}} = ($split[3], $split[4]);
}

#print Dumper(@{$snps{AM933172}{'3484'}});
close VCF;

my $gff = shift @ARGV;
my %cds;
my %genome;
# have to read through all the way through, to gather CDS info then the sequence
open GFF, $gff;
my $gen = 0;
my $contig;
while (<GFF>) {
    if ($_ =~ /^>(\S+)\s/) {
        $gen = 1;
        $contig = $1;
        next;
    }
    if ($gen==0) {
        my @split = split /\t/;
        next if ($_ !~ /CDS/);
        my $id = "missing";
        if ($_ =~ /locus_tag=(\w+\d+\w?);/) {
            $id = $1;
        } elsif ($_ =~ /gene=(\w+);/) {
            $id = $1;
        } elsif ($_ =~ /locus_tag=(\S+);/) {
            $id = $1;
        }
        # cds{locus_id} = (contig, start, end, strand)
        @{$cds{$id}} = ($split[0], $split[3], $split[4], $split[6]);
    }
    if ($gen==1) {
        chomp;
        push @{$genome{$contig}}, split("", $_);
    }
}

#print Dumper(@{$cds{SEN0033}});
close GFF;

my $out = shift @ARGV;
open FAA, "> $out.faa";
open FNA, "> $out.fna";
open LOG, "> $out.log";

my @genes = keys(%cds);
print $#genes, "\n";
#print Dumper(sort{ $b <=> $a }keys($snps{CU458896}));

#print Dumper(@{$genome{$cds{SEN0033}[0]}}[$cds{SEN0033}[1]..$cds{SEN0033}[2]]);

foreach my $seq (keys(%cds)) {
#    print $seq, " ";
    my @fasta;
    next if (!defined($cds{$seq}[3]));      # look into this later
    @fasta = @{$genome{$cds{$seq}[0]}}[($cds{$seq}[1]-1)..($cds{$seq}[2]-1)];
    if (defined($snps{$cds{$seq}[0]})) {
        foreach my $pos (sort{ $b <=> $a }keys($snps{$cds{$seq}[0]})) { # for each snp in the array
            if ($pos >= $cds{$seq}[1] & $pos <= $cds{$seq}[2]) { # if it falls between the coordinates of the CDS
                #            print Dumper(@{$snps{$cds{$seq}[0]}{$pos}});
                splice @fasta, $pos-($cds{$seq}[1]), length($snps{$cds{$seq}[0]}{$pos}[0]) or die "$seq $pos-($cds{$seq}[1])-1 $#fasta $cds{$seq}[1]-1)..($cds{$seq}[2]-1"; # take out the old nucleotides
                splice @fasta, $pos-($cds{$seq}[1]), 0, split("",$snps{$cds{$seq}[0]}{$pos}[1]); # sub in the new ones
                print LOG $seq, "\t", $pos, "\t", $snps{$cds{$seq}[0]}{$pos}[1], "\n";
            }
        }
    }
    if ($cds{$seq}[3] eq "+") {
        print FNA ">$seq\n", join("", @fasta), "\n";
        my $seq_obj = Bio::Seq->new(-seq => join("", @fasta), alphabet => 'dna');
        my $prot = $seq_obj->translate(-complete => 1);
        print FAA ">$seq\n", $prot->seq, "\n";
    } else {
#        print "HERE\n";
        my $seq_obj = Bio::Seq->new(-seq => join("", @fasta), alphabet => 'dna');
        my $so2 = $seq_obj->revcom;
        my $prot = $so2->translate(-complete => 1);     #-complete => 1 tells you if you have a disrupted CDS
        print FNA ">$seq\n", $so2->seq, "\n";
        print FAA ">$seq\n", $prot->seq, "\n";
    }
}
close FAA;

