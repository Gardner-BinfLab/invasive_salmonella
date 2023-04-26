#!/bin/bash

sample=$1
echo $sample > strain

# convert the assemblies to simulated reads
fastaq to_perfect_reads $sample reads 150 1 5 100
# map the simulated reads to the genes of interest in the ML model
bowtie2 --sensitive-local -x topgenes -q reads > bowtieout.sam
samtools view -b -S bowtieout.sam -t topgenes_nt.fasta.fai > samtoolsout.bam
samtools sort samtoolsout.bam -o samtoolsout.sort.bam
samtools index samtoolsout.sort.bam

bcftools mpileup -BQ0 -a DP,ADF,ADR -A -x -f topgenes_nt.fasta samtoolsout.sort.bam > samtoolsout.mpileup
bcftools call -P 0.001 -O b -A -M -S sample.ploidy -c samtoolsout.mpileup > bcftools.bcf
bcftools view  -i  'MIN(FMT/DP)>0 & MIN(INFO/MQ)>15' -c 1 bcftools.bcf > bcftools.vcf

# check for Ns in the assembly for our genes of interest
cat topgenes | while read i; do grep $i bowtieout.sam | grep "NNNNN" | wc -l >> NsinAssembly; done

# identify genes to be skipped in the sequence editing phase due to very few bases being covered - these genes will come up as NA
cat topgenes | while read i; do grep $i samtoolsout.mpileup | grep "DP=" | wc -l >> coverage; done
awk '$1 < 100 {print NR}' coverage > missinglines
cat missinglines | while read i; do sed -n ${i}p topgenes >> missinggenes; done

./vcfandfasta_to_edited.pl bcftools.vcf topgenes_nt.fasta sample_edited
translate6frames.sh frames=1 in=sample_edited.fna out=sample.faa overwrite=true

./trimtruncations.pl sample.faa
hmmsearch --domtblout hmmsearch models.hmm sample.faa.trim > tmp.txt
./parse_hmmsearch.pl hmmsearch

R CMD BATCH --no-save invasiveness_index.R
cat votes.csv >> votelog

rm reads sample.faa hmmsearch hmmsearch.parse sample.faa.trim strain tmp.txt votes.csv bowtieout.sam samtoolsout.bam samtoolsout.sort.bam* samtoolsout.mpileup NsinAssembly bcftools.bcf bcftools.vcf sample_edited.fna coverage missinggenes missinglines
