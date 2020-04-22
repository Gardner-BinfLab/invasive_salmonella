#!/bin/bash

sample=$1
echo $sample > strain
echo $sample

bowtie2 --sensitive-local -x topgenes_nt -1 ${sample}_1.fastq.gz -2 ${sample}_2.fastq.gz -p 20 > bowtieout.sam
samtools view -b -S bowtieout.sam -t topgenes_nt.fasta.fai > samtoolsout.bam
samtools sort samtoolsout.bam -o samtoolsout.sort.bam
samtools index samtoolsout.sort.bam

bcftools mpileup -BQ0 -a DP,ADF,ADR -A -x -f topgenes_nt.fasta samtoolsout.sort.bam > samtoolsout.mpileup
bcftools call -P 0.001 -O b -A -M -S sample.ploidy -c samtoolsout.mpileup > bcftools.bcf
bcftools view  -i 'MAX(FMT/DP)<1' -c 1 bcftools.bcf > missing.vcf
bcftools view -i 'MIN(FMT/DP)>5 & MIN(INFO/MQ)>15' -c 1 bcftools.bcf > bcftools.vcf

cat topgenes | while read i; do grep $i samtoolsout.mpileup | grep -P 'DP=' | wc -l >> coverage; done
awk '$1 < 100 {print NR}' coverage > missinglines
cat missinglines | while read i; do sed -n ${i}p topgenes >> missinggenes; done

./vcfandfasta_to_edited.pl bcftools.vcf topgenes_nt.fasta sample_edited
translate6frames.sh frames=1 in=sample_edited.fna out=sample.faa overwrite=true

./trimtruncations.pl sample.faa
hmmsearch --domtblout hmmsearch models.hmm sample.faa.trim > tmp.txt
./parse_hmmsearch.pl hmmsearch

R CMD BATCH --no-save invasiveness_index.R
cat votes.csv >> votelog

rm hmmsearch hmmsearch.parse sample.faa sample.faa.trim votes.csv missinggenes missing.vcf strain coverage missinglines tmp.txt invasiveness_index.Rout
rm bowtieout.sam samtoolsout.bam samtoolsout.sort* samtoolsout.mpileup bcftools.bcf bcftools.vcf sample_edited.fna
