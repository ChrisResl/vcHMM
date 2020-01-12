#!/bin/bash

samtools faidx ref.fa

for i in {1..100}
do
wgsim -1 70 -2 70 -r 0.01 -R 0.15 -X 0.3 -e 0.0175 -N 5000 -S $i ref.fa read1_$i.fasta read2_$i.fasta 1 > variants_$i.txt

bwa index -p simref -a is ref.fa

bwa aln -t 4 simref read1_$i.fasta > read1_$i.sai
bwa aln -t 4 simref read2_$i.fasta > read2_$i.sai

bwa sampe simref read1_$i.sai read2_$i.sai read1_$i.fasta read2_$i.fasta > read_$i.sam


#samtools sort read_$i.sam

samtools view -bS read_$i.sam > read_$i.bam
samtools sort read_$i.sam > read_$i.bam

#FREEBAYES
freebayes -f ref.fa read_$i.bam > freebayes_$i.vcf

#SAMTOOLS
samtools mpileup -g -f ref.fa read_$i.bam > samtools_$i.bcf
bcftools call -c -v samtools_$i.bcf > samtool_$i.vcf

#VARSCAN
samtools mpileup -f ref.fa read_$i.bam > varscan_$i.mpileup
java -jar /Users/chris/bin/varscan/VarScan.v2.3.9.jar  mpileup2snp varscan_$i.mpileup --output-vcf 1 > varscan_$i.vcf


done
