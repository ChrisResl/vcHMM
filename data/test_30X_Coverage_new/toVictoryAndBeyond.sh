#!/bin/bash
for i in {1..100}
do

python vcHMM.py -i ref.fa -r read_$i.sam -o vcHMM_new_$i.vcf
done
