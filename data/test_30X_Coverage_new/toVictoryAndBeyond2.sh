#!/bin/bash

for i in {51..100}
do

python vcHMM.py -i ref.fa -r read_$i.sam -o vcHMM_$i.vcf
done
