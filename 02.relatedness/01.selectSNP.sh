#!bin/bash

for i in {1..24}
do
awk -F "\t" '$6 > 0.05 && $6 < 0.95 && $8 > 0.8' /project/richards/restricted/ukb-general/scratch/genetic-data/genome/imputed.v3/bgen.stats/$i.txt >> /home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/SNPlist.txt 
done
