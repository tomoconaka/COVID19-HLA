
path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/


plink2 --bgen /project/richards/restricted/ukb-general/scratch/genetic-data/genome/imputed.v3/bgen/9.bgen 'ref-first' \
--sample /project/richards/restricted/ukb-27449/scratch/old-storage/full_release/v3/imputed/w27449_20180503/sample/1-22.sample \
--keep $path/WB.couples \
--snps rs505922, rs8176719, rs8176746 \
--recode A --out $path/WB.couples.ABO
