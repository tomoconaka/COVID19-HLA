
path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/

plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs \
--allow-no-sex \
--keep $path/realcouples/WBcouples.sample \
--maf 0.05 \
--geno 0.01 \
--chr 6 \
--from-kb 28477 \
--to-kb 33449 \
--export A \
--out $path/realcouples/WB.couples_HLA

plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs \
--allow-no-sex \
--keep $path/randompairs/WBcouples.sample \
--maf 0.05 \
--geno 0.01 \
--chr 6 \
--from-kb 28477 \
--to-kb 33449 \
--export A \
--out $path/randompairs/WB.couples_HLA

for sample in randompairs realcouples
do
for chr in {1..22}
do
plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs \
--allow-no-sex \
--keep $path/${sample}/WBcouples.sample \
--indep-pairwise 5000 kb 5 0.1 \
--maf 0.05 \
--geno 0.01 \
--chr ${chr} \
--out $path/${sample}/WB.couples_${chr}

plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs \
--allow-no-sex \
--keep $path/${sample}/WBcouples.sample \
--extract $path/${sample}/WB.couples_${chr}.prune.in \
--export A \
--out $path/${sample}/WB.couples_${chr}

done
done
