
path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/

plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs \
--allow-no-sex \
--keep <(awk '{print $1, $1}' /home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt) \
--maf 0.01 \
--geno 0.01 \
--exclude range ~/my_project/repo/BQC19_genotype_pipeline/data/LdRegion-AbecasisHg19.txt \
--make-bed \
--out $path/tmp1

plink --bfile $path/tmp1 \
--indep-pairwise 5000 kb 5 0.1 --out $path/tmp2

plink --bfile $path/tmp1 --extract $path/tmp2.prune.in \
--make-bed \
--keep $path/WB.couples \
--out $path/WB.couples

king -b $path/WB.couples.bed --related --degree 13 --rplot --prefix $path/WB.couples-pruned-related-degree13

head -1 $path/WB.couples-pruned-related-degree13.kin0 > $path/WB.couples-pruned-related-degree13.kin0_filtered
awk '(FNR==NR){m[$2]=$2; n[$2]=$3; next}($2 in m && $4 == n[$2]){print $0}($4 in m && $2 == n[$4]){print $0}' \
$path/WBcoupleswithHLA $path/WB.couples-pruned-related-degree13.kin0 >> $path/WB.couples-pruned-related-degree13.kin0_filtered


