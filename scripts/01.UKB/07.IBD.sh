
path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/
plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs  \
--allow-no-sex \
--keep $path/WB.couples \
--maf 0.01 \
--mac 5 \
--geno 0.05 \
--hwe 1e-6 \
--exclude range ~/my_project/repo/BQC19_genotype_pipeline/data/LdRegion-AbecasisHg19.txt \
--snps-only just-acgt \
--indep-pairwise 1000 kb 5 0.1 \
--out $path/WB.couples

plink --bfile ~/09.COVID19/src/01.UKBB/01.GWAS/data/ukb_cal_allChrs  \
--allow-no-sex \
--keep $path/WB.couples \
--extract $path/WB.couples.prune.in \
--make-bed \
--out $path/WB.couples

plink2 --bfile $path/WB.couples \
--make-king-table \
--king-table-subset $path/WBcouples.king.new \
--memory 1000 \
--out $path/WBcouples.king
