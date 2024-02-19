WIR_DIR=/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA

gawk -f ~/my_project/bin/fetch_tab.awk ${WIR_DIR}/household \
~/projects/richards/restricted/ukb-27449/scratch/old-storage/dataset/ukb27449_20688.tab | gzip -c > ${WIR_DIR}/ukb27449_20688_household.tab.gz

#gawk -f ~/my_project/bin/fetch_tab.awk ../data/data.liver \
#~/projects/richards/restricted/ukb-27449/scratch/old-storage/dataset/ukb27449_20688.tab | gzip -c > ${WIR_DIR}/ukb27449_20688_count.tab.gz
