# Created by Courtney Smith - HLA - HESS Local Gen Cor
# Started on 1-11-2023
# Goal of script: Reformat to bed file with relevant columns
# Called by snakefile in step step3_to_bed

cat $1 | cut -f 1,2,3,6,4,7 | awk 'OFS="\t" {if (NR > 1) $1="chr"$1; print}' | sed '1d' > $2

