gen_cor_dir=/scratch/groups/pritch/courtrun/hla/gen_cor/
output=/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/gencor/results/combined_gen_cor_sumstats.txt

cd $gen_cor_dir
# add header (can be from any existing gen_cor.log file)
awk '/p1/{{print $0}}' ${gen_cor_dir}AB1_BACT_INTEST_OTH-AB1_BACT_INTEST_OTH_genetic_cor.log > ${output}

for file in *
do
        #echo $file
        # to get the actual sum stats
        awk '/gz  /{print $0}' $file >> ${output} # select the line of interest bc it is the only one with "gz  " in it
done
