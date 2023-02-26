### Courtney Smith - Local Genetic Correlation Analysis - Make pathway file
### Started on 2-25-2023

# Define the header row
header="CHROM\tPOS_start\tPOS_end\tPathway\tGene\tPathway_Group"

# Define the first row of data
data="chr6\t28477897\t33448355\tHLA_REGION\tHLA\tHLA"

# Write the header row and first row of data to a new file
echo -e "$header\n$data" > /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt

# Add REACTOME pathway w/ boundaries extended 100kb
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\tREACTOME"}' /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION" >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\tREACTOME"}' /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt

# Add just class 2 HLA genes at gene boundaries
grep "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "HLA-DR" | awk '{print $1"\t"$2+100000"\t"$3-100000"\tCLASS_II\t"$4"\tHLA"}' >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt
grep "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "HLA-DP" | awk '{print $1"\t"$2+100000"\t"$3-100000"\tCLASS_II\t"$4"\tHLA"}' >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt
grep "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "HLA-DQ" | awk '{print $1"\t"$2+100000"\t"$3-100000"\tCLASS_II\t"$4"\tHLA"}' >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt

# Add just class 1 HLA genes at gene boundaries
grep "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION" /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "HLA-A" | awk '{print $1"\t"$2+100000"\t"$3-100000"\tCLASS_I\t"$4"\tHLA"}' >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt
grep "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION" /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "HLA-B" | awk '{print $1"\t"$2+100000"\t"$3-100000"\tCLASS_I\t"$4"\tHLA"}' >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt
grep "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION" /oak/stanford/groups/pritch/users/courtrun/general_data/allpathways_msig_within100kbOfGenes_wpathways.bed | grep "HLA-C" | awk '{print $1"\t"$2+100000"\t"$3-100000"\tCLASS_I\t"$4"\tHLA"}' >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt
