### Courtney Smith - Local Genetic Correlation Analysis - Calculate local pathway values
### Started on 2-25-2023

ml biology
ml bedtools

cd /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/

pathwayfile=/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt

# Pull out the cumulative values for each pathway group
for path in `grep -v Pathway_Group $pathwayfile | cut -f 6 | sort | uniq`; do
    for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step3/*/step3.h2_snps.bed; do
        intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -w $path)  -wa | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12` -v path=$path '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
    done
done | tr ' ' '\t' >> pathway_group_rhog.txt

# Pull out the cumulative values for everywhere in genome
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step3/*/step3.h2_snps.bed; do
    awk -v trait=`echo $trait | cut -d'/' -f12` -v path=GENOME '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}' $trait
done | tr ' ' '\t' >> pathway_group_rhog.txt

# Pull out the cumulative values for each individual pathway
for path in `grep -v Pathway_Group $pathwayfile | cut -f 4 | sort | uniq`; do
    for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step3/*/step3.h2_snps.bed; do
        intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -w $path)  -wa | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12` -v path=$path '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
    done
done | tr ' ' '\t' >> pathway_group_rhog.txt

# Pull out the cumulative values for LD blocks covered by all pathways together in the pathway file
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step3/*/step3.h2_snps.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile)  -wa | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12` -v path=ALL '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_rhog.txt

# Pull out the cumulative values for everywhere in genome except the LD blocks included in the pathway file
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step3/*/step3.h2_snps.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile) -v -wa | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12` -v path=NONE '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_rhog.txt

# Pull out the cumulative values for everywhere in genome except the LD blocks included in the pathway file that are part of the pathway group HLA
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step3/*/step3.h2_snps.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -v REACTOME) -v -wa | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12` -v path=NONE_HLA '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_rhog.txt

###
pathwayfile=/oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/pathway.txt

# Pull out the cumulative values for each pathway group
for path in `grep -v Pathway_Group $pathwayfile | cut -f 6 | sort | uniq`; do
    for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait1.bed; do
        intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -w $path) -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f1` -v path=$path '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
    done
done | tr ' ' '\t' >> pathway_group_h2.txt

# Pull out the cumulative values for everywhere in genome
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait1.bed; do
    awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f1` -v path=GENOME '(NF == 6) {n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}' $trait
done | tr ' ' '\t' >> pathway_group_h2.txt

# Pull out the cumulative values for each individual pathway
for path in `grep -v Pathway_Group $pathwayfile | cut -f 4 | sort | uniq`; do
    for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait1.bed; do
        intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -w $path)  -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f1` -v path=$path '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
    done
done | tr ' ' '\t' >> pathway_group_h2.txt

# Pull out the cumulative values for LD blocks covered by all pathways together in the pathway file
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait1.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile)  -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f1` -v path=ALL '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_h2.txt

# Pull out the cumulative values for everywhere in genome except the LD blocks included in the pathway file
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait1.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile) -v -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f1` -v path=NONE '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_h2.txt

# Pull out the cumulative values for everywhere in genome except the LD blocks included in the pathway file
for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait1.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -v REACTOME) -v -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f1` -v path=NONE_HLA '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_h2.txt

# Repeat the above but for second trait
for path in `grep -v Pathway_Group $pathwayfile | cut -f 6 | sort | uniq`; do
    for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait2.bed; do
        intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -w $path) -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f2` -v path=$path '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
    done
done | tr ' ' '\t' >> pathway_group_h2.txt

for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait2.bed; do
    awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f2` -v path=GENOME '(NF == 6) {n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}' $trait
done | tr ' ' '\t' >> pathway_group_h2.txt

for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait2.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile)  -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f2` -v path=ALL '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_h2.txt

for path in `grep -v Pathway_Group $pathwayfile | cut -f 4 | sort | uniq`; do
    for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait2.bed; do
        intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -w $path)  -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f2` -v path=$path '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
    done
done | tr ' ' '\t' >> pathway_group_h2.txt

for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait2.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile) -v -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f2` -v path=NONE '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_h2.txt

for trait in /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hess_local_gc/hg19/step2/*/step2.h2_snps_trait2.bed; do
    intersectBed -a $trait -b <(grep -v Pathway_Group $pathwayfile | grep -v REACTOME) -v -wa | awk 'NF == 6' | sortBed -i | uniq | awk -v trait=`echo $trait | cut -d'/' -f12 | cut -d'-' -f2` -v path=NONE_HLA '{n+=$4; h2+=$5; se+=sqrt($6); i+=1} END {print trait, path, n, h2, se, i}'
done | tr ' ' '\t' >> pathway_group_h2.txt

# Make sure only one row for each trait for each pathway group (since some slight variations in calculation that are insignificant but can't just uniq out
cat pathway_group_h2.txt | sort | uniq | awk '!seen[$1,$2]++' | wc -l
cat pathway_group_h2.txt | cut -f1 | sort | uniq | wc -l
cat pathway_group_h2.txt | cut -f2 | sort | uniq | wc -l
mv pathway_group_h2.txt temp.txt
cat temp.txt | sort | uniq | awk '!seen[$1,$2]++' > pathway_group_h2.txt
rm temp.txt
