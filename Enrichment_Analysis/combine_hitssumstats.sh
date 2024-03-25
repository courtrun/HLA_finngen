### Courtney Smith - Enrichment Analysis - Prep genome sumstats hits files
### Started on 12-7-2023
### Goal: Make a file that has a list of all traits w/ hit outside HLA, combine sumstats for the hits

# All traits with hit MAF > 1% and p < 1e-6 in genome outside HLA and the sumstats of these hits
cd /oak/stanford/groups/pritch/users/strausz/finngen_R10_gwas_hits/filt/

rm temp.txt
ls *.out | sed 's/^\.\///' >> temp.txt

while read p; do
    lines=$(wc -l < "$p")
    if [ "$lines" -gt 1 ]; then
        tail -n +2 "$p" >> all_hits_sumstats.txt
    fi
done < temp.txt

rm temp.txt
