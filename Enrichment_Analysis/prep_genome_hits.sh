### Courtney Smith - Enrichment Analysis - Prep genome hits files
### Started on 3-5-2023
### Goal: Make a file that has a list of all traits w/ hit outside HLA, count number of hits

# Filtering to finemapped hits that are MAF > 1% and p < 1e-6

cd /oak/stanford/groups/pritch/users/strausz/finngen_R10_gwas_hits
mkdir filt

# Get it to work on all that start with the disease tag
rm temp.txt
ls | sed 's/^\.\///' >> temp.txt

while read p; do
awk -F'\t' 'NR==1 || ($10 <= 1e-6 && $20 < 0.99 && $20 > 0.01){print $2,$3,$6,$7,$8,$9,$18,$19,$10,$20}' $p >> filt/$p
done <temp.txt

rm temp.txt

cd /oak/stanford/groups/pritch/users/strausz/finngen_R10_gwas_hits/filt/

rm temp.txt
ls *.out | sed 's/^\.\///' >> temp.txt

echo "trait non_HLA_hits" > all_hits.txt

while read p; do
lines=$(wc -l < "$p")
if [ "$lines" -gt 1 ]; then
	echo "${p%%.*} $(tail -n +2 "$p" | wc -l)">> all_hits.txt
fi
done <temp.txt

rm temp.txt
