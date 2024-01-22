# selecting 4123 (10% of the 41235 snps in HLA region w/ MAF > 1%)
cat /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/HLA_reg_snps_maf_ID.txt | head -n 1 > /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_HLAhitsand10percentbackground_refalt.txt
cat /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/HLA_reg_snps_maf_ID.txt | tail -n +2 | shuf -n 4123 >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_HLAhitsand10percentbackground_refalt.txt

# Adding inthe 428 hits
cat /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_HLAhits_refalt.txt >> /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_HLAhitsand10percentbackground_refalt.txt

# Sort and keep unique (4517 snps)
mv /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_HLAhitsand10percentbackground_refalt.txt temp.txt
{ cat temp.txt | head -n 1; cat temp.txt | tail -n +2 | sort | uniq; } > /oak/stanford/groups/pritch/users/strausz/finngen_R10_sumstats/hits/R10_HLAhitsand10percentbackground_refalt.txt
rm temp.txt
