# HLA_finngen

Scripts for HLA analysis in FinnGen data by Satu Strauz and Courtney Smith

Directory Overview
- Genetic_Correlations: Running LDSC genetic correlation, identifying pairs of non-redundant (rg < 0.95) traits
  Snakemake_gencor, filter_indeptraits.R, combine_gen_cor.sh, plot_gen_cor.R

- Enrichment_Analysis: Quantifying enrichment in HLA region
  prep_HLA_combine.R, prep_genome_hits.R, plot_enrich.R, combine_hitssumstats.R, manhattan_hitsacrosstraits.R

- Pleiotropy_Matrix: Generation of matrix of disease hits (summary stats of all HLA GWAS hits combined across all traits by all traits), creating heatmap
  Snakefile_fullhits, make_matrix_fullhits.R, make_heatmap.R, make_heatmap_wspacing

- Haplotype_Regression: Identifying haplotypes, clustering into haplotype groups, performing haplotype regression analysis and adjacent analyses
  haplo_blockdefining.R (generated list of snps in each block)
  hapclusterandregress.R (Defines haplotype groups and runs regression)
  hapclusterandregress_alltraits.R (Repeats about but for all traits)
  hap_regress_plotting.R (Makes the dendrogram plots and haplotype genotype heatmaps for fig 5)
  annotategenes.R (Makes annotate snps with genes for fig 5)
  regress_enrich.R (Haplotype Group Burden Analysis)
  comparecorrresults.R (Haplotype Regression Trait Pair Correlation Measure analyses for fig 6)
  clusterregadjallele.R (Regression analysis with adjusting for relevant classical HLA alleles (VIF < 5) for each block + Plots comparing the z-scores for regression analysis before and after adjustment)
  allele_regress.R (Allele regression analyses w/ two approaches: one allele at a time and with all alleles together (after removing colinear) + plot z-score heatmap results)
  compareregmethods.R (Compare different regression methods, generate data tables for supplement)
  regresults_traitadjtrait.R (Trait adj for trait regression analysis)
  plot_clusterallelecor.R, check_cluster_correlations.R (Makes plot for correlations between alleles and haplotype groups)

- Figure_Generation: Generation of figures to visualize results
  trait_categories.R (Makes barplot with breakdown of trait categories for the FinnGen traits/classification, and barplot manual classification of HLA-associated triats by pathophys then by organ block)
  genes_hlahits_clean.R (Makes gene plot for fig 1)
  dist_hitstraits.R (Makes ridgeline plot showing distribution by trait group)
  R10_HLA_figures.R (Makes hits binned into genes barplot, makes plot with horizontal LD lines)
