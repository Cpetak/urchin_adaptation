# Logbook
# Process vcf and make into plink format
cut -d" " -f2 temp_vcf > proc_temp_vcf
cat vcf_header_temp proc_temp_vcf > temp_vcf_ready
sed -i 's/_//g' temp_vcf_ready
sed -i 's/NA/./g' temp_vcf_ready
./plink --vcf temp_vcf_ready --recode --allow-extra-chr --out temp_vcf_plink #now we have pad and map
./plink --file temp_vcf_plink --make-bed --allow-extra-chr --out afterQC

# Calculate r2 values (correlation between SNPs)
# limited window - just nearby SNPs - default behaviour
./plink --bfile afterQC --r2 --allow-extra-chr --out resultLD1

# all pairs also below LD 0.2 (default threshold), ld-window-r2 0 means minimum LD to display is 0
./plink --bfile afterQC --r2 --allow-extra-chr --ld-window-r2 0 --out resultLD2

# adjust the number of SNPs and inter-SNP distances for which you want to compute LD
./plink --bfile afterQC --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 100 --ld-window-kb 2000 --out resultLD3
# ld-window is max number of SNPs, ld-window-kb is max basepairs distance between SNPs

# LD saved in a matrix of numbers
# for neat heatmap
./plink --bfile afterQC --r2 square --allow-extra-chr --out resultLD4
awk -F "\t" '{print $2}' temp_vcf_ready > locs_info

# R script to visualise decay of LD over distances -> LD_decay_vis.R
# R script to visualise heatmap -> LDheatmap.R

########################################

# LD pruning -  remove SNPs with high LD with each other (removes one from each pair)

# replace --nonfounders with --make-founders!
plink --bfile afterQC --cow --make-founders --indep-pairwise 50 5 0.7 --out afterQC
# keep only selected markers for your data
plink --bfile afterQC --cow --exclude afterQC.prune.out --make-bed --out prunedSet
