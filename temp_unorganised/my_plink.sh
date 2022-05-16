# PLINK 1.9
# PLINK v1.90b6.25 64-bit (5 Mar 2022)
# Process vcf and select SNPs to do LD decay and heatmap on
# cp filtered_vcf from LFMM, rename as temp_vcf
cut -d" " -f2 temp_vcf > proc_temp_vcf

#awk -F " " '{print>$1}' proc_temp_vcf # This writes each line to a file named after the first column = chromosome
#Biggest chromosome: NW_022145612.1 -> let's try this first

sed -i 's/_//g' vcf_header_temp
tail -10000 proc_temp_vcf > temp_vcf2 # downsample as it is too big
cat vcf_header_temp temp_vcf2 > temp_vcf_ready
sed -i 's/NA/./g' temp_vcf_ready

########################################
# LD decay

# Let's try with 500 SNPs
head -500 temp_vcf_ready > vcf_input # they were all on the same chromosome, NW_022145612.1

./plink --vcf vcf_input --recode --allow-extra-chr --out temp_vcf_plink #now we have pad and map
./plink --file temp_vcf_plink --make-bed --allow-extra-chr --out afterQC

# Calculate r2 values (correlation between SNPs)
# limited window - just nearby SNPs - default behaviour
#./plink --bfile afterQC --r2 --allow-extra-chr --out resultLD1

# all pairs also below LD 0.2 (default threshold), ld-window-r2 0 means minimum LD to display is 0
#./plink --bfile afterQC --r2 --allow-extra-chr --ld-window-r2 0 --out resultLD2

# adjust the number of SNPs and inter-SNP distances for which you want to compute LD
./plink --bfile afterQC --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 10000 --ld-window-kb 500000 --out resultLD
# ld-window is max number of SNPs, ld-window-kb is max basepairs distance between SNPs

# R script to visualise decay of LD over distances -> LD_decay_vis.R

########################################
# Create heatmap
# it calculates LD regardless of chromosome, but not a problem, just check that it's all 1 chromosome

# LD saved in a matrix of numbers
# for neat heatmap
./plink --bfile afterQC --r2 square --allow-extra-chr --out resultLD_heat
awk -F "\t" '{print $2}' vcf_input > locs_info

# R script to visualise heatmap -> LDheatmap.R

########################################
# LD pruning

# Split file into multiple files for each chromosome
# do stuff below for each in paralell -> makes it faster -> didn't need to do it in the end

# LD pruning -  remove SNPs with high LD with each other (removes one from each pair, or sliding window)
# What I want: if distance between 2 SNPs is < 2500 AND LD > 0.2 (?) -> prune one of the two

# This can be achieved via two commands: --indep which prunes based on the variance inflation factor (VIF),
# which recursively removes SNPs within a sliding window; second, --indep-pairwise which is similar,
# except it is based only on pairwise genotypic correlation

# The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold.
# The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.
# That is, this considers the correlations between SNPs but also between linear combinations of SNPs. A VIF of 10 is often taken to represent near
# collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). A VIF of 1 would imply that the SNP is completely independent
# of all other SNPs. Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the
# window size is too large, too many SNPs may be removed.

# The second generates the same output files as the first version; the only difference is that a simple pairwise threshold is used.
# The first two parameters (50 and 5) are the same as above (window size and step); the third parameter represents the r^2 threshold.
# Note: this represents the pairwise SNP-SNP metric now, not the multiple correlation coefficient; also note, this is based on the genotypic correlation,
# i.e. it does not involve phasing.

# To give a concrete example: the command above that specifies 50 5 0.5 would a) consider a window of 50 SNPs, b) calculate LD between each pair of SNPs
# in the window, b) remove one of a pair of SNPs if the LD is greater than 0.5, c) shift the window 5 SNPs forward and repeat the procedure.

# Example:
plink --file data --indep-pairwise 50 5 0.5

# Mine:
cat vcf_header_temp proc_temp_vcf > temp_vcf_ready # no downsampling this time
sed -i 's/NA/./g' temp_vcf_ready

./plink --vcf temp_vcf_ready --recode --allow-extra-chr --out temp_vcf_plink
./plink --file temp_vcf_plink --make-bed --allow-extra-chr --set-missing-var-ids @:#[b37]\$1,\$2 --out afterQC # my SNPs didn't have IDs

./plink --bfile afterQC --allow-extra-chr --indep-pairwise 50 5 0.5 # test this first, see how many prune
# -> Pruning complete.  215,911 of 991,430 variants removed. # this seems ok... 20%

# Also tried:
./plink --bfile afterQC --allow-extra-chr --indep-pairwise 50 5 0.2
# -> Pruning complete.  393,874 of 991,430 variants removed... 39.7%

# The output of either of these commands is two lists of SNPs: those that are pruned out and those that are not.
# A separate command using the --extract or --exclude option is necessary to actually perform the pruning.
# Each is a simple list of SNP IDs; both these files can subsequently be specified as the argument for a --extract or --exclude command.

# keep only selected markers for your data
#./plink --bfile afterQC --allow-extra-chr --exclude plink.prune.out --make-bed --out prunedSet

# convert back to vcf while doing that
./plink --bfile afterQC --allow-extra-chr --exclude plink.prune.out --recode vcf --out prunedvcf

# pruned vcf = prunedvcf.vcf

cp -r LFMM LFMM_pruned
#rm all stuff not needed, only keep .R, .sh, and cp here pruncedvcf.vcf_input
# repeat LFMM steps
sed -i '/^#/d' prunedvcf.vcf
# --> about same number of qval and pval outliers...
# --> increase strictness of pruning, eg 0.2??
# --> that basically didn't change results either
