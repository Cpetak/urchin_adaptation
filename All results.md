# All results

## Sequencing/mapping quality check:

[Multiqc Report nicely displayed is available here](https://htmlpreview.github.io/?https://github.com/Cpetak/urchin_adaptation/blob/main/images/multiqc_report.html) 

[File with all mapping stats](data/all_mapping_stats.csv)
In all 3 images below, x axis is the 140 individuals

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/coverage_fig.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/flagstat_fig.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/mapping_stat.png" width="400" />

Histogram of average coverage per individual:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/hist_coverage.png" width="400" />

PCA of average coverage:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_cov.png" width="400" />
     
Again, no clustering is visible.

## Evidence for non population structure:

### PCA

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_1.png" width="400" />
	
3 individuals seem to be very different from the other 137 individuals. Thus, these 3 were dropped from further analysis. New PCA with 137 individuals (ANGSD was rerun with only 137 individuals):

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/PCA_2.png" width="400" />

No clustering by population can be seen. Holds to every combination of PCs, tested up to PC 5.

### Pairwise global Fst

2 methods used, they correlate with each other:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/realsfs_avefst.png" width="400" />

And pairwise global Fst does not correlate with distance regardless of method:

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/realsfs_dist.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/avefst_dist.png" width="400" />

### Bayenv matrix

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/bay_meanof5_corr.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/global_pairwise_fst.png" width="400" />

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/fst_vs_bayenv.png" width="400" />

Mantel test to compare Bayenv cov matrix with the global pairwise Fst matrix: with 10,000 permutations, p=0.0402, r=0.623.

### LFMM PCA

<img src="https://github.com/Cpetak/urchin_adaptation/blob/main/images/LFMM/PCA.jpg" width="400" />