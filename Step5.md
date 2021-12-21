# Step 5: running PCAngsd

folders:

WGS/my_pcangsd/persite

Run similar Angsd code for all individuals as before, but less strict (-setMinDepthInd 3 instead of 4) as follows:

```bash
./angsd -b /users/c/p/cpetak/WGS/all_rmdups_noout.txt # all individuals except for the 3 outliers
-ref /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-anc /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-out /users/c/p/cpetak/WGS/angsd_noout/allpopsdepth3_angsd_polysites 
-nThreads 16 
-remove_bads 1 
-C 50 
-baq 1 
-minMapQ 30 
-minQ 20 
-minInd 119 
-setMinDepthInd 3 
-skipTriallelic 1 
-GL 1 
-doCounts 1 
-doMajorMinor 1 
-doMaf 1 
-doGlf 2 
-SNP_pval 1e-6
```

Followed by pcangsd as before:

```bash
python /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/my_pcangsd/persite/allpopsdepth3_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/my_pcangsd/persite/PCangsd_selection -selection -sites_save
```

TODO repeat with -pcadapt instead of -selection

Still running...

Tutorial I was following: http://www.popgen.dk/software/index.php/PCAngsdTutorial