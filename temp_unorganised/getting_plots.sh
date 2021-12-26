#!/bin/sh

dir=cpetak@vacc-user1.uvm.edu:WGS/angsd_new/pairwise_fst/outs_folded

while read line ; do
        pop1=$(cut -d '.' -f1 <<< $line)
        pop2=$(cut -d '.' -f2 <<< $line)
        echo $pop1
        echo $pop2
	scp $dir/res_${pop1}_${pop2}/manhattan_th.png man_plots/${pop1}_${pop2}_manplot.png
        sleep 0.5
done < $1
