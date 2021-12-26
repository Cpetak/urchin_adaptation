# Hardcalling specific regions

Folders: WGS/hard_called_genotypes

for all:

```bash
bcftools mpileup -f $ref $f | bcftools call -mv -Ob -o ${f}_all_data.bcf
bcftools view ${f}_all_data.bcf > ${f}_all_data.vcf
```

for specific region:

```bash
bcftools mpileup -r CHR:FROM-TO -f $ref $f | bcftools call -mv -Ob -o ${f}_spec_data.bcf
```

For example, carbonic anhydrase (LOC579101)

```bash
bcftools mpileup -r NW_022145605.1:5673513-5691489 -f $ref $f | bcftools call -mv -Ob -o ${f}_LOC579101.bcf
bcftools view ${f}_LOC579101.bcf > ${f}_LOC579101.vcf
#mv all .bcf and vcf to designated folder
```

but we need to call haplotypes!

and for that, we need to merge all individual vcfs into 1 vcf

```bash
spack load bcftools@1.10.2
for i in $(ls|grep .vcf); do bgzip $i; done
for i in $(ls|grep .vcf); do bcftools index $i; done
bcftools merge --merge all *.vcf.gz -O v -o merged.vcf
```

then:

vcf_phase.py --vcf merged.vcf --phase-algorithm beagle

TODO this output (phased) doesnt make any sense but merged.vcf is correct!

to visualise vcf files:

```python

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


def vcf2name(vcf):
    return "_".join(vcf.name.split("_")[:2])


gcolors = {"0/1": 1,"1/1": 2}
gcolors_inv = {v: k for k, v in gcolors.items()}
gcolors_inv[0]="N/A"
data = []
inds = list(Path("vcfs").iterdir())

names = [vcf2name(vcf) for vcf in inds]
sortedname2id = {name: i for i, name in enumerate(sorted(names))}
id2name = {v: k for k, v in sortedname2id.items()}

pops = []
for name in sorted(names):
    pop = name[:3]
    if pop not in pops:
        pops.append(pop)

for vcf in inds:
    with open(vcf, "r") as f:
        lines = f.readlines()
        lines = [line for line in lines if line[0] != "#"]

    ID = vcf2name(vcf)

    for line in lines:
        ch, pos, _, ref, alt, _, _, _, _, geno = line.split("\t")
        gtype = geno.split(":")[0]
        if len(ref) == 1 and len(alt) == 1:
            data.append((ID, int(pos), gcolors[gtype]))


max_pos = max([el[1] for el in data])
min_pos = min([el[1] for el in data])
min_pos = 12346747
print(min_pos, max_pos, max_pos - min_pos)

mat = np.zeros((len(inds), max_pos - min_pos + 1))

for name, pos, gtype in data:
    rel_pos = pos - min_pos
    if rel_pos >= 0:
        mat[sortedname2id[name], pos - min_pos] = gtype



reps = mat.shape[1] // mat.shape[0]
tmp = []
for row in mat:
    for i in range(reps):
        tmp.append(row)
tmp = np.array(tmp)


fig,ax = plt.subplots(figsize=(10,10))
# from https://stackoverflow.com/a/7230921
# make a color map of fixed colors
cmap = colors.ListedColormap(["black", "white", "red"])
bounds = [0, 1, 2, 3]
norm = colors.BoundaryNorm(bounds, cmap.N)
# tell imshow about color map so that only set colors are used
img = plt.imshow(tmp, interpolation="nearest", cmap=cmap, norm=norm)
labels = [gcolors_inv[i] for i in range(len(gcolors)+1)]
# from https://joseph-long.com/writing/colorbars/
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.2)
cbar = plt.colorbar(cmap=cmap, norm=norm, boundaries=bounds, ticks=np.arange(5) + 0.5, cax=cax)
cbar.ax.set_yticklabels(labels)

# separate individuals
for i in range(1,len(inds)):
    ax.axhline(reps*i-0.5,color="grey",lw=0.15)
# separate populations
for i in range(1,7):
    ax.axhline(20*reps*i-0.5,color="green",lw=0.5)

step = reps*20
yticks = step/2 + np.arange(7)*step
ax.set_yticks(yticks)
ax.set_yticklabels(pops)

xticks = np.array(ax.get_xticks())
ax.set_xticklabels([int(x) for x in (xticks + min_pos)])
ax.set_xlabel("Position")
plt.tight_layout()
plt.savefig("hardwork.pdf",dpi=100)
plt.show()
```

