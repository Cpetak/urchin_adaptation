# Step 8: Per-site nucleotide diversity

Folders: WGS/angsd_new/persite_diversity

Copied in all correct saf files, so LOM, TER BOD, KIB, SAN from angsd_new and CAP and FOG from angsd_noout into /persite_diversity

```
find . -name '*lapino*' -exec bash -c ' mv $0 ${0/\lapino/div}' {} \;
```

I followed this tutorial: http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests

Then, calculating sfs for each pop, reusing code from doing it per population pair (Step 6)

```bash
./realSFS ${dir}/${line}_angsd_allsites.saf.idx -P 16 -fold 1 > ${dir}/${line}_allsites.sfs
```

Then, calling saf2theta

```bash
./realSFS saf2theta ${dir}/${line}_angsd_allsites.saf.idx -sfs ${dir}/${line}_allsites.sfs -outname ${dir}/${line}_thetaout_div
```

Then, to calculate Tajimas D

```bash
./thetaStat do_stat ${dir}/${line}_thetaout_div.thetas.idx
.thetaStat do_stat out.thetas.idx -win 100 -step 100  -outnames theta.thetasWindow.gz
```

then repeated but with -win 1000 -step 1000 -> eg BOD_1000.thetasWindow.gz.pestPG

To clean the pestPG files:

```python
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import sys

for pop in ["LOM","TER","BOD","KIB","CAP","FOG"]:
    filename=pop + "_1000.thetasWindow.gz.pestPG"
    print(filename)
    df=pd.read_csv(filename,sep="\t", names=["window","chrome","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites"])
    df = df.iloc[1: , :]
    print("read file")
    df[['indicies','withdata','wind']] = df['window'].str.split('\)\(',expand=True)
    df = df.replace('\(','', regex=True)
    df = df.replace('\)','', regex=True)
    df[['Start','Stop']] = df['withdata'].str.split(',',expand=True) #look up this, it could be that I need to use another column for start and stop
    print("cleaned columns")
    df["tP"]=pd.to_numeric(df["tP"])
    df["numSites"]=pd.to_numeric(df["numSites"])
    df["Start"]=pd.to_numeric(df["Start"])
    df["Stop"]=pd.to_numeric(df["Stop"])
    df["tP/nSites"]=df["tP"]/df["numSites"]
    newdf=df[["chrome","Start","Stop","tP/nSites","wincenter"]]
    print("created new df")
    outfilename=pop + "_1000_cleaned.thetasWindow.gz.pestPG"
    newdf.to_csv(outfilename)
```

Then to make a figure:

```python
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import sys

mychr=str(sys.argv[1])
mypos=int(sys.argv[2])
sight_bef=int(sys.argv[3])
sight_aft=int(sys.argv[4])

print(mychr)
print(mypos)

for pop in ["LOM","TER","BOD","KIB","CAP","FOG"]:
    filename=pop + "_cleaned.thetasWindow.gz.pestPG"
    print(filename)
    newdf=pd.read_csv(filename, names=["chrome","Start","Stop","tP/nSites","wincenter"])
    newdf = newdf.iloc[1:,:]
    newdf["Start"]=pd.to_numeric(newdf["Start"])
    newdf["Stop"]=pd.to_numeric(newdf["Stop"])
    newdf["tP/nSites"]=pd.to_numeric(newdf["tP/nSites"])
    newdf["wincenter"]=pd.to_numeric(newdf["wincenter"])
    print("read file")

    print(newdf.head())

    mywindow = newdf[(newdf["chrome"]==mychr) & (newdf["Stop"]>=mypos) & (newdf["Start"]<mypos)]
    print(mywindow)

    if len(mywindow) < 1:
        print("ERROR: no position in dataset")
        abort

    before=newdf[(newdf["chrome"]==mychr) & (newdf["Stop"]<mypos)].tail(sight_bef)
    after=newdf[(newdf["chrome"]==mychr) & (newdf["Start"]>mypos)].head(sight_aft)

    #this is the output
    pi_bef=before["tP/nSites"].to_numpy()
    pi_aft=after["tP/nSites"].to_numpy()
    pi_mywindow=mywindow["tP/nSites"].to_numpy()
    pos_bef=before["wincenter"].to_numpy()
    pos_aft=after["wincenter"].to_numpy()
    pos_mywindow=mywindow["wincenter"].to_numpy()
    positions=np.concatenate((pos_bef, pos_mywindow,pos_aft), axis=None)
    pies=np.concatenate((pi_bef, pi_mywindow,pi_aft), axis=None)
    plt.xticks(np.arange(0, len(positions)+1, len(positions)/10))
    plt.axvline(x=len(before),color="red")
    plt.plot(positions,pies)
    plt.xlabel("Position in the genome, window (100bp) centers")
    plt.ylabel("Theta-Pi")
    figname="figures/"+mychr + "_"+ str(mypos) + "_"+str(sight_bef) +"_"+ str(sight_aft)+"_"+pop+".jpg"
    plt.savefig(figname)
```

