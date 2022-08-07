


## BD vs Controls analysis 
Running differential methylation analysis (DMA)-bipolar vs controls group using MethylKit package 

## Workflow of the pre-processing steps

- Excluding the samples from run #13 due to the low quality of the run. 
- Excluding the samples that have less than 93% or NA value for the bisulfite conversion rate.  
- Summarizing the technical replicates
- Identifying the cytosines that were covered in all samples, with at least have read coverage equal 1 
- Normalizing the read coverage using the median read coverage as scaling factor
- Filtering the CpG sites that have less than 10 reads, and CpG sites that have more than 99.9 percentile of coverage in each sample.

## Workflow of the pre-processing steps-continued

- Combine the reads for F and R strands for each CpG site. Now, we have the total methylation for each CpG site.
- Running PCA to check which principal components are associated with the batch ID.  
-	Correcting the data by removing the PC1 and PC2
- Conducting the differential methylation analysis 

## Installing the required pacakges

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylKit")

 library (methylKit)
 
## Reading the methylation calls files 
               
```{r echo= TRUE, results= "hide",message = FALSE  }
myobjgroup = methRead(file.list,
               sample.id = sample.id,
               header = TRUE, 
               sep = "\t", 
               assembly="hg38",
               treatment = group,
               context="CpG",resolution = "base", 
               mincov = 0)#reading the files 
```
## Filtering any CpG sites that have 0 read coverage 
``` {r echo= TRUE, results = "hide", message = FALSE}
filtered.nocoverage.myobjgroup=filterByCoverage(myobjgroup,lo.count=1,
                                           lo.perc=NULL,
                                           hi.count=NULL,hi.perc=NULL)
```
## Normalizing the read coverage 


``` {r echo= TRUE, results= "hide" , message = FALSE}
newObjgroup = normalizeCoverage(filtered.nocoverage.myobjgroup,
                                method="median")
```

## Filtering the low read coverage CpG sites 

Filtering the CpG sites that have less than 10 reads, and CpG sites that have more than 99.9 percentile of coverage in each sample.


``` {r echo= TRUE, results= "hide" ,message = FALSE }
filtered.newobjgroup = filterByCoverage(newObjgroup,lo.count=10,
                                lo.perc=NULL,
                                 hi.count=NULL,hi.perc=99.9)
```
## Combining the reads from the two strands

``` {r echo= TRUE, results="hide" , message = FALSE}
meth.filteredgroup <- methylKit::unite (filtered.newobjgroup,
                             destrand=TRUE )
```
## Running the PCA to detect which PC is associated with the batch ID

``` {r echo= TRUE, message = FALSE}
as.all.samplesgroup=assocComp(mBase=meth.filteredgroup, metadata)# PCA after all filtration to test the association between PC and batch ID variable 
as.all.samplesgroup$association[, 1:10]
 write.csv(as.all.samplesgroup$association, "R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/associationbetweenPCandcovariates.csv")
```
## Removing the PC 1 and PC2 to account for the batch effect 

``` {r echo= TRUE, message = FALSE}
new.meth.filterdgroup = removeComp(meth.filteredgroup,comp=c (1, 2))
as.filtered <-assocComp(mBase=new.meth.filterdgroup, metadata)
write.csv(as.all.samplesgroup$association, "R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/associationbetweenPCandcovariates-afterbatcheffectcorrection.csv")
methylationall <- percMethylation (new.meth.filterdgroup)
rownames (methylationall) <- paste(getData(new.meth.filterdgroup)[, 1], ".", getData(new.meth.filterdgroup)[,2], sep ="")
write.csv (methylationall, "R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/BDvsCO/total-methylation-allsamples.csv")
```

## DMA at the level of one CpG sites

``` {r echo= TRUE, message = FALSE}
covariates0 <-  metadata[, c("Race", "Sex", "Age", "BMI",  "Smoking.History")]
my.diffMeth0 <-calculateDiffMeth(new.meth.filterdgroup,
                               covariates=covariates0,
                               overdispersion="MN",test="Chisq",mc.cores=1)
write.csv(my.diffMeth0,"R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/BDvsCO/DMA-BDvscon.csv" )
```

## DMA at the level of once CpG sites-Hypermetylated CpG sites

Hypermethylated CpG sites

``` {r echo= TRUE, message = FALSE }
myDiff25p0.hyper=getMethylDiff(my.diffMeth0,difference=0,
                              qvalue=0.05,type="hyper")
myDiff25p0.hyper
```

## DMA at the level of one CpG sites-Hypomethylated CpG sites

Hypomethylated CpG sites

``` {r echo= TRUE, message = FALSE }
myDiff25p0.hypo = getMethylDiff(my.diffMeth0,difference=0,
                             qvalue=0.05,type="hypo")
myDiff25p0.hypo
```

## DMA at the level of one CpG site-Signficant differentially methyalted CpG 
``` {r echo= TRUE, message = FALSE }
myDiff25p0 =getMethylDiff(my.diffMeth0,difference=0,qvalue=0.05)
myDiff25p0
```
## visualize the distribution of hypo/hyper-methylated bases per chromosome

``` {r echo= TRUE, message = FALSE}
bedgraph(myDiff25p0, col.name = "meth.diff", file.name = "diff_cpg_25p.bed")
```
## DMA at the level of region

``` {r echo= TRUE, message = FALSE }
tiles = tileMethylCounts(new.meth.filterdgroup,win.size=300,step.size=300,cov.bases = 2)
head (tiles)
```

## DMA at the level of region

``` {r echo= TRUE, message = FALSE }
mydiffdmr0 <- calculateDiffMeth(tiles,
                               covariates=covariates0,
                               overdispersion="MN",test="Chisq",mc.cores=1)
head (mydiffdmr0)
bedgraph(myDiff25p0, col.name = "meth.diff", file.name = "diff_cpg_25p.bed")
write.csv(mydiffdmr0,"R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/BDvsCO/DMA-BDvscon-DMR.csv" )
```
