
library("MatrixEQTL")
library (hrbrthemes )
library (ggplot2)
library(lmtest)
base.dir = "R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/MatrixeQTL/Matrix eQTL"
useModel = modelLINEAR
SNP_file_name = paste(base.dir, "/Genoptypes.csv", sep="");
genotypes <- read.csv (SNP_file_name)
rownames (genotypes) <- genotypes[,1]
genotypes <- genotypes [,-1]
#methylation file
methylation_file_name = paste(base.dir, "/ME.csv", sep="");
methylation <- read.csv (methylation_file_name, header = TRUE)
rownames (methylation)<- methylation[,1]
methylation <- methylation [,-1]
#A separate file may be provided with extra covariates. In case of no covariates set the variable covariates_file_name to character().

covariates_file_name = paste(base.dir, "/final-covariates.csv", sep="");
covariates <- read.csv(covariates_file_name)
rownames (covariates) <- covariates [,1]
covariates <- covariates[, -1]
# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

#SNPs location file name
SNPsloc_file_name = paste(base.dir, "/snploc.txt", sep="")

#Gene location file name
Geneloc_file_name = paste(base.dir, "/geneloc.txt", sep="")
# Error covariance matrix
errorCovariance = numeric()
cvrt = SlicedData$new()
# Distance for local gene-SNP pairs
cisDist = 1e6;
## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(methylation_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "./."; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
pvOutputThreshold_cis =1
pvOutputThreshold_tra = 1
snpspos = read.table(SNPsloc_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(Geneloc_file_name, header = TRUE, stringsAsFactors = FALSE);
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
## Results:

message('Analysis done in: ', me$time.in.sec, ' seconds');
message('Detected eQTLs:');
show(me$trans$eqtls);

## Plot the histogram of all p-values
plot(me)
#look at the number and type of eQTL
me$all$neqtls
eqtlcis <- me$cis$eqtls
eqtlcissig <- eqtlcis[ which (eqtlcis$pvalue < 0.05),]
eqtltrans <- me$trans$eqtls
eqtltranssig <- eqtltrans[ which (eqtltrans$FDR < 0.05),]

write.csv(eqtlcis,'eqtl-modelLINEAR_ME_SNPs_cis.csv')
write.csv(eqtltrans,"eqtl-modelLINEAR_ME_SNPs_trans.csv" )
## Plot the Q-Q plot of local and distant p-values


for (i in 1:nrow(eqtltranssig)){
  df.figure <-  as.data.frame (methylation [paste (eqtltranssig$gene[i]), ]) 
  colnames (df.figure) <- colnames(gene)# getting the columns CpG site in total we have 90
  df.figure <- data.frame(t(df.figure))
  df.figure$samples.id <- rownames (df.figure) # adding the sample ID as one column
  df.figure$SNP <- as.factor(t(genotypes[paste (eqtltranssig$snps[i]), ]))
  df.figure$Group <- t(covariates["final-group", ])
  df.figure$Group <- as.factor (ifelse (df.figure$Group == 2, "BDSB", ifelse(df.figure$Group == 1, "BDNSB", "Control"))) # setting the level of suicide history 0,1 
  df.figure$SNP[df.figure$SNP == "exm768107"] = "rs75981117"
  df.figure$SNP[df.figure$SNP == "exm-rs300774"] = "rs300774"
  file_name <- paste ("R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/MatrixeQTL/Matrix eQTL/figurestrans/", paste (eqtltranssig$gene[i], eqtltranssig$snp [i]),".tiff", sep ="")
  tiff(file_name, units="in", width=8, height=5, res=300)
  
  print (ggplot(data = df.figure,  aes(x= SNP, y= df.figure[,1], fill=Group))+geom_boxplot() +
           scale_fill_brewer(palette="RdBu")+
           theme_ipsum() +
           theme(
             legend.position="bottom",
             plot.title = element_text(size=14)) + ggtitle(paste ("Methylation level of", paste (eqtltranssig$gene[i])," vs.", paste (eqtltranssig$snp [i]),   "Genotypes" )) +
           xlab("SNP Genotypes") + ylab("Methylation level (log (methylation %)")+
           geom_point() )
  
  dev.off()}
for (i in 1:nrow(eqtlcissig)){
  df.figure <-  as.data.frame (methylation [paste (eqtlcissig$gene[i]), ]) 
  colnames (df.figure) <- colnames(gene)# getting the columns CpG site in total we have 90
  df.figure <- data.frame(t(df.figure))
  df.figure$samples.id <- rownames (df.figure) # adding the sample ID as one column
  df.figure$SNP <- as.factor(t(genotypes[paste (eqtlcissig$snps[i]), ]))
  df.figure$Group <- t(covariates["final-group", ])
  df.figure$Group <- as.factor (ifelse (df.figure$Group == 2, "BDSB", ifelse(df.figure$Group == 1, "BDNSB", "Control"))) # setting the level of suicide history 0,1 

  file_name <- paste ("R:/Aysheh/Gaine Lab/BD-SB/Analysis/Final/MatrixeQTL/Matrix eQTL/figurescis/", paste (eqtlcissig$gene[i], eqtlcissig$snp [i]),".tiff", sep ="")
  tiff(file_name, units="in", width=8, height=5, res=300)
  
  print (ggplot(data = df.figure,  aes(x= SNP, y= df.figure[,1], fill=SNP))+geom_boxplot() +
           scale_fill_brewer(palette="RdBu")+
           theme_ipsum() +
           theme(
             legend.position="bottom",
             plot.title = element_text(size=14)) + ggtitle(paste ("Methylation level of", paste (eqtlcissig$gene[i])," vs.", paste (eqtlcissig$snp [i]),   "Genotypes" )) +
           xlab("SNP Genotypes") + ylab("Methylation level (log (methylation %)")+
           geom_point() )
  
  dev.off()}