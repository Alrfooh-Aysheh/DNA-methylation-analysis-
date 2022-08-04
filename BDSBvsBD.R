library(DSS)
require(bsseq)
library(edgeR)
library (tidyr)
library (ggplot2)
library(viridis)
library(ggplot2)
library(ggplot2)
library(dplyr)
library (hrbrthemes)
library(ggplot2)
library(gridExtra)
library(extrafont)
library (methylKit)
library(ggplot2)
library(gtable)
library(gridExtra)
library(grid)
library(ggsignif)
library(MASS)
library(car)
library(dummies)
library(mlr)
library(emmeans)
library(tidyverse)
library(broom)
library(ggfortify)
library(RColorBrewer)
library(outliers)
library(corrplot)
library(stargazer)
library(foreign)  
library(car)      
library(gvlma)   
library(effects) 
library(sjPlot)  
library(MASS) 
library(relaimpo) 
library(ggplot2)
library(gplots)
library ("ben-laufer/DMRichR")
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
BiocManager::install("ben-laufer/DMRichR")
####Reading the files#####

#reading the file that has the list of covriates
Covariates <- read.csv("R:\\Aysheh\\Gaine Lab\\BD-SB\\Analysis\\Final\\Covariates.csv")

#reading the file that has the list of the methylation data 
meth <- read.csv( "R:\\Aysheh\\Gaine Lab\\BD-SB\\Analysis\\Final\\methdata.csv")
meth <- meth[,-1]
#reading the list of the sample IDs

sample_id <- read.csv ( "R:\\Aysheh\\Gaine Lab\\BD-SB\\Analysis\\Final\\sample_id.csv")
sample_id <- as.data.frame (sample_id[,2])
colnames (sample_id)<- "sample_id"

#####Converting the datframe meth to list ####
DSSlist= list()
x <- 0
for (i in seq (from =5, to = length (colnames(meth)), by=4))
{
  column5= i
  Column6 = i+1
  column7= i+2
  column8 = i+3
  df= data.frame (meth [c (1:4,column5,Column6,column7, column8 )])
  colnames (df)<- c("chrBase", "chr", "pos","S",  "N", "X","T", "P" )
  dat <- data.frame(df$chr, df$pos, df$N, df$X)
  colnames (dat)<- c("chr", "pos",  "N", "X" )
  x <- x+1
  DSSlist[[x]] <- dat
  print (x)
}
names (DSSlist) <- c(sample_id$sample_id)
#####merging two dataframes by ID####
match(sample_id$sample_id, Covariates$sample_id)
total <- merge (sample_id, Covariates, by ="sample_id")
#####extracting the covraites for the first group#######
sample = as.vector(total$sample_id[total$group == "2"]) # we need to the list of the sample for the group 2 to exclude from DSSlist 
total <- subset (total, group == 1 )#)# including only one group (Bipolar)
dim (total) 
colnames (total)
variable  <- data.frame (as.factor(total$history_of_suicide_attempt), as.factor(total$Gender), total$Age, as.factor (total$batch_id), as.factor( (total$Ethnicity)))# reading the variable in the cvs file for batch effect correction 
dim (variable)
colnames(variable) <- c("history_of_suicide_attempt", "gender","age", "batch_id", "Ethnicity")

#######Differentiall methylation analysis##########
## make BSseq objects

DSSlist[sample] <- NULL
length (DSSlist)
BSobj <- makeBSseqData( DSSlist,
                        as.vector (unlist (total$sample_id)) )
BSobj
sampleNames(BSobj)
#  x = model.matrix (~history_of_suicide_attempt+ gender+age+ batch_id, covariates)
# contrast= matrix (c(0,0,0,0,1), ncol=1)
DMLfit= DMLfit.multiFactor(BSobj, variable,formula= ~age+ batch_id+ Ethnicity +history_of_suicide_attempt, smoothing=TRUE)
## hypothesis testing
DMLtest.suicide = DMLtest.multiFactor(DMLfit, coef="history_of_suicide_attempt1")
# DMLtest.suicide.term = DMLtest.multiFactor(DMLfit, term ="batch_id")
write.csv(DMLtest.suicide, "DSS_DML_Results-total.csv")
dmrs = callDMR(DMLtest.suicide, p.threshold = 0.05)
write.csv (dmrs, "DSS_DMR_Results.csv")


## look at distributions of test statistics and p-values
# eg. for tiff()
# eg. for tiff()
par(mar=c(1,1,1,1))
tiff(filename =  "qq.tiff",
     res = 300,                                                 # the margin error.
     width = 5, height = 4, units = 'in',                       # fixed
     compression = c("lzw") )


 showOneDMR(dmrs[1,], BSobj)


dev.off()

