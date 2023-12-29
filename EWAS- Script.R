#This script demonstrates one of the EWAS pipelines I used in combination with 
#the RODAM study. This will not include the sex specific EWAS. 50% of this pipeline
#was provided for me by the department's bioinformaticians the other 50% was self-written.

#################################################################################################
## THE METHYLAID AND VISUALISATION CODE

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

library("MethylAid")
library("IlluminaHumanMethylation450kmanifest")
library(minfi)
library(minfiData)
library(limma)
library(affy)
library(IlluminaHumanMethylationEPICmanifest)

options(stringsAsFactors = F)

idats_RODAM = "~/lkg/Rodam/idats/"

idats <- read.metharray.sheet("~/lkg/Rodam/idats/")
data <- summarize(idats)
visualize(data,thresholds=list(epic = list(MU = 10, OP = 12, BS = 11.75, HC = 12.75, DP = 0.95)) )

#################################################################################################
## SEX PREDICTION CODE

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

##For using the test sample 
idats_RODAM = "~/lkg/Rodam/idats/"
idatspath="~/lkg/Rodam/idats_test"
targets = read.csv("Rodam_Data_Simplified.csv", sep=",",header=T)
targets = read.csv("Testtarget.csv", sep=",", header=T)
prob = read.csv("~/lkg/Rodam/STUDENT_CHARGE_TRAINING_DATA/Charge/subsetSexProbes.txt",stringsAsFactors  = F, col.names = 1)[,1]

BasePath = idats_RODAM
BasePath=idatspath
RGset=read.metharray.exp(BasePath,targets)
RGset2=preprocessRaw(RGset)


##NON-VISUAL METHYLAID FUNCTION
methylAid_nonvisual <- function(RGset,RGset2, targets, MU_threshold = 10, NP_threshold = 12, BS_threshold = 11.75, HC_threshold = 12.75, DP_threshold = 0.95){
  
  resultsMethyl=data.frame(matrix(ncol = 5, nrow = length(targets$Basename)))
  rownames(resultsMethyl)=targets$Basename
  colnames(resultsMethyl)=c("MU","OP","BS","HC","DP")
  resultsMethyl[is.na(resultsMethyl)]<- "pass"
  
  #Calculate detection p-value and frequency of probe passing threshold
  DP <- detectionP(RGset)
  DPfreq <- colSums(DP < 0.01, na.rm=TRUE)/nrow(DP)
  
  #Medians of the overall methylated and unmethylated intensities
  MU.full <- matrix(0.0, nrow = 2, ncol = ncol(RGset))
  M.full <- getMeth(RGset2)
  U.full <- getUnmeth(RGset2)
  MU.full[1,] <- colMedians(M.full, na.rm = TRUE)
  MU.full[2,] <- colMedians(U.full, na.rm = TRUE)
  colnames(MU.full) <- colnames(RGset)
  rownames(MU.full) <- c("Methylated", "Unmethylated")
  
  #Find the red and green intensities of the control probes
  data(hm450.controls, package="FDb.InfiniumMethylation.hg19", envir=environment())
  Red.full <- getRed(RGset)
  Green.full <- getGreen(RGset)
  id <- intersect(hm450.controls$Address, rownames(Red.full))
  Red.controls <- Red.full[rownames(Red.full) %in% id,]
  Green.controls <- Green.full[rownames(Green.full) %in% id,]
  hm450.controls <- hm450.controls[hm450.controls$Address %in% id,]
  hm450.controls <- hm450.controls[order(hm450.controls$Address), ]
  
  Red.controls.log <- log2(Red.controls)
  Green.controls.log <- log2(Green.controls)
  
  hm450.controls <- hm450.controls[!(hm450.controls$Type %in% c("NORM_A", "NORM_G", "NORM_C", "NORM_T")), ]
  control.data <- data.frame(Address=rep(rownames(Red.controls.log), ncol(Red.controls.log)), Samples=rep(colnames(Red.controls.log), each=nrow(Red.controls.log)), IntRed=as.vector(Red.controls.log),IntGrn=as.vector(Green.controls.log))
  hm450.control.data <- merge(hm450.controls, control.data)
  
  #Methylated/Unmethylated
  MU.full.log2 <- t(as.data.frame(log2(MU.full)))
  MU.outliers <- names(which(MU.full.log2[,1] <= MU_threshold))
  
  #Non-Polymorphic
  NP.control <- hm450.control.data[grepl("^NON-POLYMORPHIC$", hm450.control.data$Type),]
  NP.green <- NP.control[NP.control$Name %in% c("NP (C)", "NP (G)"), c(1:5,7)]
  NP.avg.green <- tapply(NP.green$IntGrn, NP.green$Samples, mean)
  NP.outliers <- names(which(NP.avg.green <= NP_threshold))
  
  #Bisulfite Conversion I
  BS.control <- hm450.control.data[grepl("^BISULFITE CONVERSION I$", hm450.control.data$Type),]
  BS.green <- BS.control[grepl("C1|C2|C3", BS.control$Name), c(1:5,7)]
  BS.avg.green <- tapply(BS.green$IntGrn, BS.green$Samples, mean)
  BS.outliers <- names(which(BS.avg.green <= BS_threshold))
  
  #Hybridization
  HC.control <- hm450.control.data[grepl("^HYBRIDIZATION$", hm450.control.data$Type),]
  HC.control <- HC.control[order(HC.control$Samples),]
  HC.green <- 0.5*(HC.control$IntGrn[grepl("High", HC.control$Name)] + HC.control$IntGrn[grepl("Low", HC.control$Name)])
  names(HC.green) <- HC.control$Samples[grepl("High", HC.control$Name)]
  HC.outliers <- names(which(HC.green <= HC_threshold))
  
  #Detection P-values
  DP.outliers <- names(which(DPfreq <= DP_threshold))
  
  for( i in MU.outliers){
    resultsMethyl[i,"MU"]="FAIL"
    
  }
  for( i in NP.outliers){
    resultsMethyl[i,"OP"]="FAIL"
    
  }
  for( i in HC.outliers){
    resultsMethyl[i,"HC"]="FAIL"
    
  }
  for( i in BS.outliers){
    resultsMethyl[i,"BS"]="FAIL"
    
  }
  for( i in DP.outliers){
    resultsMethyl[i,"DP"]="FAIL"
    
  }
  
  #Return the corrected target sheet
  return( resultsMethyl= resultsMethyl)
}


metresults=methylAid_nonvisual(RGset,RGset2,targets)
difprob=RGset2 [which(featureNames(RGset2) %in% prob), ]
beta_x=getBeta(difprob)
Averagex=colMeans(na.omit(beta_x))

metresults[,"PredictedSex"]=NA
count=1
for(i in Averagex){
  if (i>0.3){
    metresults[count,"PredictedSex"]="F"
    count=count+1
  }
  else if(i<0.3){
    metresults[count,"PredictedSex"]="M"
    count=count+1
  }}

datares=cbind(metresults,targets)
write.csv(datares, "Sex Prediction File.csv")

#################################################################################################
## CODE FOR BETA AND M VALUES

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

idat_phenofile = read.csv("PhenoRodam15032106Analyse.csv", sep=",",header=T)  #Eva: replaced "PhenoRodam15032106Analyse1.csv" by "PhenoRodam15032106Analyse.csv"
targets = read.csv("Rodam_Data_Simplified.csv", sep=",",header=T)

raw_setRODAM = read.metharray.exp(idats_RODAM, targets, force=TRUE) 

raw_betaRODAM = getBeta(raw_setRODAM)
normal_beta_setRODAM = preprocessFunnorm(raw_setRODAM, nPCs=2, sex = c("M","F"), bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)
pdRODAM <- pData(normal_beta_setRODAM)
anno_RODAM <- getAnnotation(normal_beta_setRODAM)
probe.features=as.matrix(anno_RODAM)
indices <- which(anno_RODAM$chr == "chrX"| anno_RODAM$chr == "chrY")
GMset_no_sex <- normal_beta_setRODAM[-indices,]
GMset_culled1 <-dropLociWithSnps(GMset_no_sex, snps = c("CpG", "SBE"), maf = 0.01, snpAnno = NULL)

probe_removal <- read.csv("crosshybEnSNPSRODAM.txt", stringsAsFactors  = F, col.names = 1)[,1]

GMset_culled <- GMset_culled1 [which(!featureNames(GMset_culled1) %in% probe_removal), ]

m_val = getM(GMset_culled)
beta = getBeta(GMset_culled)
x=which(m_val == -Inf, arr.ind = T)
xx=sort(rownames(x))
g=m_val
zz1=g[!rownames(g) %in% xx, ]
mvalue=zz1
betavalue = beta[!rownames(beta) %in% xx, ]

write.table(betavalue,"Test RODAM Beta.txt", sep="\t")
write.table(raw_betaRODAM, "Test Raw Beta Values.txt",sep="\t")
write.table(mvalue,"Test RODAM M value.txt",sep="\t")

#################################################################################################
## CODE FOR THE CELL TYPE DISTRIBUTION WITHIN THE RODAM DATA

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library("FlowSorted.Blood.450k")
library("FlowSorted.Blood.EPIC")
options(stringsAsFactors = F)

raw_setRODAM = read.metharray.exp(idats_RODAM, targets,force=TRUE)
pdRODAM <- pData(normal_beta_setRODAM)

cell_distribution=estimateCellCounts(raw_setRODAM, compositeCellType = "Blood",
     cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
     returnAll = T, meanPlot = F, verbose = F)
cell_counts=cell_distribution$counts

rownames(cell_counts)=pdRODAM$Basename
write.table(cell_counts,"Cell Type Table.txt",sep="\t") 
cell_counts = read.table("Cell Type Table.txt", sep="\t", header = T)
targets = merge(targets, cell_counts, by="Basename")
targets = subset(targets, select = -c(Basename, X))

#################################################################################################
## CODE FOR DETECTING DMP IN THE POPULATION
 
setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(minfiData)
library(limma)
library(affy)
library(expss)

options(stringsAsFactors = F) 

mvalue = read.table("RODAM M value.txt", sep="\t", header=T)
betavalue = read.table("RODAM Beta.txt", sep="\t", header=T)

densityPlot(raw_setRODAM, sampGroups = pd$group, main = "Beta", xlab = "Beta")
par(oma=c(2,10,1,1))
densityBeanPlot(raw_set, sampGroups = pd$group, sampNames = pd$Sample_ID)
controlStripPlot(raw_setRODAM, controls ="BISULFITE CONVERSION II", sampNames = pd$Sample_ID)

#Normalised density plots
densityPlot(beta, sampGroups = pd$group, main = "Beta", xlab = "Beta")
par(oma=c(2,10,1,1))
densityBeanPlot(beta, sampGroups = pd$group, sampNames = pd$Sample_ID)


Location = factor(targets$Location)
sex = factor(targets$Sex)
age = targets$Age

leg = targets$Leg_Length
sitting = targets$Sitting_Calculation
trunk = targets$Trunk_Length
whr = targets$WHR

hyperten = factor(targets$Hyperten_Meds.2)
diabetes = factor(targets$DiabetesYN)
obesity = targets$BMI

cd8t = targets$CD8T
cd4t = targets$CD4T
mono= targets$Mono
nk = targets$NK
bcell= targets$Bcell
gran= targets$Gran

batch = targets$Batch_BSR.x
batchA = targets$BATCH_arry.x
pos = targets$Sentrix_Position.x


#DMP model matrices 

leg_modela= model.matrix(~leg+cd8t+cd4t+mono+nk+bcell+gran+batch+batchA+pos+age
    +sex)
leg_modelb= model.matrix(~leg+cd8t+cd4t+mono+nk+bcell+gran+batch+batchA+pos+age
    +sex+Location)
sit_modela= model.matrix(~sitting+cd8t+cd4t+mono+nk+bcell+gran+batch+batchA+pos+age
    +sex)
sit_modelb= model.matrix(~sitting+cd8t+cd4t+mono+nk+bcell+gran+batch+batchA+pos+age
    +sex+Location)
trunk_modela= model.matrix(~trunk+cd8t+cd4t+mono+nk+bcell+gran+batch+batchA+pos+age
     +sex)
trunk_modelb= model.matrix(~trunk+cd8t+cd4t+mono+nk+bcell+gran+batch+batchA+pos+age
     +sex+Location)

#Produce fit models for M values

fit_leg2 = lmFit(mvalue, leg_modela)
fit_leg3 = lmFit(mvalue, leg_modelb)
fit_sit2 = lmFit(mvalue, sit_modela)
fit_sit3 = lmFit(mvalue, sit_modelb)
fit_trunk2 = lmFit(mvalue, trunk_modela)
fit_trunk3 = lmFit(mvalue, trunk_modelb)

fit_leg2e = eBayes(fit_leg2)
fit_leg3e = eBayes(fit_leg3)
fit_sit2e = eBayes(fit_sit2)
fit_sit3e = eBayes(fit_sit3)
fit_trunk2e = eBayes(fit_trunk2)
fit_trunk3e = eBayes(fit_trunk3)

#Produce FDR models top tables

leg2_topTableFDR= topTable(fit_leg2e, number=nrow(mvalue), adjust="fdr")
leg3_topTableFDR= topTable(fit_leg3e, number=nrow(mvalue), adjust="fdr")
sit2_topTableFDR= topTable(fit_sit2e, number=nrow(mvalue), adjust="fdr")
sit3_topTableFDR= topTable(fit_sit3e, number=nrow(mvalue), adjust="fdr")
trunk2_topTableFDR= topTable(fit_trunk2e, number=nrow(mvalue), adjust="fdr")
trunk3_topTableFDR= topTable(fit_trunk3e, number=nrow(mvalue), adjust="fdr")

leg2_topTable = topTable(fit_leg2e, number=nrow(mvalue), adjust="none")
leg3_topTable = topTable(fit_leg3e, number=nrow(mvalue), adjust="none")
sit2_topTable = topTable(fit_sit2e, number=nrow(mvalue), adjust="none")
sit3_topTable = topTable(fit_sit3e, number=nrow(mvalue), adjust="none")
trunk2_topTable = topTable(fit_trunk2e, number=nrow(mvalue), adjust="none")
trunk3_topTable = topTable(fit_trunk3e, number=nrow(mvalue), adjust="none")

#Calculate P Value inflation
t_stats <- fit2e$t
t_stats = bacon(t_stats)
t_stats <- pval(t_stats)

leg2_tstats <- fit_leg2e$t
leg2_tstats = bacon(leg2_tstats)
leg2_tstats <- pval(leg2_tstats)

leg3_tstats <- fit_leg3e$t
leg3_tstats = bacon(leg3_tstats)
leg3_tstats <- pval(leg3_tstats)

sit2_tstats <- fit_sit2e$t
sit2_tstats = bacon(sit2_tstats)
sit2_tstats <- pval(sit2_tstats)

sit3_tstats <- fit_sit3e$t
sit3_tstats = bacon(sit3_tstats)
sit3_tstats <- pval(sit3_tstats)

trunk2_tstats <- fit_trunk2e$t
trunk2_tstats = bacon(trunk2_tstats)
trunk2_tstats <- pval(trunk2_tstats)

trunk3_tstats <- fit_trunk3e$t
trunk3_tstats = bacon(trunk3_tstats)
trunk3_tstats <- pval(trunk3_tstats)

#Make your base data frame for your DMP
detach("package:bacon", unload=TRUE)

x_leg2 =data.frame(probe.features[match(row.names(leg2_topTableFDR), row.names(probe.features)),])
x_leg3 =data.frame(probe.features[match(row.names(leg3_topTableFDR), row.names(probe.features)),])

x_sit2 =data.frame(probe.features[match(row.names(sit2_topTableFDR), row.names(probe.features)),])
x_sit3 =data.frame(probe.features[match(row.names(sit3_topTableFDR), row.names(probe.features)),])

x_trunk2 =data.frame(probe.features[match(row.names(trunk2_topTableFDR), row.names(probe.features)),])
x_trunk3 =data.frame(probe.features[match(row.names(trunk3_topTableFDR), row.names(probe.features)),])



z_leg2 = data.frame(leg2_topTable[match(row.names(leg2_topTableFDR),row.names(leg2_topTable)),])
z_leg3 = data.frame(leg3_topTable[match(row.names(leg3_topTableFDR),row.names(leg3_topTable)),])

z_sit2 = data.frame(sit2_topTable[match(row.names(sit2_topTableFDR),row.names(sit2_topTable)),])
z_sit3 = data.frame(sit3_topTable[match(row.names(sit3_topTableFDR),row.names(sit3_topTable)),])

z_trunk2 = data.frame(trunk2_topTable[match(row.names(trunk2_topTableFDR),row.names(trunk2_topTable)),])
z_trunk3 = data.frame(trunk3_topTable[match(row.names(trunk3_topTableFDR),row.names(trunk3_topTable)),])

###THIS IS WHERE THE INFLATED P VALUE ERROR OCCURS! 1. THE P VALUE AND THE ADJUSTED 
##P VALUE ARE OVER INFLATED 2. THERE IS NO VARIABLE PRODUCED CALLED "LOGFC" WHICH
#WE NEED FOR THE NEXT SECTION
#Bind the data frames with top FDR tables and then rename columns


y_leg2 = z_leg2$logFC
leg2_topTableFDRannot = cbind(leg2_topTableFDR, y_leg2)
colnames(leg2_topTableFDRannot)[7]="Mvalue_Delta"
colnames(leg2_topTableFDRannot)[1]="Mvalue_logFC"

y_leg3 = z_leg3$logFC
leg3_topTableFDRannot = cbind(leg3_topTableFDR, y_leg3)
colnames(leg3_topTableFDRannot)[7]="Mvalue_Delta"
colnames(leg3_topTableFDRannot)[1]="Mvalue_logFC"

y_sit2 = z_sit2$logFC
sit2_topTableFDRannot = cbind(sit2_topTableFDR, y_sit2)
colnames(sit2_topTableFDRannot)[7]="Mvalue_Delta"
colnames(sit2_topTableFDRannot)[1]="Mvalue_logFC"

y_sit3 = z_sit3$logFC
sit3_topTableFDRannot = cbind(sit3_topTableFDR, y_sit3)
colnames(sit3_topTableFDRannot)[7]="Mvalue_Delta"
colnames(sit3_topTableFDRannot)[1]="Mvalue_logFC"

y_trunk2 = z_trunk2$logFC
trunk2_topTableFDRannot = cbind(trunk2_topTableFDR, y_trunk2)
colnames(trunk2_topTableFDRannot)[7]="Mvalue_Delta"
colnames(trunk2_topTableFDRannot)[1]="Mvalue_logFC"

y_trunk3 = z_trunk3$logFC
trunk2_topTableFDRannot = cbind(trunk3_topTableFDR, y_trunk3)
colnames(trunk3_topTableFDRannot)[7]="Mvalue_Delta"
colnames(trunk3_topTableFDRannot)[1]="Mvalue_logFC"


################## REPEAT EVERYTHING WITH BETA VALUES ##########################
fit_leg2BETA = lmFit(betavalue, leg_modela)
fit_leg3BETA = lmFit(betavalue , leg_modelb)
fit_sit2BETA = lmFit(betavalue, sit_modela)
fit_sit3BETA = lmFit(betavalue, sit_modelb)
fit_trunk2BETA = lmFit(betavalue, trunk_modela)
fit_trunk3BETA = lmFit(betavalue, trunk_modelb)

fit_leg2eBETA = eBayes(fit_leg2BETA)
fit_leg3eBETA = eBayes(fit_leg3BETA)
fit_sit2eBETA = eBayes(fit_sit2BETA)
fit_sit3eBETA = eBayes(fit_sit3BETA)
fit_trunk2eBETA = eBayes(fit_trunk2BETA)
fit_trunk3eBETA = eBayes(fit_trunk3BETA)


leg2_topTableFDRBETA = topTable(fit_leg2eBETA, number=nrow(betavalue), adjust="fdr")
leg3_topTableFDRBETA = topTable(fit_leg3eBETA, number=nrow(betavalue), adjust="fdr")
sit2_topTableFDRBETA = topTable(fit_sit2eBETA, number=nrow(betavalue), adjust="fdr")
sit3_topTableFDRBETA = topTable(fit_sit3eBETA, number=nrow(betavalue), adjust="fdr")
trunk2_topTableFDRBETA = topTable(fit_trunk2eBETA, number=nrow(betavalue), adjust="fdr")
trunk3_topTableFDRBETA = topTable(fit_trunk3eBETA, number=nrow(betavalue), adjust="fdr")

leg2_topTableBETA = topTable(fit_leg2eBETA, number=nrow(betavalue), adjust="none")
leg3_topTableBETA = topTable(fit_leg3eBETA, number=nrow(betavalue), adjust="none")
sit2_topTableBETA = topTable(fit_sit2eBETA, number=nrow(betavalue), adjust="none")
sit3_topTableBETA = topTable(fit_sit3eBETA, number=nrow(betavalue), adjust="none")
trunk2_topTableBETA = topTable(fit_trunk2eBETA, number=nrow(betavalue), adjust="none")
trunk3_topTableBETA = topTable(fit_trunk3eBETA, number=nrow(betavalue), adjust="none")


x_leg2BETA = data.frame(probe.features[match(row.names(leg2_topTableFDRBETA), row.names(probe.features)),])
x_leg3BETA = data.frame(probe.features[match(row.names(leg3_topTableFDRBETA), row.names(probe.features)),])

x_sit2BETA = data.frame(probe.features[match(row.names(sit2_topTableFDRBETA), row.names(probe.features)),])
x_sit3BETA = data.frame(probe.features[match(row.names(sit3_topTableFDRBETA), row.names(probe.features)),])

x_trunk2BETA = data.frame(probe.features[match(row.names(trunk2_topTableFDRBETA), row.names(probe.features)),])
x_trunk3BETA = data.frame(probe.features[match(row.names(trunk3_topTableFDRBETA), row.names(probe.features)),])


z_leg2BETA = data.frame(leg2_topTableBETA[match(row.names(leg2_topTableFDRBETA),row.names(leg2_topTableBETA)),])
z_leg3BETA = data.frame(leg3_topTableBETA[match(row.names(leg3_topTableFDRBETA),row.names(leg3_topTableBETA)),])

z_sit2BETA = data.frame(sit2_topTableBETA[match(row.names(sit2_topTableFDRBETA),row.names(sit2_topTableBETA)),])
z_sit3BETA = data.frame(sit3_topTableBETA[match(row.names(sit3_topTableFDRBETA),row.names(sit3_topTableBETA)),])

z_trunk2BETA = data.frame(trunk2_topTableBETA[match(row.names(trunk2_topTableFDRBETA),row.names(trunk2_topTableBETA)),])
z_trunk3BETA = data.frame(trunk3_topTableBETA[match(row.names(trunk3_topTableFDRBETA),row.names(trunk3_topTableBETA)),])


y_leg2BETA = z_leg2BETA$logFC
leg2_topTableFDRannotBETA = cbind(leg2_topTableFDRBETA, y_leg2BETA, x_leg2BETA)
colnames(leg2_topTableFDRannotBETA)[7]="Beta_Delta"
colnames(leg2_topTableFDRannotBETA)[1]="Beta_logFC"

y_leg3BETA = z_leg3BETA$logFC
leg3_topTableFDRannotBETA = cbind(leg3_topTableFDRBETA, y_leg3BETA, x_leg3BETA)
colnames(leg3_topTableFDRannotBETA)[7]="Beta_Delta"
colnames(leg3_topTableFDRannotBETA)[1]="Beta_logFC"

y_sit2BETA = z_sit2BETA$logFC
sit2_topTableFDRannotBETA = cbind(sit2_topTableFDRBETA, y_sit2BETA, x_sit2BETA)
colnames(sit2_topTableFDRannotBETA)[7]="Beta_Delta"
colnames(sit2_topTableFDRannotBETA)[1]="Beta_logFC"

y_sit3BETA = z_sit3BETA$logFC
sit3_topTableFDRannotBETA = cbind(sit3_topTableFDRBETA, y_sit3BETA, x_sit3BETA)
colnames(sit3_topTableFDRannotBETA)[7]="Beta_Delta"
colnames(sit3_topTableFDRannotBETA)[1]="Beta_logFC"

y_trunk2BETA = z_trunk2BETA$logFC
trunk2_topTableFDRannotBETA = cbind(trunk2_topTableFDRBETA, y_trunk2BETA, x_trunk2BETA)
colnames(trunk2_topTableFDRannotBETA)[7]="Beta_Delta"
colnames(trunk2_topTableFDRannotBETA)[1]="Mvalue_logFC"

y_trunk3BETA = z_trunk3BETA$logFC
trunk2_topTableFDRannotBETA = cbind(trunk3_topTableFDRBETA, y_trunk3BETA, x_trunk3BETA)
colnames(trunk3_topTableFDRannotBETA)[7]="Beta_Delta"
colnames(trunk3_topTableFDRannotBETA)[1]="Beta_logFC"



#Finish your DMP with the "z" dataframe, correct for pvalues, and then combine all 
##dataframe

leg2_corrected_p = z_leg2 = data.frame(leg2_topTableFDRannotBETA[match(row.names
     (leg2_topTableFDRannot),row.names(leg2_topTableFDRannotBETA)),])
leg2_bcp.value <- data.frame(leg2_tstats[match(row.names(leg2_corrected_p), row.names(leg2_tstats)),])
leg2_bcadjust <- as.matrix(p.adjust(as.matrix(leg2_tstats), "fdr", n=length(leg2_tstats)))
leg2_tstats <- leg2_tstats[, 1]
row.names(leg2_bcadjust) = row.names(leg2_tstats)
colnames(leg2_bcp.value)[1]="BaconPvalue_Mvalues"
colnames(leg2_bcadjust)[1]="Baconadjust.pval_Mvalues"
leg2_TOPTABLE=cbind(leg2_topTableFDRannot, leg2_bcp.value, leg2_bcadjust, leg2_corrected_p)

leg3_corrected_p = z_leg3 = data.frame(leg3_topTableFDRannotBETA[match(row.names
     (leg3_topTableFDRannot),row.names(leg3_topTableFDRannotBETA)),])
leg3_bcp.value <- data.frame(leg3_tstats[match(row.names(leg3_corrected_p), row.names(leg3_tstats)),])
leg3_bcadjust <- as.matrix(p.adjust(as.matrix(leg3_tstats), "fdr", n=length(leg3_tstats)))
leg3_tstats <- leg3_tstats[, 1]
row.names(leg3_bcadjust) = row.names(leg3_tstats)
colnames(leg3_bcp.value)[1]="BaconPvalue_Mvalues"
colnames(leg3_bcadjust)[1]="Baconadjust.pval_Mvalues"
leg3_TOPTABLE=cbind(leg3_topTableFDRannot, leg3_bcp.value, leg3_bcadjust, leg3_corrected_p)

sit2_corrected_p = z_sit2 = data.frame(sit2_topTableFDRannotBETA[match(row.names
      (sit2_topTableFDRannot),row.names(sit2_topTableFDRannotBETA)),])
sit2_bcp.value <- data.frame(sit2_tstats[match(row.names(sit2_corrected_p), row.names(sit2_tstats)),])
sit2_bcadjust <- as.matrix(p.adjust(as.matrix(sit2_tstats), "fdr", n=length(sit2_tstats)))
sit2_tstats <- sit2_tstats[, 1]
row.names(sit2_bcadjust) = row.names(sit2_tstats)
colnames(sit2_bcp.value)[1]="BaconPvalue_Mvalues"
colnames(sit2_bcadjust)[1]="Baconadjust.pval_Mvalues"
sit2_TOPTABLE=cbind(sit2_topTableFDRannot, sit2_bcp.value, sit2_bcadjust, sit2_corrected_p)

sit3_corrected_p = z_sit3 = data.frame(sit3_topTableFDRannotBETA[match(row.names
     (sit3_topTableFDRannot),row.names(sit3_topTableFDRannotBETA)),])
sit3_bcp.value <- data.frame(sit3_tstats[match(row.names(sit3_corrected_p), row.names(sit3_tstats)),])
sit3_bcadjust <- as.matrix(p.adjust(as.matrix(sit3_tstats), "fdr", n=length(sit3_tstats)))
sit3_tstats <- sit3_tstats[, 1]
row.names(sit3_bcadjust) = row.names(sit3_tstats)
colnames(sit3_bcp.value)[1]="BaconPvalue_Mvalues"
colnames(sit3_bcadjust)[1]="Baconadjust.pval_Mvalues"
sit3_TOPTABLE=cbind(sit3_topTableFDRannot, sit3_bcp.value, sit3_bcadjust, sit3_corrected_p)

trunk2_corrected_p = z_trunk2 = data.frame(trunk2_topTableFDRannotBETA[match(row.names
     (trunk2_topTableFDRannot),row.names(trunk2_topTableFDRannotBETA)),])
trunk2_bcp.value <- data.frame(trunk2_tstats[match(row.names(trunk2_corrected_p), row.names(trunk2_tstats)),])
trunk2_bcadjust <- as.matrix(p.adjust(as.matrix(trunk2_tstats), "fdr", n=length(trunk2_tstats)))
trunk2_tstats <- trunk2_tstats[, 1]
row.names(trunk2_bcadjust) = row.names(trunk2_tstats)
colnames(trunk2_bcp.value)[1]="BaconPvalue_Mvalues"
colnames(trunk2_bcadjust)[1]="Baconadjust.pval_Mvalues"
trunk2_TOPTABLE=cbind(trunk2_topTableFDRannot, trunk2_bcp.value, trunk2_bcadjust, trunk2_corrected_p)
trunk3_corrected_p = z_trunk3 = data.frame(trunk3_topTableFDRannotBETA[match(row.names
       (trunk3_topTableFDRannot),row.names(trunk3_topTableFDRannotBETA)),])
trunk3_bcp.value <- data.frame(trunk3_tstats[match(row.names(trunk3_corrected_p), row.names(trunk3_tstats)),])
trunk3_bcadjust <- as.matrix(p.adjust(as.matrix(trunk3_tstats), "fdr", n=length(trunk3_tstats)))
trunk3_tstats <- trunk3_tstats[, 1]
row.names(trunk3_bcadjust) = row.names(trunk3_tstats)
colnames(trunk3_bcp.value)[1]="BaconPvalue_Mvalues"
colnames(trunk3_bcadjust)[1]="Baconadjust.pval_Mvalues"
trunk3_TOPTABLE=cbind(trunk3_topTableFDRannot, trunk3_bcp.value, trunk3_bcadjust, trunk3_corrected_p)



write.table(leg2_TOPTABLE,"Leg Length DMP.txt", sep="\t",col.names=T,row.names=T) 
write.table(leg3_TOPTABLE,"Leg Length DMP (location).txt", sep="\t",col.names=T,row.names=T) 

write.table(sit2_TOPTABLE,"Sitting Height DMP.txt", sep="\t",col.names=T,row.names=T) 
write.table(sit3_TOPTABLE,"Sitting Height DMP (location).txt", sep="\t",col.names=T,row.names=T) 

write.table(trunk2_TOPTABLE,"Trunk Length DMP.txt", sep="\t",col.names=T,row.names=T) 
write.table(trunk3_TOPTABLE,"Trunk Length DMP (location).txt", sep="\t",col.names=T,row.names=T) 



####################################################################################################################################
## DMP FILTERS

leg2mini_file = leg2_TOPTABLE[1:1000,]
leg2mini_sort = cbind(leg2mini_file$P.Value, leg2mini_file$Beta_Delta, leg2mini_file$Mvalue_Delta)
leg2mini_sort=data.frame(leg2mini_sort)
colnames(leg2mini_sort)=c("P.Value", "Beta_Delta", "Mval_Delta")

leg3mini_file = leg3_TOPTABLE[1:1000,]
leg3mini_sort = cbind(leg3mini_file$P.Value, leg3mini_file$Beta_Delta, leg3mini_file$MvalueDelta)
leg3mini_sort=data.frame(leg3mini_sort)
colnames(leg3mini_sort)=c("P.Value", "Beta_Delta", "Mval_Delta")

sit2mini_file = sit2_TOPTABLE[1:1000,]
sit2mini_sort = cbind(sit2mini_file$P.Value, sit2mini_file$Beta_Delta, sit2mini_file$MvalueDelta)
sit2mini_sort=data.frame(sit2mini_sort)
colnames(sit2mini_sort)=c("P.Value", "Beta_Delta", "Mval_Delta")

sit3mini_file = sit3_TOPTABLE
sit3mini_sort = cbind(sit3mini_file$P.Value, sit3mini_file$Beta_Delta, sit3mini_file$MvalueDelta)
sit3mini_sort=data.frame(sit3mini_sort)
colnames(sit3mini_sort)=c("P.Value", "Beta_Delta", "Mval_Delta")

trunk2mini_file = trunk2_TOPTABLE[1:1000,]
trunk2mini_sort = cbind(trunk2mini_file$P.Value, trunk2mini_file$Beta_Delta, trunk2mini_file$MvalueDelta)
trunk2mini_sort=data.frame(trunk2mini_sort)
colnames(trunk2mini_sort)=c("P.Value", "Beta_Delta", "Mval_Delta")

trunk3mini_file = trunk3_TOPTABLE[1:1000,]
trunk3mini_sort = cbind(trunk3mini_file$P.Value, trunk3mini_file$Beta_Delta, trunk3mini_file$MvalueDelta)
trunk3mini_sort=data.frame(trunk3mini_sort)
colnames(trunk3mini_sort)=c("P.Value", "Beta_Delta", "Mval_Delta")


####################################################################################################################################
## DMP FILTERS

leg2_hypo=subset(leg2mini_file, as.numeric(Beta_Delta)<0)
length(rownames(leg2_hypo))
head(rownames(leg2_hypo))
leg2_hyper=subset(leg2mini_file, as.numeric(Beta_Delta)>0)
length(rownames(leg2_hyper))
head(rownames(leg2_hyper))

leg3_hypo=subset(leg3mini_file, as.numeric(Beta_Delta)<0)
length(rownames(leg3_hypo))
head(rownames(leg3_hypo))
leg3_hyper=subset(leg3mini_file, as.numeric(Beta_Delta)>0)
length(rownames(leg3_hyper))
head(rownames(leg3_hyper))

sit2_hypo=subset(si2tmini_file, as.numeric(Beta_Delta)<0)
length(rownames(sit2_hypo))
head(rownames(sit2_hypo))
sit2_hyper=subset(si2tmini_file, as.numeric(Beta_Delta)>0)
length(rownames(sit2_hyper))
head(rownames(sit2_hyper))

si3t_hypo=subset(sit3mini_file, as.numeric(Beta_Delta)<0)
length(rownames(sit3_hypo))
head(rownames(sit3_hypo))
sit3_hyper=subset(sit3mini_file, as.numeric(Beta_Delta)>0)
length(rownames(sit3_hyper))
head(rownames(sit3_hyper))

trunk2_hypo=subset(trunk2mini_file, as.numeric(Beta_Delta)<0)
length(rownames(trunk2_hypo))
head(rownames(trunk2_hypo))
trunk2_hyper=subset(trunk2mini_file, as.numeric(Beta_Delta)>0)
length(rownames(trunk2_hyper))
head(rownames(trunk2_hyper))

trunk3_hypo=subset(trunk3mini_file, as.numeric(Beta_Delta)<0)
length(rownames(trunk3_hypo))
head(rownames(trunk3_hypo))
trunk3_hyper=subset(trunk3mini_file, as.numeric(Beta_Delta)>0)
length(rownames(trunk3_hyper))
head(rownames(trunk3_hyper))


colnames(leg2mini_sort)[1]="P.Value"
leg2.p.value = subset(leg2mini_sort, as.numeric(P.Value)<0.05)
head(leg2.p.value)
length(leg2.p.value)

write.table(leg2_hyper, "Leg Length hypermethylatedDMPs.txt", sep="\t")
write.table(leg2_hypo, "Leg Length hypomethylatedDMPs.txt", sep="\t")

colnames(leg2mini_sort)[1]="P.Value"
leg3.p.value = subset(leg3mini_sort, as.numeric(P.Value)<0.05)
head(leg3.p.value)
length(leg3.p.value)

write.table(leg3_hyper, "Leg Length (location) hypermethylatedDMPs.txt", sep="\t")
write.table(leg3_hypo, "Leg Length (location) hypomethylatedDMPs.txt", sep="\t")

colnames(sit2mini_sort)[1]="P.Value"
sit2.p.value = subset(sit2mini_sort, as.numeric(P.Value)<0.05)
head(sit2.p.value)
length(sit2.p.value)

write.table(sit2_hyper, "Sitting Height hypermethylatedDMPs.txt", sep="\t")
write.table(sit2_hypo, "Sitting Height hypomethylatedDMPs.txt", sep="\t")

colnames(si3tmini_sort)[1]="P.Value"
sit3.p.value = subset(sit3mini_sort, as.numeric(P.Value)<0.05)
head(sit3.p.value)
length(sit3.p.value)

write.table(sit3_hyper, "Sitting Height (location) hypermethylatedDMPs.txt", sep="\t")
write.table(sit3_hypo, "Sitting Height (location) hypomethylatedDMPs.txt", sep="\t")

trunk2.p.value = subset(trunk2mini_sort, as.numeric(P.Value)<0.05)
head(trunk2.p.value)
length(trunk2.p.value)

write.table(trunk2_hyper, "Trunk Height hypermethylatedDMPs.txt", sep="\t")
write.table(trunk2_hypo, "Trunk Height hypomethylatedDMPs.txt", sep="\t")

colnames(trunk2mini_sort)[1]="P.Value"
trunk3.p.value = subset(trunk3mini_sort, as.numeric(P.Value)<0.05)
head(trunk3.p.value)
length(trunk3.p.value)

write.table(trunk3_hyper, "Trunk Height (location) hypermethylatedDMPs.txt", sep="\t")
write.table(trunk3_hypo, "Trunk Height (location) hypomethylatedDMPs.txt", sep="\t")


leg2_sortedtop= cbind(leg2_TOPTABLE$P.Value, leg2_TOPTABLE$adj.P.Val,
     leg2_TOPTABLE$pos, leg2_TOPTABLE$UCSC_RefGene_Name, 
     leg2_TOPTABLE$UCSC_RefGene_Accession, leg2_TOPTABLE$UCSC_RefGene_Group, 
     leg2_TOPTABLE$Relation_to_Island)
dim(leg2_sortedtop)
dim(leg2_TOPTABLE)
colnames(leg2_sortedtop)=c("Pvalue","FDRadj.pval","chr","pos","UCSCGenename",
      "UCSCGeneID","UCSCGENEgroup","Relation_to_Island" )
rownames(leg2_sortedtop)=rownames(leg2_sortedtop)
leg2_sortedtop=data.frame(leg2_sortedtop)
write.table(leg2_sortedtop, "Leg Length sortedDMP.txt", sep="\t")

leg3_sortedtop= cbind(leg3_TOPTABLE$P.Value, leg3_TOPTABLE$adj.P.Val,
       leg3_TOPTABLE$pos, leg3_TOPTABLE$UCSC_RefGene_Name, 
       leg3_TOPTABLE$UCSC_RefGene_Accession, leg3_TOPTABLE$UCSC_RefGene_Group, 
       leg3_TOPTABLE$Relation_to_Island)
dim(leg3_sortedtop)
dim(leg3_TOPTABLE)
colnames(leg3_sortedtop)=c("Pvalue","FDRadj.pval","chr","pos","UCSCGenename",
        "UCSCGeneID","UCSCGENEgroup","Relation_to_Island" )
rownames(leg3_sortedtop)=rownames(leg3_sortedtop)
leg3_sortedtop=data.frame(leg3_sortedtop)
write.table(leg3_sortedtop, "Leg Length sortedDMP2.txt", sep="\t")

sit2_sortedtop= cbind(sit2_TOPTABLE$P.Value, sit2_TOPTABLE$adj.P.Val,
       sit2_TOPTABLE$pos, sit2_TOPTABLE$UCSC_RefGene_Name, 
       sit2_TOPTABLE$UCSC_RefGene_Accession, sit2_TOPTABLE$UCSC_RefGene_Group, 
       sit2_TOPTABLE$Relation_to_Island)
dim(sit2_sortedtop)
dim(sit2_TOPTABLE)
colnames(sit2_sortedtop)=c("Pvalue","FDRadj.pval","chr","pos","UCSCGenename",
      "UCSCGeneID","UCSCGENEgroup","Relation_to_Island" )
rownames(sit2_sortedtop)=rownames(sit2_sortedtop)
sit2_sortedtop=data.frame(sit2_sortedtop)
write.table(sit2_sortedtop, "Sitting Height sortedDMP.txt", sep="\t")

sit3_sortedtop= cbind(sit3_TOPTABLE$P.Value, sit3_TOPTABLE$adj.P.Val,
      sit3_TOPTABLE$pos, sit3_TOPTABLE$UCSC_RefGene_Name, 
      sit3_TOPTABLE$UCSC_RefGene_Accession, sit3_TOPTABLE$UCSC_RefGene_Group, 
      sit3_TOPTABLE$Relation_to_Island)
dim(sit3_sortedtop)
dim(sit3_TOPTABLE)
colnames(sit3_sortedtop)=c("Pvalue","FDRadj.pval","chr","pos","UCSCGenename",
      "UCSCGeneID","UCSCGENEgroup","Relation_to_Island" )
rownames(sit3_sortedtop)=rownames(sit3_sortedtop)
sit3_sortedtop=data.frame(sit3_sortedtop)
write.table(sit3_sortedtop, "Sitting Height sortedDMP2.txt", sep="\t")

trunk2_sortedtop= cbind(trunk2_TOPTABLE$P.Value, trunk2_TOPTABLE$adj.P.Val,
      trunk2_TOPTABLE$pos, trunk2_TOPTABLE$UCSC_RefGene_Name, 
      trunk2_TOPTABLE$UCSC_RefGene_Accession, trunk2_TOPTABLE$UCSC_RefGene_Group, 
      trunk2_TOPTABLE$Relation_to_Island)
dim(trunk2_sortedtop)
dim(trunk2_TOPTABLE)
colnames(trunk2_sortedtop)=c("Pvalue","FDRadj.pval","chr","pos","UCSCGenename",
      "UCSCGeneID","UCSCGENEgroup","Relation_to_Island" )
rownames(trunk2_sortedtop)=rownames(trunk2_sortedtop)
trunk2_sortedtop=data.frame(trunk2_sortedtop)
write.table(trunk2_sortedtop, "Trunk Height sortedDMP.txt", sep="\t")

trunk3_sortedtop= cbind(trunk3_TOPTABLE$P.Value, trunk3_TOPTABLE$adj.P.Val,
       trunk3_TOPTABLE$pos, trunk3_TOPTABLE$UCSC_RefGene_Name, 
       trunk3_TOPTABLE$UCSC_RefGene_Accession, trunk3_TOPTABLE$UCSC_RefGene_Group, 
       trunk3_TOPTABLE$Relation_to_Island)
dim(trunk3_sortedtop)
dim(trunk3_TOPTABLE)
colnames(trunk3_sortedtop)=c("Pvalue","FDRadj.pval","chr","pos","UCSCGenename",
        "UCSCGeneID","UCSCGENEgroup","Relation_to_Island" )
rownames(trunk3_sortedtop)=rownames(trunk3_sortedtop)
trunk3_sortedtop=data.frame(trunk3_sortedtop)
write.table(trunk3_sortedtop, "Trunk Height sortedDMP2.txt", sep="\t")



####################################################################################################################################
## CODE FOR DETECTING DMR USING BUMPHUNTER

library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(minfiData)
library(limma)
library(affy)
library(DMRcate)
library(bumphunter)

options(stringsAsFactors = F)

rownames(leg_modela) <- targets$IDs
head(leg_modela)
leg.anno = getAnnotation(raw_set)
leg.anno.range = makeGRangesFromDataFrame(leg.anno, keep.extra.columns = T,
    start.field = "pos", end.field = "pos")
bump_nostat <- bumphunter(GMset.culled, leg_modela, coef=2 ,cutoff=0.1)
bump <- bump_nostat$table
head(bump)
bump1 = subset(bump, bump$L > 2)
length(bump1)
bumps.range <- makeGRangesFromDataFrame(bump1, keep.extra.columns = TRUE)
bumps.anno = leg.anno.range[nearest(x = bumps.range, subject = leg.anno.range),]
bump2 = cbind(bump1, Gene_start = as.data.frame(ranges(bumps.anno))$start, 
      Gene_end = as.data.frame(ranges(bumps.anno))$end, Gene_Name = 
      bumps.anno$UCSC_RefGene_Name)
write.table(bump2,"Leg_bumphunterDMRs.txt",sep="\t")

rownames(sit_modela) <- targets$IDs
head(sit_modela)
sit.anno = getAnnotation(raw_set)
sit.anno.range = makeGRangesFromDataFrame(sit.anno, keep.extra.columns = T,
    start.field = "pos", end.field = "pos")
bump_nostat <- bumphunter(GMset.culled, sit_modela, coef=2 ,cutoff=0.1)
bump <- bump_nostat$table
head(bump)
bump1 = subset(bump, bump$L > 2)
length(bump1)
bumps.range <- makeGRangesFromDataFrame(bump1, keep.extra.columns = TRUE)
bumps.anno = sit.anno.range[nearest(x = bumps.range, subject = sit.anno.range),]
bump2 = cbind(bump1, Gene_start = as.data.frame(ranges(bumps.anno))$start, 
     Gene_end = as.data.frame(ranges(bumps.anno))$end, Gene_Name = 
     bumps.anno$UCSC_RefGene_Name)
write.table(bump2,"Sitting_bumphunterDMRs.txt",sep="\t")

rownames(trunk_modela) <- targets$IDs
head(trunk_modela)
trunk.anno = getAnnotation(raw_set)
trunk.anno.range = makeGRangesFromDataFrame(trunk.anno, keep.extra.columns = T,
      start.field = "pos", end.field = "pos")
bump_nostat <- bumphunter(GMset.culled, trunk_model, coef=2 ,cutoff=0.1)
bump <- bump_nostat$table
head(bump)
bump1 = subset(bump, bump$L > 2)
length(bump1)
bumps.range <- makeGRangesFromDataFrame(bump1, keep.extra.columns = TRUE)
bumps.anno = trunk.anno.range[nearest(x = bumps.range, subject = trunk.anno.range),]
bump2 = cbind(bump1, Gene_start = as.data.frame(ranges(bumps.anno))$start, 
      Gene_end = as.data.frame(ranges(bumps.anno))$end, Gene_Name = 
      bumps.anno$UCSC_RefGene_Name)
write.table(bump2,"Trunk_bumphunterDMRs.txt",sep="\t")
