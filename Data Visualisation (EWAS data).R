
## THE GRAPHS AND PLOTS GENERATED FROM THE EWAS
#As the EWAS script is proving to be a very long document I have decided to record
#my graph and plot scripts separately as to keep the EWAS script just the
#necessary code

#################################################################################################
## CODE THE MULTIPLE VARIABLES INTO SEVERAL PCAS AND PCA CORRELATIONS

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data/1 EWAS Output Files")

library(ggplot2)
library(limma)
library(vsn)
library(affy)
library(fields)
library(preprocessCore)
library(lumi)
library(svd)
library (reshape)
library(ggplot2)
library(Cairo)


idats_RODAM = "~/lkg/Rodam/idats/"
idat_phenofile = read.csv("PhenoRodam15032106Analyse1.csv", sep=",",header=T)
targets = read.csv("Rodam_Data_Simplified.csv", sep=",",header=T)
raw_betaRODAM = read.table("Raw Beta Values.txt",sep="\t", header = T)
cell_counts = read.table("Cell Type Table.txt",sep="\t", header = T)
mvalue = read.table("RODAM M value.txt", sep="\t", header = T)


targets = cbind(targets, cell_counts)
attach(targets)

bmigcpopulation = length(rownames(targets))

M.culled.centered <- mvalue-rowMeans(mvalue)
M.svd <- svd(na.omit(M.culled.centered))

a = svd(na.omit(raw_betaRODAM))
b=a$v 

PCA_data <- prcomp(raw_betaRODAM)
summary(PCA_data)

cd8t = targets$CD8T
cd4t = targets$CD4T
mono= targets$Mono
nk = targets$NK
bcell= targets$Bcell
gran= targets$Gran

pdf("PCA_gender.pdf",
    for (i in 1:3){
      for (j in 1:3){
        plot (b[,i],b[,j],main ="PCA", xlab=i, ylab=j, pch=1, col="white")
        for (l in 1:population ){
          if (targets$gender[l]== "1"){
            text(b[l,i],b[l,j],targets$Sample_ID[l], col="blue")
          }
          else if (targets$gender[2]== "2"){
            text(b[l,i],b[l,j],targets$Sample_ID[l], col="pink")
          }
          else{
            text(b[l,i],b[l,j],targets$Sample_ID[l], col="grey")
          }}}})
dev.off()


pdf("PCA_location.pdf")
for (i in 1:3){
  for (j in 1:3){
    plot (b[,i],b[,j],main ="PCA", xlab=i, ylab=j, pch=1, col="white")
    for (l in 1:population ){  
      if (targets$Location[1]== "Amsterdam"){
        text(b[l,i],b[l,j],targets$Sample_ID[l], col="red")
        else if (targets$Location[2]== "Berlin"){
          text(b[l,i],b[l,j],targets$Sample_ID[l], col="orange")
          else if (targets$Location[3]== "London"){
            text(b[l,i],b[l,j],targets$Sample_ID[l], col="yellow")
            else if (targets$Location[4]== "Urban Ghana"){
              text(b[l,i],b[l,j],targets$Sample_ID[l], col="blue")
              else if (targets$Location[5]== "Rural Ghana"){
                text(b[l,i],b[l,j],targets$Sample_ID[l], col="purple")
              }}}}}}}}
dev.off()

pdf("PCA_Diabetes.pdf")
for (i in 1:3){
  for (j in 1:3){
    plot (b[,i],b[,j],main ="PCA", xlab=i, ylab=j, pch=1, col="white")
    for (l in 1:population ){
      if (targets$DiabetesYN[l]== "Yes"){
        text(b[l,i],b[l,j],targets$Sample_ID[l], col="blue")
      }
      else if (targets$DiabetesYN[2]== "No"){
        text(b[l,i],b[l,j],targets$Sample_ID[l], col="yellow")
      }
      else if (targets$DiabetesYN[3]== "Unknown"){
        text(b[l,i],b[l,j],targets$Sample_ID[l], col="green")
      }}}}
dev.off()

pdf("PCA_Hypertension.pdf")
for (i in 1:3){
  for (j in 1:3){
    plot (b[,i],b[,j],main ="PCA", xlab=i, ylab=j, pch=1, col="white")
    for (l in 1:population ){
      if (targets$Hyperten_Meds.2[l]== "Yes"){
        text(b[l,i],b[l,j],targets$Sample_ID[l], col="green")
      }
      else if (targets$Hyperten_Meds.2[0]== "No"){
        text(b[l,i],b[l,j],targets$Sample_ID[l], col="orange")
      }}}}
dev.off()

#Code for the PCA correlations which will be divided by demographic variables,
##physical proxies, secondary conditions, and then all the variables together.

age.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~age)
  return(summary(fit)$adj.r.squared)
})

gender.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~as.factor(sex))
  return(summary(fit)$adj.r.squared)
})

locate.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~location)
  return(summary(fit)$adj.r.squared)
})


leg_length.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~leg)
  return(summary(fit)$adj.r.squared)
})

SH.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~sitting)
  return(summary(fit)$adj.r.squared)
})

trunk.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~trunk)
  return(summary(fit)$adj.r.squared)
})

WHR.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~whr)
  return(summary(fit)$adj.r.squared)
})

diabetes.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~diabetes)
  return(summary(fit)$adj.r.squared)
})

hyper.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~hyperten)
  return(summary(fit)$adj.r.squared)
})

obe.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~obesity)
  return(summary(fit)$adj.r.squared)
})



#Cell Type
cd8t.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~cd8t)
  return(summary(fit)$adj.r.squared)
})


cd4t.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~cd4t)
  return(summary(fit)$adj.r.squared)
})

nk.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~nk)
  return(summary(fit)$adj.r.squared)
})

bcell.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~bcell)
  return(summary(fit)$adj.r.squared)
})

mono.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~mono)
  return(summary(fit)$adj.r.squared)
})

gran.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~gran)
  return(summary(fit)$adj.r.squared)
})

pos.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~pos)
  return(summary(fit)$adj.r.squared)
})


batchArry.corr <- sapply(1:ncol(M.svd$v), function(i){
  fit <- lm(M.svd$v[,i]~batchA)
  return(summary(fit)$adj.r.squared)
})

var_per_PC.melt <- melt(M.svd$d^2/sum(M.svd$d^2))
var_per_PC.melt <- cbind( var_per_PC.melt, (1:length(M.svd$d)))
colnames(var_per_PC.melt) <- c("Variance", "PC")

#Generate the PCA correlation plot and Variance Plot

Cairo(file = "~/lkg/Rodam/GALATEA/Data/Test Variance_PC.jpeg", type = "jpeg", units = "px", width = 800, height = 800, dpi = 90, bg = "white")
ggplot(var_per_PC.melt, aes(x = PC, y = Variance, ymax = 0.15)) +
  geom_point(size = 1) +
  theme_bw() +
  ggtitle("Variance explained per PC") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
dev.off

Cairo(file = "Test demographic_PCA.jpeg", type = "jpeg", units = "px",
      width = 800, height = 800, dpi = 110, bg = "gray97")
plot(head(gender.corr, n=8L),pch=8, col="gold",ylim=c(-1,1),xlab="PC's",
     ylab="correlation")
points(head(age.corr, n=8L),pch=5, col="pink")
points(head (locate.corr, n=8L),pch=2, col="darkorchid4")
legend("bottomleft", c("Gender", "Age", "Location"),col=c("gold", "plum1", 
                                                          "darkorchid4"), pch =c(8, 5, 2))
dev.off()

Cairo(file = "Test proxy_PCA.jpeg", type = "jpeg", units = "px",
      width = 800, height = 800, dpi = 110, bg = "gray97")
plot(head(leg_length.corr, n=8L),pch=8, col="seagreen1",ylim=c(-1,1),xlab="PC's",
     ylab="correlation")
points(head(WHR.corr, n=8L),pch=5, col="cadetblue1")
points(head(SH.corr, n=8L),pch=5, col="turquoise4")
points(head (trunk.corr, n=8L),pch=2, col="dodgerblue4")
legend("bottomleft", c("Leg Length", "Waist-to-Hip Ratio", "Sitting Height", 
                       "Trunk Length"),col=c("seagreen1", "cadetblue1", "turquoise4", "dodgerblue4"), pch =c(8, 9, 5, 2))
dev.off()

Cairo(file = "Test health_PCA.jpeg", type = "jpeg", units = "px",
      width = 800, height = 800, dpi = 110, bg = "gray97")
plot(head(diabetes.corr, n=8L),pch=8, col="firebrick1",ylim=c(-1,1),xlab="PC's",
     ylab="correlation")
points(head(obe.corr, n=8L),pch=5, col="plum1")
points(head (hyper.corr, n=8L),pch=2, col="dodgerblue4")
legend("bottomleft", c("Diabetes Diagnosis", "Obesity", "Hypertension"),
       col=c("firebrick1", "plum1", "dodgerblue4"), pch =c(8, 5, 2))
dev.off()

Cairo(file = "Test All Variables_PCA.jpeg", type = "jpeg", units = "px",
      width = 800, height = 800, dpi = 110, bg = "white")
plot(head(gender.corr, n=8L),pch=1, col="gold",ylim=c(-0.2,0.2),xlab="PC's",
     ylab="correlation") 
points(head(age.corr, n=8L),pch=2, col="plum1") 
points(head (locate.corr, n=8L),pch=3, col="darkorchid4") 
points(head(diabetes.corr, n=8L),pch=4, col="firebrick1") 
points(head(hyper.corr, n=8L),pch=5, col="sienna1") 
points(head(obesity.corr, n=8L),pch=6, col="deeppink3") 
points(head(cd8t.corr, n=8L),pch=7, col="darkgoldenrod1") 
points(head(cd4t.corr, n=8L),pch=8, col="red3") 
points(head(nk.corr, n=8L),pch=9, col="forestgreen") 
points(head(bcell.corr, n=8L),pch=1, col="dodgerblue1") 
points(head(mono.corr, n=8L),pch=2, col="purple3") 
points(head(batchArry.corr, n=8L),pch=3, col="violetred1") 
points(head(pos.corr, n=8L),pch=4, col="springgreen3") 
legend("right", c("Gender", "Age", "Location", "Diabetes Diagnosis", "Hypertension", 
      "Obesity", "CD8T", "CD4T", "Natural Killer", "B", "Mono", "Batch Array", "Plate Position"), 
      col=c("gold", "plum1", "darkorchid4", "firebrick1","sienna1", "deeppink3", "darkgoldenrod1",
      "red3", "forestgreen",  "dodgerblue1", "purple3", "violetred1", "springgreen3"), 
       pch =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4))
dev.off()

Cairo(file = "Cell_Type_PCA.jpeg", type = "jpeg", units = "px",
      width = 800, height = 800, dpi = 110, bg = "gray97")
plot(head(cd8t.corr, n=8L),pch=1, col="darkgoldenrod1",ylim=c(-1,1),xlab="PC's",
     ylab="correlation")
points(head(cd4t.corr, n=8L),pch=2, col="red3")
points(head(nk.corr, n=8L),pch=5, col="forestgreen")
points(head(bcell.corr, n=8L),pch=6, col="dodgerblue1")
points(head(mono.corr, n=8L),pch=8, col="purple3")
legend("bottomleft", c("CD8T", "CD4T", "Natrual Killer cell", "B Cell", "Mono"),
       col=c("darkgoldenrod1", "red3", "forestgreen", "dodgerblue1", "purple3"), pch =c(1, 2, 5, 6, 8))
dev.off()

Cairo(file = "Compotents.jpeg", type = "jpeg", units = "px",
      width = 800, height = 800, dpi = 110, bg = "gray97")
plot(head(batchArry.corr, n=8L),pch=8, col="violetred1",ylim=c(-1,1),xlab="PC's",
     ylab="correlation")
points(head(pos.corr, n=8L),pch=5, col="springgreen3")
legend("bottomleft", c("Batch Array", "Position"),
       col=c("violetred1", "springgreen3"), pch =c(8, 5, 2))
dev.off()

#################################################################################################
## VOLCANO PLOTS

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

library(ggplot2)
install.packages("viridis")
library(viridis)

leg_TOPTABLE = read.table("Leg Length DMP.txt", sep="\t", header = T)
sit_TOPTABLE = read.table("Sitting Height DMP.txt", sep="\t", header = T)
trunk_TOPTABLE = read.table("Trunk Length DMP.txt", sep="\t", header = T)
whr_TOPTABLE = read.table("WHR with BMI DMP.txt", sep="\t", header = T)

#Each individual plot

volplot <- data.frame(logFC = leg_TOPTABLE$Beta.logF, Padj = leg_TOPTABLE$P.Value)
rownames(volplot) <- rownames(leg_TOPTABLE)
volplot$threshold <- (volplot$Padj < 10E-6 & abs(volplot$logFC) > 0.015)*1
volplot$threshold[which(volplot$threshold == 1)] <- "Significant"
volplot$threshold[which(volplot$threshold == 0)] <- "Non significant"
volplot$threshold <- as.factor(volplot$threshold)
volplot$threshold <- relevel(volplot$threshold, ref= "Significant")
jpeg("Vol-Leg Length.jpeg")
ggplot(volplot, aes(x = logFC, y = -log10(Padj), color = threshold)) + 
  geom_point(alpha = 0.4, size = 2) +
  geom_vline(xintercept = -0.002, linetype = "longdash") +
  geom_vline(xintercept = 0.002, linetype = "longdash") +
  geom_hline(yintercept = -log10(10E-6), linetype = "longdash") +
  
  theme_bw() +
  #xlim(c(-0.002, 0.002)) +
  ylab("-log10(pval)") +
  xlab("Mean effect size (Beta)") +
  scale_color_jcolors(palette = "default") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()


volplot <- data.frame(logFC = sit_TOPTABLE$Beta.logFC, Padj = sit_TOPTABLE$P.Value)
volplot$threshold <- (volplot$Padj < 10E-6 & abs(volplot$logFC) > 0.015)*1
volplot$threshold[which(volplot$threshold == 1)] <- "Significant"
volplot$threshold[which(volplot$threshold == 0)] <- "Non significant"
volplot$threshold <- as.factor(volplot$threshold)
volplot$threshold <- relevel(volplot$threshold, ref = 1)
jpeg("Vol-Sitting Height.jpg")
ggplot(volplot, aes(x = logFC, y = -log10(Padj), color = threshold)) + 
  geom_point(alpha = 0.4, size = 1.5) +
  geom_vline(xintercept = -0.002, linetype = "longdash") +
  geom_vline(xintercept = 0.002, linetype = "longdash") +
  geom_hline(yintercept = -log10(10E-6), linetype = "longdash") +
  
  theme_bw() +
  #xlim(c(-0.002, 0.002)) +
  ylab("-log10(pval)") +
  xlab("Mean effect size (Beta)") +
  scale_color_jcolors(palette = "default") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()


volplot <- data.frame(logFC = trunk_TOPTABLE$Beta.logF, Padj = trunk_TOPTABLE$P.Value)
rownames(volplot) <- rownames(trunk_TOPTABLE)
volplot$threshold <- (volplot$Padj < 10E-6 & abs(volplot$logFC) > 0.015)*1
volplot$threshold[which(volplot$threshold == 1)] <- "Significant"
volplot$threshold[which(volplot$threshold == 0)] <- "Non significant"
volplot$threshold <- as.factor(volplot$threshold)
volplot$threshold <- relevel(volplot$threshold, ref = 1)
jpeg("Vol-Trunk Height.jpg")
ggplot(volplot, aes(x = logFC, y = -log10(Padj), color = threshold)) + 
  geom_point(alpha = 0.4, size = 2) +
  geom_vline(xintercept = -0.002, linetype = "longdash") +
  geom_vline(xintercept = 0.002, linetype = "longdash") +
  geom_hline(yintercept = -log10(10E-6), linetype = "longdash") +
  
  theme_bw() +
  #xlim(c(-0.002, 0.002)) +
  ylab("-log10(pval)") +
  xlab("Mean effect size (Beta)") +
  scale_color_jcolors(palette = "default") +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()


#################################################################################################
## ManhattanPLOTS

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

library(ggplot2)

Obe = leg_TOPTABLE
Obe$chr=gsub("chr","",Obe$chr)
Obe$chr=as.numeric(Obe$chr)
Obe$pos=as.numeric(Obe$pos)
Obe <- Obe[!is.na(Obe$pos), ]
Obe <- Obe[Obe$pos > 0, ]
Obe <- Obe[is.element(Obe$chr, 1:22), ]
cols <- c()
coor <- c()
.col <- c("seagreen3", "skyblue3")
.oldMaxCoor <- 0
ticks <- c()
p <- c()
for (chr in 1:22) {
  newRes <- Obe[which(Obe$chr == chr), ]
  p <- c(p, newRes$P.Value)
  coor <- c(coor, newRes$pos+ .oldMaxCoor)      
  cols <- c(cols, rep(.col[chr %% 2 + 1], nrow(newRes)))
  ticks <- c(ticks, mean(c(.oldMaxCoor, max(coor[!is.na(coor)]))))
  .oldMaxCoor <- max(coor[!is.na(coor)]) #*
}

jpeg("Manhatta-LEG.jpeg", height = 600,
    width = 1600, res = 96)
plot(coor,
     -log10(p), type = "p",
     col = cols, pch = 20,
     xlab = "Chromosome", ylab = "-log10(Pvalue)",
     cex.lab=1.5, cex.axis=1.5,
     ylim = c(0, 15), xaxt = "n")
axis(1, ticks, 1:22)
abline(a= NULL, b= NULL, h = -log10(1E-7), col = "red", lty = 3)
dev.off()

Obe = sit_TOPTABLE
Obe$chr=gsub("chr","",Obe$chr)
Obe$chr=as.numeric(Obe$chr)
Obe$pos=as.numeric(Obe$pos)
Obe <- Obe[!is.na(Obe$pos), ]
Obe <- Obe[Obe$pos > 0, ]
Obe <- Obe[is.element(Obe$chr, 1:22), ]
cols <- c()
coor <- c()
.col <- c("seagreen3", "skyblue3")
.oldMaxCoor <- 0
ticks <- c()
p <- c()
for (chr in 1:22) {
  newRes <- Obe[which(Obe$chr == chr), ]
  p <- c(p, newRes$P.Value)
  coor <- c(coor, newRes$pos+ .oldMaxCoor)      
  cols <- c(cols, rep(.col[chr %% 2 + 1], nrow(newRes)))
  ticks <- c(ticks, mean(c(.oldMaxCoor, max(coor[!is.na(coor)]))))
  .oldMaxCoor <- max(coor[!is.na(coor)]) #*
}

jpeg("Manhattan-SIT.jpeg", height = 600,
     width = 1600, res = 96)
plot(coor,
     -log10(p), type = "p",
     col = cols, pch = 20,
     xlab = "Chromosome", ylab = "-log10(Pvalue)",
     cex.lab=1.5, cex.axis=1.5,
     ylim = c(0, 15), xaxt = "n")
axis(1, ticks, 1:22)
abline(a= NULL, b= NULL, h = -log10(1E-7), col = "red", lty = 3)
dev.off()

Obe = trunk_TOPTABLE
Obe$chr=gsub("chr","",Obe$chr)
Obe$chr=as.numeric(Obe$chr)
Obe$pos=as.numeric(Obe$pos)
Obe <- Obe[!is.na(Obe$pos), ]
Obe <- Obe[Obe$pos > 0, ]
Obe <- Obe[is.element(Obe$chr, 1:22), ]
cols <- c()
coor <- c()
.col <- c("seagreen3", "skyblue3")
.oldMaxCoor <- 0
ticks <- c()
p <- c()
for (chr in 1:22) {
  newRes <- Obe[which(Obe$chr == chr), ]
  p <- c(p, newRes$P.Value)
  coor <- c(coor, newRes$pos+ .oldMaxCoor)      
  cols <- c(cols, rep(.col[chr %% 2 + 1], nrow(newRes)))
  ticks <- c(ticks, mean(c(.oldMaxCoor, max(coor[!is.na(coor)]))))
  .oldMaxCoor <- max(coor[!is.na(coor)]) #*
}

jpeg("Manhattan-TRUNK.jpeg", height = 600,
     width = 1600, res = 96)
plot(coor,
     -log10(p), type = "p",
     col = cols, pch = 20,
     xlab = "Chromosome", ylab = "-log10(Pvalue)",
     cex.lab=1.5, cex.axis=1.5,
     ylim = c(0, 15), xaxt = "n")
axis(1, ticks, 1:22)
abline(a= NULL, b= NULL, h = -log10(1E-7), col = "red", lty = 3)
dev.off()


#################################################################################################
## QQ PLOTS

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

x=function (pvector, ...) 
{
  if (!is.numeric(pvector)) 
    stop("Input must be numeric.")
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & 
           is.finite(pvector) & pvector < 1 & pvector > 0]
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  def_args <- list(pch = 20, xlim = c(0, max(e)), ylim = c(0, 
            max(o)), xlab = expression(Expected ~ ~-log[10](italic(p))), 
            ylab = expression(Observed ~ ~-log[10](italic(p))))
  dotargs <- list(...)
  tryCatch(do.call("plot", c(list(x = e, y = o), def_args[!names(def_args) %in% 
            names(dotargs)], dotargs)), warn = stop)
  abline(0, 1, col = "#B41A21", lwd=3)
}

qqT2D = leg_TOPTABLE
pdf("Leg QQplot.jpeg")
par(mai=c(1.6,2.5,0.5,0.5),mgp=c(7,2,0))
x(qqT2D$P.Value, pch = 2, frame.plot=F, cex = 0.5, cex.lab=1, cex.axis=1, bg="white", ylim=c(0, 12))

par(mai=c(1.6,2.5,0.5,0.5),mgp=c(7,2,0))
x(qqT2D$BaconPvalue_Mvalues, pch = 2, frame.plot=F, cex = 0.5, cex.lab=1, cex.axis=1, bg="white", ylim=c(0, 12))
dev.off()

qqT2D = sit_TOPTABLE
pdf("Sitting Height QQplot.jpeg")
par(mai=c(1.6,2.5,0.5,0.5),mgp=c(7,2,0))
x(qqT2D$P.Value, pch = 2, frame.plot=F, cex = 0.5, cex.lab=1, cex.axis=1, bg="white", ylim=c(0, 12))

par(mai=c(1.6,2.5,0.5,0.5),mgp=c(7,2,0))
x(qqT2D$BaconPvalue_Mvalues, pch = 2, frame.plot=F, cex = 0.5, cex.lab=1, cex.axis=1, bg="white", ylim=c(0, 12))
dev.off()

qqT2D = trunk_TOPTABLE
pdf("Trunk Length QQplot.jpeg")
par(mai=c(1.6,2.5,0.5,0.5),mgp=c(7,2,0))
x(qqT2D$P.Value, pch = 2, frame.plot=F, cex = 0.5, cex.lab=1, cex.axis=1, bg="white", ylim=c(0, 12))

par(mai=c(1.6,2.5,0.5,0.5),mgp=c(7,2,0))
x(qqT2D$BaconPvalue_Mvalues, pch = 2, frame.plot=F, cex = 0.5, cex.lab=1, cex.axis=1, bg="white", ylim=c(0, 12))
dev.off()

gc()
