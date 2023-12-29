
#This script demonstrates one of the preparation of the RODAM data for my EWAS pipelines.
#all code was self-written and not drafted from a template pipeline.

#################################################################################################
## THE METHYLAID AND VISUALISATION CODE

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

library(foreign)
library(dplyr)
library(tidyverse)
library(tidyr)
library(haven)

setwd("/mnt/smb/lkg/gaswart/Rodam/GALATEA/Data")

RODAMdata <- read.sav("Data.sav")
attach(RODAMdata)
dim(RODAMdata)


data_simple <- RODAMdata %>% drop_na(Epigenetics)
summary(data_simple)
dim(data_simple)

##################################################################################################
#Include all your epigenetic phenofile information to your phenotype file

RODAMdata <- read_sav("Data.sav")
data_simple <- RODAMdata %>% drop_na(Epigenetics)

idat_phenofile = read.csv("PhenoRodam15032106Analyse.csv", sep=",",header=T)

DATA <-merge(data_simple, idat_phenofile, by="RodamID")

write.csv(DATA, "Epigenetic Raw Rodam Data.csv")

data_simple <- read.csv("Epigenetic Raw Rodam Data.csv", sep=",",header=T)

##################################################################################################
#Drop all unnecessary columns and variables 
data_simple = subset(data_simple, select = -c(X, Illnessd, Illnesse, Illnessf, Illnessg, 
   Illnessh, Illnessi, Illnessj, BPtreat, BPmed, BPdiet, DiabGest, DiabAge, DiabTreat, DiabTabl, 
   DiabDiet, DiabInsImm, DiabIns6m, Reasonmigr1, Reasonmigr2, Reasonmigr3, Reasonmigr4, 
   Reasonmigr5, Reasonmigr6, Reasonmigr7, Reasonmigr8, Reasonmigr9, ReasonmigrX, LiveGHA,
   LiveCityGHAX, LiveTime, GroupGHAX, Countryfat, CountryfatX))
dim(data_simple)


##################################################################################################
#I am going to rename my columns to help me knowing what variable is what
colnames(data_simple) [2] <- "Location"
colnames(data_simple) [3] <- "Sex"
colnames(data_simple) [4] <- "Age"
colnames(data_simple) [5] <- "Diabetes_Meds."
colnames(data_simple) [6] <- "Hyperten_Meds."
colnames(data_simple) [7] <- "Diuretics"
colnames(data_simple) [8] <- "Beta_Blockers"
colnames(data_simple) [9] <- "Calcium_Blockers"
colnames(data_simple) [10] <-"ACE/ARB"
colnames(data_simple) [11] <-"Hyperten_Meds.2"
colnames(data_simple) [12] <-"Sitting_Height"
colnames(data_simple) [13] <-"Sitting_Height2"
colnames(data_simple) [14] <-"Sitting_Height3" 
colnames(data_simple) [15] <-"Sit._Difference"
colnames(data_simple) [16] <-"Sittin_Height4"
colnames(data_simple) [17] <-"Sitting_Height_Remarks"
colnames(data_simple) [18] <-"Sitting_Calculation"
colnames(data_simple) [19] <-"Leg_Length"
colnames(data_simple) [20] <-"Trunk_Length"
colnames(data_simple) [21] <-"Height_Y/N"
colnames(data_simple) [22] <-"Height1"
colnames(data_simple) [23] <-"Height2"
colnames(data_simple) [24] <-"Height_Difference"
colnames(data_simple) [25] <-"Height3"
colnames(data_simple) [26] <-"Height Remarks"
colnames(data_simple) [27] <-"Height_Calculation"
colnames(data_simple) [28] <-"Weight_Calculation"
colnames(data_simple) [29] <-"Waist_Y/N"
colnames(data_simple) [30] <-"Hip_Y/N"
colnames(data_simple) [31] <-"Waist_Measure.1"
colnames(data_simple) [32] <-"Hip_Measure.1"
colnames(data_simple) [33] <-"Waist_Measure.2"
colnames(data_simple) [34] <-"Hip_Measure.2"
colnames(data_simple) [35] <-"Waist_Difference"
colnames(data_simple) [36] <-"Hip_Difference"
colnames(data_simple) [37] <-"Waist_Measure.3"
colnames(data_simple) [38] <-"Hip_Measure.3"
colnames(data_simple) [39] <-"Waist_Calculation"
colnames(data_simple) [40] <-"Hip_Calculation"
colnames(data_simple) [41] <-"WHR"
colnames(data_simple) [42] <-"General_Health"
colnames(data_simple) [43] <-"CVD_Family"
colnames(data_simple) [44] <-"CVD_Father"
colnames(data_simple) [45] <-"Age_CVD_Father"
colnames(data_simple) [46] <-"CVD_Mother"
colnames(data_simple) [47] <-"Age_CVD_Mother"
colnames(data_simple) [48] <-"CVD_Brother"
colnames(data_simple) [49] <-"Age_CVD_Brother"
colnames(data_simple) [50] <-"CVD_Sister"
colnames(data_simple) [51] <-"Age_Sister_CVD"
colnames(data_simple) [52] <-"CVD_Son"
colnames(data_simple) [53] <-"Age_SON_CVD"
colnames(data_simple) [54] <-"CVD_Daughter"
colnames(data_simple) [55] <-"Age_Daughter_CVD"
colnames(data_simple) [56] <-"Stroke"
colnames(data_simple) [57] <-"Heart"
colnames(data_simple) [58] <-"Heart_Condition"
colnames(data_simple) [59] <-"Blood_PressureYN"
colnames(data_simple) [60] <-"Diabetes"
colnames(data_simple) [61] <-"DiabetesYN"
colnames(data_simple) [62] <-"Diabetes_InsYN"
colnames(data_simple) [63] <-"Diabetes_Family"
colnames(data_simple) [64] <-"Birth_Country"
colnames(data_simple) [65] <-"Second_Birth_Country"
dim(data_simple)


##################################################################################################
#The next line shows me dropping all the missing nodes from my necessary columns
data_simple <- data_simple %>%
  mutate(Leg_Length = replace(Leg_Length, Leg_Length == "NA", NA)) %>%
  mutate(Height_Calculation = replace(Height_Calculation, Height_Calculation == "NA", NA)) %>%
  mutate(Age = replace(Age, Age == "NA", NA)) %>%
  mutate(Sitting_Calculation = replace(Sitting_Calculation, Sitting_Calculation == "NA", NA)) %>%
  mutate(Trunk_Length = replace(Trunk_Length, Trunk_Length == "NA", NA)) %>%
  mutate(WHR = replace(WHR, WHR == "NA", NA)) %>%
  tidyr:: drop_na(Leg_Length, Height_Calculation, Age, Sitting_Calculation, Trunk_Length, WHR)

any(is.na(data_simple$Age))
any(is.na(data_simple$Sitting_Calculation))
any(is.na(data_simple$Leg_Length))
any(is.na(data_simple$Trunk_Length))
any(is.na(data_simple$WHR))
any(is.na(data_simple$CVD_Family))


####################################################################################
#Now that I am aware of the class type of each variable I will convert those variables
#I believe are not in their proper class. Given that many of the variables are read in
#correctly, I might also leave my data as it is. 
data_simple$Sex = as.factor(data_simple$Sex)
data_simple$Location = as.factor(data_simple$Location) 
data_simple$Diabetes_Meds. = as.factor(data_simple$Diabetes_Meds.)
data_simple$Hyperten_Meds. = as.factor(data_simple$Hyperten_Meds.)
data_simple$Diuretics = as.factor(data_simple$Diuretics)
data_simple$Beta_Blockers = as.factor(data_simple$Beta_Blockers)
data_simple$Calcium_Blockers = as.factor(data_simple$Calcium_Blockers)
data_simple$ACE/ARB = as.factor(data_simple$ACE/ARB)
data_simple$Hyperten_Meds.2 = as.factor(data_simple$Hyperten_Meds.2)
data_simple$Sitting_Height = as.factor(data_simple$Sitting_Height)
data_simple$Height_Y/N = as.factor(data_simple$Height_Y/N)
data_simple$Waist_Y/N = as.factor(data_simple$Waist_Y/N)
data_simple$Hip_Y/N = as.factor(data_simple$Hip_Y/N)
data_simple$CVD_Family = as.factor(data_simple$CVD_Family)


##################################################################################################
#I have all my missing nodes eliminated and my variables assigned to the right class which 
#means I can start looking at how my variables are distributed. First I will start with 
#QQplots and histograms of the proxy variables 
age_histrogram = hist(data_simple$Age, main="Age distribution", xlab="Age", ylab="Frequency")
qqnorm(data_simple$Age)
qqline(data_simple$Age)

Trunk_histrogram = hist(data_simple$Trunk_Length, main="Trunk Length distribution", 
     xlab="Trunk Length", ylab="Frequency", breaks = 20)
qqnorm(data_simple$Trunk_Length)
qqline(data_simple$Trunk_Length)

Leg_length_histrogram = hist(data_simple$Leg_Length, main="Leg Length distribution", 
    xlab="Leg Length", ylab="Frequency", breaks = 20)
qqnorm(data_simple$Leg_Length)
qqline(data_simple$Leg_Length)

Sitting_height_histrogram = hist(data_simple$Sitting_Calculation, main="Sitting Height distribution", 
      xlab="Sitting Height", ylab="Frequency", breaks = 20)
qqnorm(data_simple$Sitting_Calculation)
qqline(data_simple$Sitting_Calculation)

WHR_histrogram = hist(data_simple$WHR, main="Waist to Hip Ratio distribution", 
    xlab="Waist to Hip Ratio", ylab="Frequency")
qqnorm(data_simple$WHR)
qqline(data_simple$WHR)


##################################################################################################
#I have a basic idea as to which variables are skewed and which variables are normally 
#distributed. I want to now look at the percentages of my variables which are not numeric
#but designated as characters or factors
summary(data_simple$Sex)
sex_freq <- table(data_simple$Sex)/708
barplot(sex_freq, main ="Percentage of males vs females in RODAM dataset", ylim= c(0, 1),
        col="lightblue")

summary(data_simple$Location)
location_freq <- table(data_simple$Location)/708
barplot(location_freq, main ="Location Distribution", ylim= c(0, 0.5), col="lavender")

summary(data_simple$CVD_Family)
CVD_freq <- table(data_simple$CVD_Family)/708
barplot(CVD_freq, main ="Family History of CVD", ylim= c(0, 1), col="plum1")


##################################################################################################
##Filtering out the outliers
data_simple <- data_simple %>% filter(data_simple$Leg_Length > 60)
Trunk_histrogram = hist(data_simple$Trunk_Length, main="Trunk Length distribution", 
    xlab="Trunk Length", ylab="Frequency", breaks = 10)
qqnorm(data_simple$Leg_Length)
qqline(data_simple$Leg_Length)

data_simple <- data_simple %>% filter(data_simple$Trunk_Length > 60)
Sitting_height_histrogram = hist(data_simple$Sitting_Calculation, main="Sitting Height distribution", 
    xlab="Sitting Height", ylab="Frequency", breaks = 10)
qqnorm(data_simple$Trunk_Length)
qqline(data_simple$Trunk_Length)

data_simple <- data_simple %>% filter(data_simple$WHR > 0.50)
WHR_histrogram = hist(data_simple$WHR, main="Waist to Hip Ratio distribution", 
    xlab="Waist to Hip Ratio", ylab="Frequency", breaks = 10)
qqnorm(data_simple$WHR)
qqline(data_simple$WHR)

data_simple <- data_simple[!(data_simple$RodamID=="G0696"),]
dim(data_simple)

###################################################################################################
#Checking the distribution between WHR and obesity as they are often related
ggplot(targets, aes(x=obesity, y=WHR, fill=obesity)) +
  geom_bar(stat="identity")+theme_minimal()

###################################################################################################

#write the file which will be used for your entire EWAS
write.csv(data_simple, "Rodam_Data_Simplified.csv")

##################################################################################################
## TABLE 1 SCRIPT

library(tableone)

table_factors <- c("Location", "Sex", "Age")

table_var <- c("Location", "Sex", "RiskCVD", "Smoking", "Hyperten_Meds.2", "Smoking", "DiabetesYN", "Heart_Condition", "Sitting_Calculation", "Leg_Length",
               "Trunk_Length", "WHR") 

table_location <- CreateTableOne(vars = table_var[-1],
      data = data_simple,
      factorVars = table_factors, 
      strata = 'Location')
table_sex <- CreateTableOne(vars = table_var[-1],
    data = data_simple,
    factorVars = table_factors, 
    strata = 'Sex')

locationT <- print(table_location, missing=T)
sexT <- print(table_sex, missing = T)

write.table(locationT, file = "location_table.txt", sep = ",", quote = FALSE, row.names = T)
write.table(sexT, file = "sex_table.txt", sep = ",", quote = FALSE, row.names = T)

