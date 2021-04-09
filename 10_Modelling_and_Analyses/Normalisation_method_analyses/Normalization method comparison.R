#Normalisation method comparison WITHOUT BLOCK STRUCTURES
 
library(openxlsx)
library(arsenal)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)



prebio_files <- list.files(path = "/Users/monicasilva/Documents/ANALYSIS/1_PREBIOPSY", 
                           pattern = NULL, full.names = TRUE, ignore.case = FALSE)
#prebio_files <- as.data.frame(prebio_files)
mycontrols <- tableby.control(test =TRUE, total = FALSE)


# File Loading ------------------------------------------------------------
Lasso_NO_Norm_I_Prebio <- read.xlsx(prebio_files[1])
Lasso_NO_Norm_II_Prebio <- read.xlsx(prebio_files[2])
Lasso_SMOTE_Norm_I_Prebio <- read.xlsx(prebio_files[3])
Lasso_SMOTE_Norm_II_Prebio <- read.xlsx(prebio_files[4])
PLI_NO_Norm_I_Prebio <- read.xlsx(prebio_files[5])
PLI_NO_Norm_II_Prebio <- read.xlsx(prebio_files[6])
PLI_SMOTE_Norm_I_Prebio <- read.xlsx(prebio_files[7])
PLI_SMOTE_Norm_II_Prebio <- read.xlsx(prebio_files[8])
PLII_NO_Norm_I_Prebio <- read.xlsx(prebio_files[9])
PLII_NO_Norm_II_Prebio <- read.xlsx(prebio_files[10])
PLII_SMOTE_Norm_I_Prebio <- read.xlsx(prebio_files[11])
PLII_SMOTE_Norm_II_Prebio <- read.xlsx(prebio_files[12])


# Normalisation method comparison per algorithm

# Lasso I

Lasso_I <- Lasso_NO_Norm_I_Prebio
# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
Lasso_I[is.na(Lasso_I)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
Lasso_I <- Lasso_I[Lasso_I$nzcoef <= 10,] # COUNTING WITH INTERCEPT
Lasso_I <- Lasso_I[(Lasso_I$pmaxs_value != "NL"),]

#Encode and label
#(dim(data_LASSO)[2])
cols_to_factor <- colnames(Lasso_I[,c(1,2,3,15)])
cols_to_numeric <- colnames(Lasso_I[,-c(1,2,3,15)])
Lasso_I[,cols_to_factor] <- lapply(Lasso_I[,cols_to_factor],as.factor)
Lasso_I[,cols_to_numeric] <- lapply(Lasso_I[,cols_to_numeric],as.numeric)

# Group by Label

Lasso_I <- Lasso_I[-c(2,3,12:14,16)]

# summary_stats <- tableby(pmaxs_value ~ ., data = Lasso_I, control = mycontrols)
# summary_stats
# summary(summary_stats, text = TRUE, title = 'LASSO w/o Oversampling, Norm I', pfootnote = TRUE)


# Lasso II

Lasso_II <- Lasso_NO_Norm_II_Prebio
# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
Lasso_II[is.na(Lasso_II)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
Lasso_II <- Lasso_II[Lasso_II$nzcoef <= 10,] # COUNTING WITH INTERCEPT
Lasso_II <- Lasso_II[(Lasso_II$pmaxs_value != "NL"),]

#Encode and label
#(dim(data_LASSO)[2])
cols_to_factor <- colnames(Lasso_II[,c(1,2,3,15)])
cols_to_numeric <- colnames(Lasso_II[,-c(1,2,3,15)])
Lasso_II[,cols_to_factor] <- lapply(Lasso_II[,cols_to_factor],as.factor)
Lasso_II[,cols_to_numeric] <- lapply(Lasso_II[,cols_to_numeric],as.numeric)

# Group by Label

Lasso_II <- Lasso_II[-c(2,3,12:14,16)]

# summary_stats <- tableby(pmaxs_value ~ ., data = Lasso_II, control = mycontrols)
# summary_stats
# summary(summary_stats, text = TRUE, title = 'LASSO w/o Oversampling, Norm II', pfootnote = TRUE)



# create a data frame
reps = dim(Lasso_I)[1]
normalisation = rep(c("Prostate-only", "Bounding-box"), each = reps)
model = rep(c("Lasso_I", "Lasso_II"), each = reps)


df <- data.frame(Fmax = c(Lasso_I[,"F_max"], Lasso_II[,"F_max"]))
df$Normalisation = normalisation
df$Model = model
df$Algorithm = rep("Lasso", each = reps*2)

ggplot(df, aes(x=Algorithm, y=Fmax, fill=Normalisation)) + 
  geom_boxplot()

#######

# Priority-Lasso-I Norm I

PL_1_I <- PLI_NO_Norm_I_Prebio
PL_1_I <- PL_1_I[!( PL_1_I$pmaxs_value==5 | PL_1_I$pmaxs_value==6 | PL_1_I$pmaxs_value==7 | PL_1_I$pmaxs_value==8 | PL_1_I$pmaxs_value==9 ) ,]

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_1_I[is.na(PL_1_I)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_1_I <- PL_1_I[PL_1_I$total_coeffs <= 10,] # COUNTING WITH INTERCEPT
#PL_1_I <- PL_1_I[-123,]

#Encode and label
PL_1_I <- PL_1_I[,c(1,4,9,11:18)]
cols_to_factor <- colnames(PL_1_I[,c(1,2,3)])
cols_to_numeric <- colnames(PL_1_I[,-c(1,2,3)])
PL_1_I[,cols_to_factor] <- lapply(PL_1_I[,cols_to_factor],as.factor)
PL_1_I[,cols_to_numeric] <- lapply(PL_1_I[,cols_to_numeric],as.numeric)

# Group by Label

PL_1_pl1_I <- PL_1_I[PL_1_I$varNames == "pl1",]
PL_1_pl2_I <- PL_1_I[PL_1_I$varNames == "pl2",]


## 
#Priority-Lasso I Norm II

PL_1_II <- PLI_NO_Norm_II_Prebio
PL_1_II <- PL_1_II[!( PL_1_II$pmaxs_value==5 | PL_1_II$pmaxs_value==6 | PL_1_II$pmaxs_value==7 | PL_1_II$pmaxs_value==8 | PL_1_II$pmaxs_value==9 ) ,]


# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_1_II[is.na(PL_1_II)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_1_II <- PL_1_II[PL_1_II$total_coeffs <= 10,] # COUNTING WITH INTERCEPT

#Encode and label
PL_1_II <- PL_1_II[,c(1,4,9,11:18)]
cols_to_factor <- colnames(PL_1_II[,c(1,2,3)])
cols_to_numeric <- colnames(PL_1_II[,-c(1,2,3)])
PL_1_II[,cols_to_factor] <- lapply(PL_1_II[,cols_to_factor],as.factor)
PL_1_II[,cols_to_numeric] <- lapply(PL_1_II[,cols_to_numeric],as.numeric)


# Group by Label
# PL_1_pl1_II <- PL_1_II[PL_1_II$varNames == "pl1",]
# PL_1_pl2_II <- PL_1_II[PL_1_II$varNames == "pl2",]


# create a data frame
reps = dim(PL_1_I)[1]
normalisation = rep(c("Prostate-only", "Bounding-box"), each = reps)
model = rep(c("PL_1_I", "PL_1_II"), each = reps)


df2 <- data.frame(Fmax = c(PL_1_I[,"max_F"], PL_1_II[,"max_F"]))
df2$Normalisation = normalisation
df2$Model = model
df2$Algorithm = rep("Priority-Lasso-1 (2 blocks)", each = reps*2)


dftotal = rbind(df,df2)
ggplot(dftotal, aes(x=Algorithm, y=Fmax, fill=Normalisation)) + 
  geom_boxplot()


# Priority-LASSO_II -------------------------------------------------------

# Model selection ---------------------------------------------------------
# Choosing between configurations for one algorithm
PL_2_I <- PLII_NO_Norm_I_Prebio
PL_2_I <- PL_2_I[!( PL_2_I$pmaxs_value==4 | PL_2_I$pmaxs_value==5 |
                      PL_2_I$pmaxs_value==6 |  PL_2_I$pmaxs_value== 7 | PL_2_I$pmaxs_value==8 | 
                      PL_2_I$pmaxs_value== 9 | PL_2_I$pmaxs_value== "1,4,4,4" | PL_2_I$pmaxs_value==3 |
                      PL_2_I$pmaxs_value== "3,3,3,3" |
                      PL_2_I$pmaxs_value== "3,3,3,2" | PL_2_I$pmaxs_value== "3,3,2,2" | PL_2_I$pmaxs_value== "2,4,4,0" | 
                      PL_2_I$pmaxs_value== "4,4,2,0" | PL_2_I$pmaxs_value== "6,4,0,0" | PL_2_I$pmaxs_value== "1,2,3,4" | PL_2_I$pmaxs_value== "1,1,4,4" |
                      PL_2_I$pmaxs_value== "1,3,3,3"), ]  # 

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_2_I[is.na(PL_2_I)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_2_I <- PL_2_I[PL_2_I$total_coeffs < 10,] # COUNTING WITH INTERCEPT

#Encode and label
PL_2_I <- PL_2_I[,c(1,4,13,15:22)]
cols_to_factor <- colnames(PL_2_I[,c(1,2,3)])
cols_to_numeric <- colnames(PL_2_I[,-c(1,2,3)])
PL_2_I[,cols_to_factor] <- lapply(PL_2_I[,cols_to_factor],as.factor)
PL_2_I[,cols_to_numeric] <- lapply(PL_2_I[,cols_to_numeric],as.numeric)


# summary_stats_PL_2_I <- tableby(varNames ~ ., data = PL_2_I, control = mycontrols)
# summary_stats_PL_2_I
# summary(summary_stats_PL_2_I, text = TRUE, title = 'PLII w/o Oversampling, Norm I, dataI', pfootnote = TRUE)


#Priority Lasso 2 _II

PL_2_II <- PLII_NO_Norm_II_Prebio
PL_2_II <- PL_2_II[!( PL_2_II$pmaxs_value==4 | PL_2_II$pmaxs_value==5 |
                        PL_2_II$pmaxs_value==6 |  PL_2_II$pmaxs_value== 7 | PL_2_II$pmaxs_value==8 | 
                        PL_2_II$pmaxs_value== 9 | PL_2_II$pmaxs_value== "1,4,4,4" | PL_2_II$pmaxs_value==3 |
                        PL_2_II$pmaxs_value== "3,3,3,3" |
                        PL_2_II$pmaxs_value== "3,3,3,2" | PL_2_II$pmaxs_value== "3,3,2,2" | PL_2_II$pmaxs_value== "2,4,4,0" | 
                        PL_2_II$pmaxs_value== "4,4,2,0" | PL_2_II$pmaxs_value== "6,4,0,0" | PL_2_II$pmaxs_value== "1,2,3,4" | PL_2_II$pmaxs_value== "1,1,4,4" |
                        PL_2_II$pmaxs_value== "1,3,3,3"), ]  # 


# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_2_II[is.na(PL_2_II)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_2_II <- PL_2_II[PL_2_II$total_coeffs < 10,]

#Encode and label
PL_2_II <- PL_2_II[,c(1,4,13,15:22)]
cols_to_factor <- colnames(PL_2_II[,c(1,2,3)])
cols_to_numeric <- colnames(PL_2_II[,-c(1,2,3)])
PL_2_II[,cols_to_factor] <- lapply(PL_2_II[,cols_to_factor],as.factor)
PL_2_II[,cols_to_numeric] <- lapply(PL_2_II[,cols_to_numeric],as.numeric)

# summary_stats_PL_2_II <- tableby(varNames ~ ., data = PL_2_II, control = mycontrols)
# summary_stats_PL_2_II
# summary(summary_stats_PL_2_II, text = TRUE, title = 'PLII w/o Oversampling, Norm I, PL_2_II', pfootnote = TRUE)


# create a data frame
reps = dim(PL_2_I)[1]
normalisation = rep(c("Prostate-only", "Bounding-box"), each = reps)
model = rep(c("PL_2_I", "PL_2_II"), each = reps)


df3 <- data.frame(Fmax = c(PL_2_I[,"max_F"], PL_2_II[,"max_F"]))
df3$Normalisation = normalisation
df3$Model = model
df3$Algorithm = rep("Priority-Lasso-2 (4 blocks)", each = reps*2)


my_comparisons = list( c("Lasso_I", "Lasso_II"), c("PL_1_I", "PL_1_II"), c("PL_2_I", "PL_2_II") )

dftotal = rbind(df,df2,df3)

#Just to see the p-values
ggplot(dftotal, aes(x=Model, y=Fmax, fill=Normalisation)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.9), method = "wilcox.test")

#FINAL GRAPH

my_comparisons = list( c("Lasso_I", "Lasso_II"), c("PL_1_I", "PL_1_II"), c("PL_2_I", "PL_2_II") )
dftotal$Normalisation <- factor(dftotal$Normalisation,levels = c("Prostate-only", "Bounding-box"))

ggplot(dftotal, aes(x=Algorithm, y=Fmax, fill=Normalisation)) + 
  geom_boxplot() + 
  geom_signif(comparisons = my_comparisons, annotations=c("****", "n.s.", "n.s."), y_position = c(0.9), tip_length = 0.01) +
  labs(fill = "Normalisation method:") + ylab(expression(F[max])) +
  ggtitle("Overview of algorithms' performances per normalisation method") +
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(plot.title = element_text(size = 9, hjust = 0.5), legend.title = element_text(size = 8.5), legend.text = element_text(size = 8)) +
  theme(legend.position="right") +
  annotate("text", x = 1, y = 0.9, label = "****", size = 2) +
  annotate("text", x = 2, y = 0.91, label = "ns", size = 2) +
  annotate("text", x = 3, y = 0.91, label = "ns", size = 2) +
  geom_segment(aes(x = 0.8, y = 0.88, xend = 1.2, yend = 0.88), size = 0.1) +
  geom_segment(aes(x = 1.8, y = 0.88, xend = 2.2, yend = 0.88), size = 0.1) +
  geom_segment(aes(x = 2.8, y = 0.88, xend = 3.2, yend = 0.88), size = 0.1) +
  theme(axis.ticks = element_blank()) 



#Wilcoxon tests

t <- dftotal[dftotal[,"Model"] == "Lasso_I",]
t2 <- dftotal[dftotal[,"Model"] == "Lasso_II",]

wilcox.test(t$Fmax, t2$Fmax, paired = TRUE, alternative = "two.sided")

# ggplot(dftotal, aes(x=Algorithm, y=Fmax, fill=Normalisation)) + 
#   geom_boxplot() +
#   facet_wrap(~Normalisation)

# ggplot(dftotal, aes(x=Algorithm, y=Fmax, fill=Normalisation)) + 
#   geom_boxplot() +
#   facet_wrap(~Algorithm, scale = "free") +
#   theme(legend.position="none")
