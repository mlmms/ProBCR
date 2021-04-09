#Normalization method comparison WITH BLOCK STRUCTURES

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
normalization = rep(c("Prostate-only", "Bounding-box"), each = reps)
model = rep(c("Lasso_I", "Lasso_II"), each = reps)


df <- data.frame(Fmax = c(Lasso_I[,"F_max"], Lasso_II[,"F_max"]))
df$Normalization = normalization
df$Model = model
df$Algorithm = rep("Lasso", each = reps*2)

ggplot(df, aes(x=Algorithm, y=Fmax, fill=Normalization)) + 
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


#### BLOCK STRUCTURE
# Group by Label
PL_1_pl1_II <- PL_1_II[PL_1_II$varNames == "pl1",]
PL_1_pl2_II <- PL_1_II[PL_1_II$varNames == "pl2",]


# create a data frame
reps = dim(PL_1_pl1_I)[1]
#bp = rep(c("Radiomics > Clinical", "Clinical > Radiomics", "Radiomics > Clinical", "Clinical > Radiomics"), each = reps)
bp = rep(c("bp1_I", "bp2_I", "bp1_II", "bp2_II"), each = reps)
normalization = rep(c("Prostate-only", "Bounding-box"), each = reps*2)
model = rep(c("PL_1_I", "PL_1_II"), each = reps*2)


df2 <- data.frame(Fmax = c(PL_1_pl1_I[,"max_F"], PL_1_pl2_I[,"max_F"], PL_1_pl1_II[,"max_F"], PL_1_pl2_II[,"max_F"]))
df2$Normalization = normalization
df2$Block_Priority = bp
df2$Model = model
df2$Algorithm = rep("Priority-Lasso-1 (2 blocks)", each = reps*2)

# dftotal = rbind(df,df2)
# 
# ggplot(dftotal, aes(x=Algorithm, y=Fmax, fill=Normalization)) + 
#   geom_boxplot()
# 
# ggplot(df2, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
#   geom_boxplot()


my_comparisons = list( c("bp1_I", "bp1_II"), c("bp2_I", "bp2_II"), c("bp1_I", "bp2_I"), c("bp1_II", "bp2_II") )

dftotal = rbind(df,df2,df3)

# CHOSEN ONE 

ggplot(df2, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons, annotations=c("ns", "***", "***", "***")) +
  labs(fill = "Normalisation \n method:") + xlab("Block priority sequence") + ylab(expression(F[max])) +
  ggtitle("Priority-Lasso-1 (2 blocks) performance by priority sequence") +
  theme(plot.title = element_text(size = 9), legend.title = element_text(size = 9)) +
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(legend.position="right")
  # annotate("text", x = 1, y = 0.57, label = "ns") +
  # annotate("text", x = 2, y = 0.85, label = "***") +
  # annotate("text", x = 2, y = 0.93, label = "***") +
  # annotate("text", x = 2, y = 1, label = "***") +
  # geom_segment(aes(x = 0.8, y = 0.54, xend = 1.2, yend = 0.54), size = 0.1) +
  # geom_segment(aes(x = 1.8, y = 0.81, xend = 2.2, yend = 0.81), size = 0.1) +
  # geom_segment(aes(x = 1.3, y = 0.93, xend = 1.8, yend = 0.93), size = 0.1) +
  # geom_segment(aes(x = 1.8, y = 1, xend = 2.2, yend = 1), size = 0.1) +

#Alternative chosen one
my_comparisons2 = list( c("bp1_I", "bp1_II"), c("bp2_I", "bp2_II") )

dftotal = rbind(df,df2,df3)

# CHOSEN ONE 
df2$Normalization <- factor(df2$Normalization,levels = c("Prostate-only", "Bounding-box"))

pl2 <- ggplot(df2, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
        geom_boxplot() +
        theme_classic() +
        stat_compare_means(comparisons = my_comparisons2, method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif", text.size = 20 ) +
        labs(fill = "Normalisation method:") + xlab("Block priority sequence") + ylab(expression(F[max])) +
        ggtitle("Priority-Lasso-1 (2 blocks) performance by block priority sequence") +
        #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
        theme(plot.title = element_text(size = 9, hjust = 0.5), legend.title = element_text(size = 8.5), legend.text = element_text(size = 8)) +
        theme(legend.position="right")

pl2$layers[[2]]$aes_params$textsize <- 3
pl2
pl2 + theme(axis.ticks = element_blank()) + scale_x_discrete(labels=c("bp1_I" = "Radiomics >", "bp1_II" = "Clinical", "bp2_I" = "Clinical >", "bp2_II" = "Radiomics"))


# ggplot(df3, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
#   geom_boxplot() +
#   stat_compare_means(method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif" , label.y = c(rep(0.7, times = 18), rep(0.85, times = 6))) +
#   labs(fill = "Normalisation method:") + xlab("Priority Sequence") + ylab(expression(F[max])) +
#   ggtitle("Priority-Lasso-2 (4 blocks) performances by priority sequence") +
#   theme(plot.title = element_text(size = 10), legend.title = element_text(size = 9)) +
#   #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
#   theme_classic() +
#   theme(legend.position="top")

# ggplot(dftotal, aes(x=Algorithm, y=Fmax, fill=Normalization)) + 
#   geom_boxplot() +
#   facet_wrap(~Normalization)

# ALTERNATIVE
ggplot(df2, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
    geom_boxplot() +
    #facet_wrap(~Model, scale = "free") +
    stat_compare_means(comparisons = my_comparisons, label.y = c(0.77, 0.87, 0.93, 1), method = "wilcox.test")


#Getting Wilcoxon statistics

wilcox.test(PL_1_pl1_I$max_F, PL_1_pl1_II$max_F, paired = TRUE)
wilcox.test(PL_1_pl2_I$max_F, PL_1_pl2_II$max_F, paired = TRUE)


wilcox.test(PL_1_pl1_I$max_F, PL_1_pl2_I$max_F, paired = TRUE) #

wilcox.test(PL_1_pl1_II$max_F, PL_1_pl2_II$max_F, paired = TRUE)


#other metrics, normalisation I
wilcox.test(PL_1_pl1_I$max_F, PL_1_pl2_I$max_F, paired = TRUE)
wilcox.test(PL_1_pl1_I$max_P, PL_1_pl2_I$max_P, paired = TRUE)
wilcox.test(PL_1_pl1_I$max_R, PL_1_pl2_I$max_R, paired = TRUE)

wilcox.test(PL_1_pl1_I$max_F0, PL_1_pl2_I$max_F0, paired = TRUE)
wilcox.test(PL_1_pl1_I$max_R0, PL_1_pl2_I$max_R0, paired = TRUE)
wilcox.test(PL_1_pl1_I$max_P0, PL_1_pl2_I$max_P0, paired = TRUE)

#other metrics, normalisation II

wilcox.test(PL_1_pl1_II$max_F, PL_1_pl2_II$max_F, paired = TRUE)
wilcox.test(PL_1_pl1_II$max_P, PL_1_pl2_II$max_P, paired = TRUE)
wilcox.test(PL_1_pl1_II$max_R, PL_1_pl2_II$max_R, paired = TRUE)

wilcox.test(PL_1_pl1_II$max_F0, PL_1_pl2_II$max_F0, paired = TRUE)
wilcox.test(PL_1_pl1_II$max_R0, PL_1_pl2_II$max_R0, paired = TRUE)
wilcox.test(PL_1_pl1_II$max_P0, PL_1_pl2_II$max_P0, paired = TRUE)


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
# 
# summary_stats_PL_2_I <- tableby(varNames ~ max_F, data = PL_2_I, control = mycontrols)
# summary_stats_PL_2_I
# summary(summary_stats_PL_2_I, text = TRUE, title = 'PLII w/o Oversampling, Norm I, dataI', pfootnote = TRUE)

# Group by Label
PL_2_pl1_I <- PL_2_I[PL_2_I$varNames == "pl1",]
PL_2_pl2_I <- PL_2_I[PL_2_I$varNames == "pl2",]
PL_2_pl3_I <- PL_2_I[PL_2_I$varNames == "pl3",]
PL_2_pl4_I <- PL_2_I[PL_2_I$varNames == "pl4",]
PL_2_pl5_I <- PL_2_I[PL_2_I$varNames == "pl5",]
PL_2_pl6_I <- PL_2_I[PL_2_I$varNames == "pl6",]
PL_2_pl7_I <- PL_2_I[PL_2_I$varNames == "pl7",]
PL_2_pl8_I <- PL_2_I[PL_2_I$varNames == "pl8",]
PL_2_pl9_I <- PL_2_I[PL_2_I$varNames == "pl9",]
PL_2_pl10_I <- PL_2_I[PL_2_I$varNames == "pl10",]
PL_2_pl11_I <- PL_2_I[PL_2_I$varNames == "pl11",]
PL_2_pl12_I <- PL_2_I[PL_2_I$varNames == "pl12",]
PL_2_pl13_I <- PL_2_I[PL_2_I$varNames == "pl13",]
PL_2_pl14_I <- PL_2_I[PL_2_I$varNames == "pl14",]
PL_2_pl15_I <- PL_2_I[PL_2_I$varNames == "pl15",]
PL_2_pl16_I <- PL_2_I[PL_2_I$varNames == "pl16",]
PL_2_pl17_I <- PL_2_I[PL_2_I$varNames == "pl17",]
PL_2_pl18_I <- PL_2_I[PL_2_I$varNames == "pl18",]
PL_2_pl19_I <- PL_2_I[PL_2_I$varNames == "pl19",]
PL_2_pl20_I <- PL_2_I[PL_2_I$varNames == "pl20",]
PL_2_pl21_I <- PL_2_I[PL_2_I$varNames == "pl21",]
PL_2_pl22_I <- PL_2_I[PL_2_I$varNames == "pl22",]
PL_2_pl23_I <- PL_2_I[PL_2_I$varNames == "pl23",]
PL_2_pl24_I <- PL_2_I[PL_2_I$varNames == "pl24",]

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

summary_stats_PL_2_II <- tableby(varNames ~ max_F, data = PL_2_II, control = mycontrols)
summary_stats_PL_2_II
summary(summary_stats_PL_2_II, text = TRUE, title = 'PLII w/o Oversampling, Norm I, PL_2_II', pfootnote = TRUE)

summary_stats_PL_2_II <- tableby(varNames ~ ., data = PL_2_II, control = mycontrols)
summary_stats_PL_2_II
summary(summary_stats_PL_2_II, text = TRUE, title = 'PLII w/o Oversampling, Norm I, PL_2_II', pfootnote = TRUE)


# NEW IN

#### BLOCK STRUCTURE
# Group by Label
PL_2_pl1_II <- PL_2_II[PL_2_II$varNames == "pl1",]
PL_2_pl2_II <- PL_2_II[PL_2_II$varNames == "pl2",]
PL_2_pl3_II <- PL_2_II[PL_2_II$varNames == "pl3",]
PL_2_pl4_II <- PL_2_II[PL_2_II$varNames == "pl4",]
PL_2_pl5_II <- PL_2_II[PL_2_II$varNames == "pl5",]
PL_2_pl6_II <- PL_2_II[PL_2_II$varNames == "pl6",]
PL_2_pl7_II <- PL_2_II[PL_2_II$varNames == "pl7",]
PL_2_pl8_II <- PL_2_II[PL_2_II$varNames == "pl8",]
PL_2_pl9_II <- PL_2_II[PL_2_II$varNames == "pl9",]
PL_2_pl10_II <- PL_2_II[PL_2_II$varNames == "pl10",]
PL_2_pl11_II <- PL_2_II[PL_2_II$varNames == "pl11",]
PL_2_pl12_II <- PL_2_II[PL_2_II$varNames == "pl12",]
PL_2_pl13_II <- PL_2_II[PL_2_II$varNames == "pl13",]
PL_2_pl14_II <- PL_2_II[PL_2_II$varNames == "pl14",]
PL_2_pl15_II <- PL_2_II[PL_2_II$varNames == "pl15",]
PL_2_pl16_II <- PL_2_II[PL_2_II$varNames == "pl16",]
PL_2_pl17_II <- PL_2_II[PL_2_II$varNames == "pl17",]
PL_2_pl18_II <- PL_2_II[PL_2_II$varNames == "pl18",]
PL_2_pl19_II <- PL_2_II[PL_2_II$varNames == "pl19",]
PL_2_pl20_II <- PL_2_II[PL_2_II$varNames == "pl20",]
PL_2_pl21_II <- PL_2_II[PL_2_II$varNames == "pl21",]
PL_2_pl22_II <- PL_2_II[PL_2_II$varNames == "pl22",]
PL_2_pl23_II <- PL_2_II[PL_2_II$varNames == "pl23",]
PL_2_pl24_II <- PL_2_II[PL_2_II$varNames == "pl24",]


# create a data frame
reps = dim(PL_2_I)[1]
normalization = rep(c("Prostate-only", "Bounding-box"), each = reps)
model = rep(c("PL_2_I", "PL_2_II"), each = reps)


df3 <- data.frame(Fmax = c(PL_2_I[,"max_F"], PL_2_II[,"max_F"]))
bpI <- rep(1:24, times = dim(PL_2_I)[1]/24)
bpII <- rep(1:24, times = dim(PL_2_I)[1]/24)
df3$Block_Priority <- factor(c(bpI, bpII))
df3$Normalization = normalization
df3$Model = model
df3$Algorithm = rep("Priority-Lasso 2", each = reps*2)
df3$Highlighted = rep(highlighted)
df3$bp_type <- c(rep(paste0("bp_I_",1:24), times = dim(PL_2_I)[1]/24), rep(paste0("bp_II_",1:24), times = dim(PL_2_I)[1]/24))


# TRYING
# my_comparisons = list( c("bp_I_1", "bp_II_1"), 
#                        c("bp_I_2", "bp_II_2"), 
#                        c("bp_I_3", "bp_II_3"),
#                        c("bp_I_4", "bp_II_4"),
#                        c("bp_I_5", "bp_II_5"),
#                        c("bp_I_6", "bp_II_6"),
#                        c("bp_I_7", "bp_II_7"),
#                        c("bp_I_8", "bp_II_8"),
#                        c("bp_I_9", "bp_II_9"),
#                        c("bp_I_10", "bp_II_10"),
#                        c("bp_I_11", "bp_II_11"),
#                        c("bp_I_12", "bp_II_12"),
#                        c("bp_I_14", "bp_II_13"))

# ggplot(df3, aes(x=bp_type, y=Fmax, fill=Normalization)) + 
#   geom_boxplot() + 
#   #geom_signif(comparisons = my_comparisons, annotations=c(rep(1:13)), y_position = c(0.93), tip_length = 0.01)
#   stat_compare_means(comparisons = my_comparisons, label.y = c(0.9), method = "wilcox.test")


#Alternative
ggplot(df3, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  facet_wrap(~Normalization)


# chosen: OVERALL PL2_I VS. PL2_II, with p value significances: 
df3$Normalization <- factor(df3$Normalization,levels = c("Prostate-only", "Bounding-box"))

pl3 <- ggplot(df3, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif" , label.y = c(rep(0.7, times = 18), rep(0.85, times = 6))) +
  labs(fill = "Normalisation method:") + xlab("Priority Sequence") + ylab(expression(F[max])) +
  ggtitle("Priority-Lasso-2 (4 blocks) performances by priority sequence") +
  theme(plot.title = element_text(size = 10), legend.title = element_text(size = 9)) +
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(legend.position="top", plot.title = element_text(size = 9, hjust = 0.5), legend.title = element_text(size = 8.5), legend.text = element_text(size = 8))

  
  pl3$layers[[2]]$aes_params$size <- 3
  pl3
  pl3 + theme(axis.ticks = element_blank()) 
  #+ scale_x_discrete(labels=c("bp1_I" = "Radiomics >", "bp1_II" = "Clinical", "bp2_I" = "Clinical >", "bp2_II" = "Radiomics"))
  

# chosen: OVERALL PL2_I VS. PL2_II, with p value significances, WILCOXON test method with alternative = greater
ggplot(df3, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", paired = TRUE, method.args = list(alternative = "less"), size = 4, label = "p.signif" , label.y = c(rep(0.7, times = 18), rep(0.85, times = 6))) +
  labs(fill = "Normalisation method:") + xlab("Priority Sequence") + ylab(expression(F[max])) +
  ggtitle("Priority-Lasso-2 (4 blocks) performances by priority sequence") +
  theme(plot.title = element_text(size = 10), legend.title = element_text(size = 9)) +
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(legend.position="top")


# BLOCK PER BLOCK NORMALISATION COMPARISON
symnumargs2 <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("hs", "****", "***", "**", "*", "ns"))

ggplot(df3, aes(x=Algorithm, y=Fmax, fill=Normalization)) +
  geom_boxplot() +
  facet_wrap(~Block_Priority) + 
  stat_compare_means(method = "wilcox.test", size = 3, label.y = c(0.8) )


# ALTERNATIVE
ggplot(df3, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  #facet_wrap(~Model, scale = "free") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.9), method = "wilcox.test")


# TRYING HIGHLIGHTING PLOT

# Add a column called 'type': do we want to highlight the group or not?"
df4 <- df3[df3$Normalization == "Prostate-Norm",]

highlighted <- as.logical(df4$Block_Priority == "19" | df4$Block_Priority == "20")

ggplot(df4, aes(x=Block_Priority, y=Fmax, fill=highlighted)) + 
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette="Dark2")
  #xlab("")




# Comparing PL19_I and II and PL20_I and II
PL_2_pl19_I <- PL_2_I[PL_2_I$varNames == "pl19",]
PL_2_pl20_I <- PL_2_I[PL_2_I$varNames == "pl20",]

PL_2_pl19_II <- PL_2_II[PL_2_II$varNames == "pl19",]
PL_2_pl20_II <- PL_2_II[PL_2_II$varNames == "pl20",]

# Statistics
wilcox.test(PL_2_pl19_I$max_F, PL_2_pl19_II$max_F, paired = TRUE)

wilcox.test(PL_2_pl19_I$max_F, PL_2_pl20_I$max_F, paired = TRUE)
wilcox.test(PL_2_pl19_I$max_F, PL_2_pl20_II$max_F, paired = TRUE)


wilcox.test(PL_2_pl20_I$max_F, PL_2_pl20_II$max_F, paired = TRUE)

# create a data frame
reps = dim(PL_2_pl19_I)[1]
bp = rep(c("bp19_I", "bp20_I", "bp19_II", "bp20_II"), each = reps)
normalization = rep(c("Prostate-Norm", "Whole-image Norm"), each = reps*2)
model = rep(c("PL_2_I", "PL_2_II"), each = reps*2)


df5 <- data.frame(Fmax = c(PL_2_pl19_I[,"max_F"], PL_2_pl20_I[,"max_F"], PL_2_pl19_II[,"max_F"], PL_2_pl20_II[,"max_F"]))
df5$Normalization = normalization
df5$Block_Priority = bp
df5$Model = model
df5$Algorithm = rep("Priority-Lasso 2 (4 blocks)", each = reps*2)


my_comparisons = list( c("bp19_I", "bp19_II"), c("bp20_I", "bp20_II"), c("bp19_I", "bp20_I"), c("bp19_II", "bp20_II")  )

#dftotal = rbind(df,df2,df3)

# CHOSEN ONE 

#changing x order
#df4$Block_Priority2 <- factor(df4$Block_Priority, levels = c("bp19_I", "bp20_I", "bp19_II", "bp20_II"))

ggplot(df5, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons, annotations=c("***", "****", "*", "***"), y_position = c(0.89, 0.89, 0.94, 1), tip_length = 0.01)

# ggplot(df5, aes(x=Block_Priority, y=Fmax, fill=Normalization)) +
#   geom_boxplot() +
#   facet_wrap(~Normalization)

# ALTERNATIVE
ggplot(df5, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  #facet_wrap(~Model, scale = "free") +
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.85, 0.85, 0.93, 1), method = "wilcox.test", paired = TRUE)


#Chosen one:
ggplot(df5, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  #facet_wrap(~Model, scale = "free") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = TRUE, label = "p.signif", size = 4) +
  labs(fill = "Normalisation \n method:") + xlab("Priority Sequence") + ylab(expression(F[max])) +
  ggtitle("Priority-Lasso-2 (4 blocks): Best performing priority sequences") +
  theme(plot.title = element_text(size = 9), legend.title = element_text(size = 9)) +
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(legend.position="right")

#Alternative graph

my_comparisons2 = list( c("bp19_I", "bp19_II"), c("bp20_I", "bp20_II")  )
#p.format
ggplot(df5, aes(x=Block_Priority, y=Fmax, fill=Normalization)) + 
  geom_boxplot() +
  #facet_wrap(~Model, scale = "free") +
  stat_compare_means(comparisons = my_comparisons2, method = "wilcox.test", paired = TRUE, label = "p.format", size = 4) +
  labs(fill = "Normalisation \n method:") + xlab("Priority Sequence") + ylab(expression(F[max])) +
  ggtitle("Priority-Lasso-2 (4 blocks): Best performing priority sequences") +
  theme(plot.title = element_text(size = 9), legend.title = element_text(size = 9)) +
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(legend.position="right")