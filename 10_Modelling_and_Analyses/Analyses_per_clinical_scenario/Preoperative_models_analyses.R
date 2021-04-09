#Choosing best hyperparameters - Pre-Biopsy

library(openxlsx)
library(arsenal)
library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)
library(PMCMRplus)
library(PMCMR)
library(ggtext)

pretreatment_files <- list.files(path = "/Users/monicasilva/Documents/ANALYSIS/2_PREOPERATIVE", 
                           pattern = NULL, full.names = TRUE, ignore.case = FALSE)
#prebio_files <- as.data.frame(pretreatment_files)
mycontrols <- tableby.control(test =TRUE, total = FALSE)


# File Loading ------------------------------------------------------------
Lasso_NO_Norm_I_Pretreat <- read.xlsx(pretreatment_files[1])
#Lasso_SMOTE_Norm_I_Pretreat <- read.xlsx(pretreatment_files[3])
PLI_NO_Norm_I_Pretreat <- read.xlsx(pretreatment_files[5])
#PLI_SMOTE_Norm_I_Pretreat <- read.xlsx(pretreatment_files[7])
PLII_NO_Norm_I_Pretreat <- read.xlsx(pretreatment_files[9])
#PLII_SMOTE_Norm_I_Pretreat <- read.xlsx(pretreatment_files[11])


# PART I ------------------------------------------------------------------

# Model selection ---------------------------------------------------------
# Choosing between configurations for one algorithm
data_I <- Lasso_NO_Norm_I_Pretreat
# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
data_I[is.na(data_I)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
data_I <- data_I[data_I$nzcoef <= 10,] # COUNTING WITH INTERCEPT
data_I <- data_I[(data_I$pmaxs_value != "NL"),]

# Friedman test
# Performing statistical tests
LASSO_pmax5 <- data_I[data_I$pmaxs_value == 5,]
LASSO_pmax6 <- data_I[data_I$pmaxs_value == 6,]
LASSO_pmax7 <- data_I[data_I$pmaxs_value == 7,]
LASSO_pmax8 <- data_I[data_I$pmaxs_value == 8,]
LASSO_pmax9 <- data_I[data_I$pmaxs_value == 9,]

LASSO_Fmax5 <- LASSO_pmax5$F_max
LASSO_Fmax6 <- LASSO_pmax6$F_max
LASSO_Fmax7 <- LASSO_pmax7$F_max
LASSO_Fmax8 <- LASSO_pmax8$F_max
LASSO_Fmax9 <- LASSO_pmax9$F_max

df <- cbind(LASSO_Fmax5, LASSO_Fmax6, LASSO_Fmax7, LASSO_Fmax8, LASSO_Fmax9)

#Friedman test
friedman.test(df)
posthoc.friedman.nemenyi.test(df)

#Post-Hoc Test, posthoc.friedman.nemenyi.test
posthoc.friedman.nemenyi.test(y=df)

#Plotting
# create a data frame
reps = dim(df)[1]
pmaxs = rep(5:9, each = reps)

df4 <- data_I[,c(1,9)]
row_sub = apply(df4, 1, function(row) all(row !=0 ))

comparisons = list( c("7","8"), c("7","9"), c("6","9"), c("5","9") )

#PLOT

pl9 <- ggplot(df4, aes(x=pmaxs_value, y=F_max, fill=pmaxs_value)) + 
          geom_boxplot() +
          scale_fill_brewer(palette = "Greens") +
          #stat_compare_means(comparisons = comparisons, method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif" ) + #, label.y = c(rep(0.7, times = 18), rep(0.85, times = 6)
          labs(fill = "") + 
          xlab("Maximum number of model variables (pmax)") + 
          ylab(expression(F[max])) +
          ggtitle("Pre-Operative Lasso-penalised models' performance \n by maximum amount of features allowed") +
          theme(plot.title = element_text(size = 10), legend.title = element_text(size = 9)) +
          #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
          theme_classic() + #geom_jitter(color="black", size=0.4, alpha=0.9) +
          theme(legend.position="none", plot.title = element_text(size = 9, hjust = 0.5), axis.title = element_text(size = 9)) +
          theme(axis.ticks = element_blank())
          theme(axis.title.x = element_markdown())

pl9$layers[[2]]$aes_params$textsize <- 3
pl9

#Descriptitive statistics
#Encode and label
#(dim(data_LASSO)[2])
cols_to_factor <- colnames(data_I[,c(1,2,3,15)])
cols_to_numeric <- colnames(data_I[,-c(1,2,3,15)])
data_I[,cols_to_factor] <- lapply(data_I[,cols_to_factor],as.factor)
data_I[,cols_to_numeric] <- lapply(data_I[,cols_to_numeric],as.numeric)

# Group by Label
data_I <- data_I[-c(2,3,12:14,16)]
mycontrols <- tableby.control(test =TRUE, total = FALSE)

summary_stats <- tableby(pmaxs_value ~ ., data = data_I, control = mycontrols)
summary_stats
summary(summary_stats, text = TRUE, title = 'LASSO w/o Oversampling, Norm I', pfootnote = TRUE)

#write2word(summary_stats, "/Users/monicasilva/Documents/Thesis_Writing/RESULTS_raw/Lasso_I", pfootnote = TRUE)


# CONCLUSION: BEST MODEL

PT_Lasso_best_model <- LASSO_pmax9

# PRIORITY-LASSO I --------------------------------------------------------

# Model selection ---------------------------------------------------------
# Choosing between configurations for one algorithm
data_I <- PLI_NO_Norm_I_Pretreat
data_I <- data_I[!( data_I$pmaxs_value==5 | data_I$pmaxs_value==6 | data_I$pmaxs_value==7 | data_I$pmaxs_value==8 | data_I$pmaxs_value==9 ) ,]

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
data_I[is.na(data_I)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
data_I <- data_I[data_I$total_coeffs <= 10,] # COUNTING WITH INTERCEPT
#data_I <- data_I[-123,]

#Encode and label
data_I <- data_I[,c(1,4,9,11:18)]
cols_to_factor <- colnames(data_I[,c(1,2,3)])
cols_to_numeric <- colnames(data_I[,-c(1,2,3)])
data_I[,cols_to_factor] <- lapply(data_I[,cols_to_factor],as.factor)
data_I[,cols_to_numeric] <- lapply(data_I[,cols_to_numeric],as.numeric)

data_pl1_I <- data_I[data_I$varNames == "pl1",]
data_pl2_I <- data_I[data_I$varNames == "pl2",]

wilcox.test(data_pl1_I$max_F, data_pl2_I$max_F, paired = TRUE)
mycontrols <- tableby.control(test =TRUE, total = FALSE)

# summary_stats_pl1_I <- tableby(pmaxs_value ~ ., data = data_pl1_I, control = mycontrols)
# summary_stats_pl1_I
# summary(summary_stats_pl1_I, text = TRUE, title = 'PLI w/o Oversampling, Norm I, pl1', pfootnote = TRUE)

#write2word(summary_stats_pl1_I, "/Users/monicasilva/Documents/Thesis_Writing/RESULTS_raw/PL-I_pl1_I", pfootnote = TRUE)


summary_stats_pl2_I <- tableby(pmaxs_value ~ ., data = data_pl2_I, control = mycontrols)
summary_stats_pl2_I
summary(summary_stats_pl2_I, text = TRUE, title = 'PLI w/o Oversampling, Norm I, pl2', pfootnote = TRUE)

#write2word(summary_stats_pl2_I, "/Users/monicasilva/Documents/Thesis_Writing/RESULTS_raw/PL-I_pl2_I", pfootnote = TRUE)


# Statistical tests -------------------------------------------------------
#NEW
# shapiro.test(data_pl1_I$max_F)
# wilcox.test(data_pl1_I$max_F, data_pl2_I$max_F, paired = TRUE, alternative = "two.sided")

# Friedman test
# Performing statistical tests
PL_I_pmax15 <- data_pl2_I[data_pl2_I$pmaxs_value == "1,5",]
PL_I_pmax17 <- data_pl2_I[data_pl2_I$pmaxs_value == "1,7",]
PL_I_pmax18 <- data_pl2_I[data_pl2_I$pmaxs_value == "1,8",]
PL_I_pmax27 <- data_pl2_I[data_pl2_I$pmaxs_value == "2,7",]
PL_I_pmax36 <- data_pl2_I[data_pl2_I$pmaxs_value == "3,6",]
PL_I_pmax54 <- data_pl2_I[data_pl2_I$pmaxs_value == "5,4",]


PL_I_Fmax15 <- PL_I_pmax15$max_F
PL_I_Fmax17 <- PL_I_pmax17$max_F
PL_I_Fmax18 <- PL_I_pmax18$max_F
PL_I_Fmax27 <- PL_I_pmax27$max_F
PL_I_Fmax36 <- PL_I_pmax36$max_F
PL_I_Fmax54 <- PL_I_pmax54$max_F


df <- cbind(PL_I_Fmax15, PL_I_Fmax17, PL_I_Fmax18, PL_I_Fmax27, PL_I_Fmax36, PL_I_Fmax54)

friedman.test(df)

#Post-Hoc Test, posthoc.friedman.nemenyi.test
posthoc.friedman.nemenyi.test(df)


# create a data frame
df_PL_I <- data_pl2_I[,c(1,6)]

comparisons_PL_I_signif = list( c("1,8","5,4") )

#Como o friedman da' significativo, vamos fazer boxplot e aplicar wilcoxon

pl10 <- ggplot(df_PL_I, aes(x=pmaxs_value, y=max_F, fill=pmaxs_value)) + 
  geom_boxplot() + 
  #stat_compare_means(comparisons = comparisons_PL_I_signif, method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif" ) + #, label.y = c(rep(0.7, times = 18), rep(0.85, times = 6)
  labs(fill = "") + xlab("Maximum number of non-zero coefficients \n max-coef = (1st block, 2nd block)") + ylab(expression(F[max])) +
  ggtitle("Pre-Operative Priority-Lasso-2-blocks models' performance \n with priority sequence 'Clinical > Radiomics'") +
  theme(plot.title = element_text(size = 8), legend.title = element_text(size = 8)) +
  #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() + #geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_brewer(palette = 'Greens') +
  theme(legend.position="none", plot.title = element_text(size = 9, hjust = 0.5), axis.title = element_text(size = 9)) +
  theme(axis.ticks = element_blank()) +
  scale_x_discrete(labels=c("1,5" = "(1,5)", "1,7" = "(1,7)", "1,8" = "(1,8)", "2,7" = "(2,7)", "3,6" = "(3,6)", "5,4" = "(5,4)"))

pl10$layers[[1]]$aes_params$textsize <- 3
pl10

# Conclusion

PT_PL_2blocks_best_model <- PL_I_pmax18

# Priority-LASSO_II -------------------------------------------------------

# Model selection ---------------------------------------------------------
# Choosing between configurations for one algorithm
data_I <- PLII_NO_Norm_I_Pretreat
data_I <- data_I[!( data_I$pmaxs_value==4 | data_I$pmaxs_value==5 |
                      data_I$pmaxs_value==6 |  data_I$pmaxs_value== 7 | data_I$pmaxs_value==8 | 
                      data_I$pmaxs_value== 9 | data_I$pmaxs_value== "1,4,4,4" | data_I$pmaxs_value==3 | data_I$pmaxs_value==2 | data_I$pmaxs_value==1 |
                      data_I$pmaxs_value== "3,3,3,3" |
                      data_I$pmaxs_value== "3,3,3,2" | data_I$pmaxs_value== "3,3,2,2" | data_I$pmaxs_value== "2,4,4,0" | 
                      data_I$pmaxs_value== "4,4,2,0" | data_I$pmaxs_value== "6,4,0,0" | data_I$pmaxs_value== "1,2,3,4" | data_I$pmaxs_value== "1,1,4,4" |
                      data_I$pmaxs_value== "1,3,3,3"), ]  # 

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
data_I[is.na(data_I)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
data_I <- data_I[data_I$total_coeffs < 10,] # COUNTING WITH INTERCEPT

#Encode and label
data_I <- data_I[,c(1,4,13,15:22)]
cols_to_factor <- colnames(data_I[,c(1,2,3)])
cols_to_numeric <- colnames(data_I[,-c(1,2,3)])
data_I[,cols_to_factor] <- lapply(data_I[,cols_to_factor],as.factor)
data_I[,cols_to_numeric] <- lapply(data_I[,cols_to_numeric],as.numeric)


summary_stats_data_I <- tableby(varNames ~ ., data = data_I, control = mycontrols)

summary_stats_data_I
summary(summary_stats_data_I, text = TRUE, title = 'PLII w/o Oversampling, Norm I, dataI', pfootnote = TRUE)

# Group by Label
data_pl1_I <- data_I[data_I$varNames == "pl1",]
data_pl2_I <- data_I[data_I$varNames == "pl2",]
data_pl3_I <- data_I[data_I$varNames == "pl3",]
data_pl4_I <- data_I[data_I$varNames == "pl4",]
data_pl5_I <- data_I[data_I$varNames == "pl5",]
data_pl6_I <- data_I[data_I$varNames == "pl6",]
data_pl7_I <- data_I[data_I$varNames == "pl7",]
data_pl8_I <- data_I[data_I$varNames == "pl8",]
data_pl9_I <- data_I[data_I$varNames == "pl9",]
data_pl10_I <- data_I[data_I$varNames == "pl10",]
data_pl11_I <- data_I[data_I$varNames == "pl11",]
data_pl12_I <- data_I[data_I$varNames == "pl12",]
data_pl13_I <- data_I[data_I$varNames == "pl13",]
data_pl14_I <- data_I[data_I$varNames == "pl14",]
data_pl15_I <- data_I[data_I$varNames == "pl15",]
data_pl16_I <- data_I[data_I$varNames == "pl16",]
data_pl17_I <- data_I[data_I$varNames == "pl17",]
data_pl18_I <- data_I[data_I$varNames == "pl18",]
data_pl19_I <- data_I[data_I$varNames == "pl19",]
data_pl20_I <- data_I[data_I$varNames == "pl20",]
data_pl21_I <- data_I[data_I$varNames == "pl21",]
data_pl22_I <- data_I[data_I$varNames == "pl22",]
data_pl23_I <- data_I[data_I$varNames == "pl23",]
data_pl24_I <- data_I[data_I$varNames == "pl24",]




df3 <- data.frame(Fmax = c(data_I[,"max_F"]))
bpI <- rep(1:24, times = dim(data_I)[1]/24)
df3$Block_Priority <- factor(bpI)

highlighted <- as.logical(df3$Block_Priority == "19" | df3$Block_Priority == "20")
df3$Highlighted = rep(highlighted)


compare_PL2 <- list(c("19", "20"))
pl11 <- ggplot(df3, aes(x=Block_Priority, y=Fmax, fill = Highlighted)) + 
          geom_boxplot() +
          scale_fill_manual(values=c("grey", "green4")) +
          scale_alpha_manual(values=c(1,0.1)) +
          stat_compare_means(comparisons= compare_PL2, method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif" ) +
          xlab("Block Priority Sequence") + ylab(expression(F[max])) +
          ggtitle("Pre-Operative Priority-Lasso-2 (4-blocks) performances by priority sequence") +
          theme(plot.title = element_text(size = 10), legend.title = element_text(size = 9)) +
          #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
          theme_classic() +
          theme(legend.position="none", plot.title = element_text(size = 9, hjust = 0.5), axis.title = element_text(size = 9)) +
          theme(axis.ticks = element_blank()) 

pl11$layers[[2]]$aes_params$textsize <- 3
pl11




Fmax_pl1_I <- data_I[data_I$varNames == "pl1",]$max_F
Fmax_pl2_I <- data_I[data_I$varNames == "pl2",]$max_F
Fmax_pl3_I <- data_I[data_I$varNames == "pl3",]$max_F
Fmax_pl4_I <- data_I[data_I$varNames == "pl4",]$max_F
Fmax_pl5_I <- data_I[data_I$varNames == "pl5",]$max_F
Fmax_pl6_I <- data_I[data_I$varNames == "pl6",]$max_F
Fmax_pl7_I <- data_I[data_I$varNames == "pl7",]$max_F
Fmax_pl8_I <- data_I[data_I$varNames == "pl8",]$max_F
Fmax_pl9_I <- data_I[data_I$varNames == "pl9",]$max_F
Fmax_pl10_I <- data_I[data_I$varNames == "pl10",]$max_F
Fmax_pl11_I <- data_I[data_I$varNames == "pl11",]$max_F
Fmax_pl12_I <- data_I[data_I$varNames == "pl12",]$max_F
Fmax_pl13_I <- data_I[data_I$varNames == "pl13",]$max_F
Fmax_pl14_I <- data_I[data_I$varNames == "pl14",]$max_F
Fmax_pl15_I <- data_I[data_I$varNames == "pl15",]$max_F
Fmax_pl16_I <- data_I[data_I$varNames == "pl16",]$max_F
Fmax_pl17_I <- data_I[data_I$varNames == "pl17",]$max_F
Fmax_pl18_I <- data_I[data_I$varNames == "pl18",]$max_F
Fmax_pl19_I <- data_I[data_I$varNames == "pl19",]$max_F
Fmax_pl20_I <- data_I[data_I$varNames == "pl20",]$max_F
Fmax_pl21_I <- data_I[data_I$varNames == "pl21",]$max_F
Fmax_pl22_I <- data_I[data_I$varNames == "pl22",]$max_F
Fmax_pl23_I <- data_I[data_I$varNames == "pl23",]$max_F
Fmax_pl24_I <- data_I[data_I$varNames == "pl24",]$max_F


df <- cbind(Fmax_pl1_I, Fmax_pl2_I, Fmax_pl3_I, Fmax_pl4_I, Fmax_pl5_I, Fmax_pl6_I, Fmax_pl7_I, Fmax_pl8_I,Fmax_pl9_I,
            Fmax_pl10_I, Fmax_pl11_I, Fmax_pl12_I, Fmax_pl13_I, Fmax_pl14_I, Fmax_pl15_I, Fmax_pl16_I,
            Fmax_pl17_I, Fmax_pl18_I, Fmax_pl19_I, Fmax_pl20_I, Fmax_pl21_I, Fmax_pl22_I,
            Fmax_pl23_I, Fmax_pl24_I)

friedman.test(df)

#Post-Hoc Test, posthoc.friedman.nemenyi.test
#posthoc.friedman.nemenyi.test(df)

wilcox.test(Fmax_pl19_I, Fmax_pl20_I, paired = TRUE)

mycontrols <- tableby.control(test =TRUE, total = FALSE)

summary_stats_pl19_I <- tableby(pmaxs_value ~ ., data = data_pl19_I, control = mycontrols)
summary_stats_pl19_I
summary(summary_stats_pl19_I, text = TRUE, title = 'PLII w/o Oversampling, Norm I, pl19', pfootnote = TRUE)

summary_stats_pl20_I <- tableby(pmaxs_value ~ ., data = data_pl20_I, control = mycontrols)
summary_stats_pl20_I
summary(summary_stats_pl20_I, text = TRUE, title = 'PLI w/o Oversampling, Norm I, pl20', pfootnote = TRUE)
# 
# summary_stats_pl21_II <- tableby(pmaxs_value ~ ., data = data_pl21_I, control = mycontrols)
# summary_stats_pl21_II
# summary(summary_stats_pl21_I, text = TRUE, title = 'PLI w/o Oversampling, Norm II, pl21', pfootnote = TRUE)


# Model selection ---------------------------------------------------------
# Choosing between configurations for one algorithm
PL_II <- PLII_NO_Norm_I_Pretreat
PL_II <- PL_II[!( PL_II$pmaxs_value==4 | PL_II$pmaxs_value==5 |
                    PL_II$pmaxs_value==6 |  PL_II$pmaxs_value== 7 | PL_II$pmaxs_value==8 | 
                    PL_II$pmaxs_value== 9 | PL_II$pmaxs_value== "1,4,4,4" | PL_II$pmaxs_value==3 | PL_II$pmaxs_value==2 | PL_II$pmaxs_value==1 | 
                    PL_II$pmaxs_value== "3,3,3,3" |
                    PL_II$pmaxs_value== "3,3,3,2" | PL_II$pmaxs_value== "3,3,2,2" | PL_II$pmaxs_value== "2,4,4,0" | 
                    PL_II$pmaxs_value== "4,4,2,0" | PL_II$pmaxs_value== "6,4,0,0" | PL_II$pmaxs_value== "1,2,3,4" | PL_II$pmaxs_value== "1,1,4,4" |
                    PL_II$pmaxs_value== "1,3,3,3"), ]  # 

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_II[is.na(PL_II)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_II <- PL_II[PL_II$total_coeffs < 10,] # COUNTING WITH INTERCEPT

#Encode and label
PL_II <- PL_II[,c(1,4,13,15:22)]
cols_to_factor <- colnames(PL_II[,c(1,2,3)])
cols_to_numeric <- colnames(PL_II[,-c(1,2,3)])
PL_II[,cols_to_factor] <- lapply(PL_II[,cols_to_factor],as.factor)
PL_II[,cols_to_numeric] <- lapply(PL_II[,cols_to_numeric],as.numeric)


summary_stats_data_I <- tableby(varNames ~ ., data = PL_II, control = mycontrols)
summary_stats_data_I
summary(summary_stats_data_I, text = TRUE, title = 'PLII w/o Oversampling, Norm I, dataI', pfootnote = TRUE)

# Group by Label

data_pl19_I <- PL_II[PL_II$varNames == "pl19",]
data_pl20_I <- PL_II[PL_II$varNames == "pl20",]

summary_stats_data_pl19_I <- tableby(pmaxs_value ~ ., data = data_pl19_I, control = mycontrols)
summary(summary_stats_data_pl19_I, text = TRUE, title = 'PLII_pl19 w/o Oversampling, Norm I, data_I', pfootnote = TRUE)

write2word(summary_stats_data_pl19_I, "/Users/monicasilva/Documents/Thesis_Writing/RESULTS_raw/PL-II_pl19_I_Kruskal", pfootnote = TRUE)

summary_stats_data_pl20_I <- tableby(pmaxs_value ~ ., data = data_pl20_I, control = mycontrols)
summary(summary_stats_data_pl20_I, text = TRUE, title = 'PLII_pl20 w/o Oversampling, Norm I, data_I', pfootnote = TRUE)


#

# Friedman test
# Performing statistical tests
#PL2_pmax1 <- data_pl19_I[data_pl19_I$pmaxs_value == "1",]
PL2_pmax1601 <- data_pl19_I[data_pl19_I$pmaxs_value == "1,6,0,1",]
PL2_pmax1701 <- data_pl19_I[data_pl19_I$pmaxs_value == "1,7,0,1",]
#PL_pmax2 <- data_pl19_I[data_pl19_I$pmaxs_value == "2",]
PL2_pmax2121 <- data_pl19_I[data_pl19_I$pmaxs_value == "2,1,2,1",]
PL2_pmax2222 <- data_pl19_I[data_pl19_I$pmaxs_value == "2,2,2,2",]
PL2_pmax3210 <- data_pl19_I[data_pl19_I$pmaxs_value == "3,2,1,0",]
PL2_pmax3220 <- data_pl19_I[data_pl19_I$pmaxs_value == "3,2,2,0",]
PL2_pmax3300 <- data_pl19_I[data_pl19_I$pmaxs_value == "3,3,0,0",]
PL2_pmax3320 <- data_pl19_I[data_pl19_I$pmaxs_value == "3,3,2,0",]
PL2_pmax3330 <- data_pl19_I[data_pl19_I$pmaxs_value == "3,3,3,0",]


#PL2_pmax1_Fmax <- PL2_pmax1$max_F
PL2_pmax1601_Fmax <- PL2_pmax1601$max_F
PL2_pmax1701_Fmax <- PL2_pmax1701$max_F
#PL_pmax2_Fmax <- PL_pmax2$max_F
PL2_pmax2121_Fmax <- PL2_pmax2121$max_F
PL2_pmax2222_Fmax <- PL2_pmax2222$max_F
PL2_pmax3210_Fmax <- PL2_pmax3210$max_F
PL2_pmax3220_Fmax <- PL2_pmax3220$max_F
PL2_pmax3300_Fmax <- PL2_pmax3300$max_F
PL2_pmax3320_Fmax <- PL2_pmax3320$max_F
PL2_pmax3330_Fmax <- PL2_pmax3330$max_F



df <- cbind(PL2_pmax1601_Fmax, PL2_pmax1701_Fmax, PL2_pmax2121_Fmax, PL2_pmax2222_Fmax,
            PL2_pmax3210_Fmax, PL2_pmax3220_Fmax, PL2_pmax3300_Fmax, PL2_pmax3320_Fmax, PL2_pmax3330_Fmax)

library(PMCMR)
friedman.test(df)
posthoc.friedman.nemenyi.test(y=df)

# library("ggpubr")
# ggboxplot(data_pl19_I, x = "pmaxs_value", y = "max_F", 
#           color = "pmaxs_value", #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#           #order = c("ctrl", "trt1", "trt2"),
#           ylab = "Fmax", xlab = "max.coefficients")

comparisons_PL_II_nemenyi_sign = list( 
  c("1,7,0,1", "2,1,2,1"), # ****
  c("1,7,0,1", "2,2,2,2"), #*
  c("1,7,0,1", "3,2,1,0"), #**
  c("1,7,0,1", "3,2,2,0"), #**
  c("1,7,0,1", "3,3,0,0"), #* 
  c("1,6,0,1", "2,1,2,1"), # ****
  c("1,6,0,1", "3,2,1,0"), #*
  c("1,6,0,1", "3,2,2,0")) #*




pl12 <- ggplot(data_pl19_I, aes(x=pmaxs_value, y=max_F, fill=pmaxs_value)) + 
          geom_boxplot() + 
          stat_compare_means(comparisons = comparisons_PL_II_nemenyi_sign, method = "wilcox.test", paired = TRUE, size = 4, label = "p.signif" ) + #, label.y = c(rep(0.7, times = 18), rep(0.85, times = 6)
          labs(fill = "") + xlab("Maximum number of non-zero coefficients for each block \n max-coef = (1st, 2nd, 3rd, 4th block)") + ylab(expression(F[max])) +
          ggtitle("Comparison of performance of Pre-Operative Priority-Lasso-4-blocks models") +
          #scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
          theme_classic() + #geom_jitter(color="black", size=0.4, alpha=0.9) +
          scale_fill_brewer(palette = 'Greens') +
          theme(legend.position="none", plot.title = element_text(size = 9, hjust = 0.5), axis.title = element_text(size = 9)) +
          theme(axis.ticks = element_blank()) +
          scale_x_discrete(labels=c("1,6,0,1" = "(1,6,0,1)", 
                                    "1,7,0,1" = "(1,7,0,1)",
                                    "2,1,2,1" = "(2,1,2,1)", 
                                    "2,2,2,2" = "(2,2,2,2)",
                                    "3,2,1,0" = "(3,2,1,0)", 
                                    "3,2,2,0" = "(3,2,2,0)",
                                    "3,3,0,0" = "(3,3,0,0)",
                                    "3,3,2,0" = "(3,3,2,0)",
                                    "3,3,3,0" = "(3,3,3,0)"))

pl12$layers[[2]]$aes_params$textsize <- 3
pl12


#Conclusion

PT_PL_4blocks_best_model <- PL2_pmax1701


# Plot
PT_Lasso_best_model
PT_PL_2blocks_best_model
PT_PL_4blocks_best_model


# PART IV ---------------------------------------------------------------

# Choosing between algorithms

best_LASSO_I <- PT_Lasso_best_model #pmax9
best_PL_I <-  PT_PL_2blocks_best_model #(1,8)
best_PL_II <- PT_PL_4blocks_best_model #(1,7,0,1)

df_best_Fmax <- cbind(best_LASSO_I$F_max, best_PL_I$max_F, best_PL_II$max_F)
colnames(df_best_Fmax) <- c("Lasso pmax=9", "Priority-Lasso-2-blocks max.coef=(1,8)", "Priority-Lasso-4-blocks max.coef=(1,7,0,1)")
df_best_Fmax <- as.data.frame(df_best_Fmax)


# create a data frame
reps = dim(best_LASSO_I)[1]
#normalization = rep(c("Prostate-Norm", "Whole-image Norm"), each = reps)
model = rep(c("Lasso pmax = 9", "Priority-Lasso-2-blocks max-coef = (1,8)", "Priority-Lasso-4-blocks max-coef = (1,7,0,1)"), each = reps)


df_best_Fmax <- data.frame(Fmax = c(best_LASSO_I[,"F_max"], best_PL_I[,"max_F"], best_PL_II[,"max_F"]))
#df$Normalization = normalization
df_best_Fmax$Model = model
#df$Algorithm = rep("Lasso", each = reps*2)


# Statistical comparison
best_LASSO_I_Fmax <- best_LASSO_I$F_max
best_PL_I_Fmax <- best_PL_I$max_F
best_PL_II_Fmax <- best_PL_II$max_F

df <- cbind(best_LASSO_I_Fmax, best_PL_I_Fmax, best_PL_II_Fmax)

friedman.test(df)
posthoc.friedman.nemenyi.test(df)



my_comparisons = list( c("Lasso pmax = 9", "Priority-Lasso-4-blocks max-coef = (1,7,0,1)"), c("Lasso pmax = 9", "Priority-Lasso-2-blocks max-coef = (1,8)"), c("Priority-Lasso-2-blocks max-coef = (1,8)", "Priority-Lasso-4-blocks max-coef = (1,7,0,1)"))

pl13 <- ggplot(df_best_Fmax, aes(x=Model, y=Fmax, fill = Model)) + 
  geom_boxplot() +
  theme(legend.position="none") +
  labs(fill = "") +
  xlab("Best performing Pre-Operative models") + 
  ylab(expression(F[max])) +
  ggtitle("Pre-Operative best performing models comparison") +
  scale_fill_brewer(palette="Greens") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired = TRUE, label = "p.signif")  +
  theme(plot.title = element_text(size = 8), legend.title = element_text(size = 8)) +
  #geom_jitter(color="black", size=0.7, alpha=0.9) +
  theme_classic() + #geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme(legend.position="none")

pl13$layers[[2]]$aes_params$textsize <- 3
pl13


mycontrols <- tableby.control(test =TRUE, total = FALSE, numeric.test = "kwt" )

summary_stats_BEST <- tableby(Model ~ ., data = df_best_Fmax, control = mycontrols)
summary_stats_BEST
summary(summary_stats_BEST, text = TRUE, title = 'I', pfootnote = TRUE)

write2word(summary_stats_BEST, "/Users/monicasilva/Documents/Thesis_Writing/RESULTS_raw/best_3_stats_Pre-biopsy" )









# SMOTE -------------------------------------------------------------------
#Selecting the best candidates, SMOTE version
# LASSO SMOTE:  -----------------------------------------------------------


LASSO_I_SMOTE <- Lasso_SMOTE_Norm_I_Prebio

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
LASSO_I_SMOTE[is.na(LASSO_I_SMOTE)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
LASSO_I_SMOTE <- LASSO_I_SMOTE[LASSO_I_SMOTE$nzcoef <= 10,] # COUNTING WITH INTERCEPT
LASSO_I_SMOTE <- LASSO_I_SMOTE[(LASSO_I_SMOTE$pmaxs_value != "NL"),]

#Encode and label
#(dim(data_LASSO)[2])
cols_to_factor <- colnames(LASSO_I_SMOTE[,c(1,2,3,15)])
cols_to_numeric <- colnames(LASSO_I_SMOTE[,-c(1,2,3,15)])
LASSO_I_SMOTE[,cols_to_factor] <- lapply(LASSO_I_SMOTE[,cols_to_factor],as.factor)
LASSO_I_SMOTE[,cols_to_numeric] <- lapply(LASSO_I_SMOTE[,cols_to_numeric],as.numeric)

# Group by Label
LASSO_I_SMOTE <- LASSO_I_SMOTE[-c(2,3,12:14,16)]
mycontrols <- tableby.control(test =TRUE, total = FALSE)

summary_stats <- tableby(pmaxs_value ~ ., data = LASSO_I_SMOTE, control = mycontrols)
summary_stats
summary(summary_stats, text = TRUE, title = 'LASSO w/ SMOTE, Norm I', pfootnote = TRUE)


# Performing statistical tests
LASSO_SMOTE_pmax5 <- LASSO_I_SMOTE[LASSO_I_SMOTE$pmaxs_value == 5,]
LASSO_SMOTE_pmax6 <- LASSO_I_SMOTE[LASSO_I_SMOTE$pmaxs_value == 6,]
LASSO_SMOTE_pmax7 <- LASSO_I_SMOTE[LASSO_I_SMOTE$pmaxs_value == 7,]
LASSO_SMOTE_pmax8 <- LASSO_I_SMOTE[LASSO_I_SMOTE$pmaxs_value == 8,]
LASSO_SMOTE_pmax9 <- LASSO_I_SMOTE[LASSO_I_SMOTE$pmaxs_value == 9,]

t.test(LASSO_pmax9$F_max, LASSO_SMOTE_pmax9$F_max, paired = FALSE, alternative = "two.sided")
t.test(LASSO_pmax9$F_max, LASSO_SMOTE_pmax9$F_max, paired = FALSE, alternative = "less")

# CONCLUSION: performance w/ SMOTE is sign. better than no SMOTE.

t.test(LASSO_I$F_max, LASSO_I_SMOTE$F_max, paired = FALSE, alternative = "two.sided")
t.test(LASSO_I$F_max, LASSO_I_SMOTE$F_max, paired = FALSE, alternative = "less")

# CONCLUSION: actually, performance w/ SMOTE is sign. better than no SMOTE; regardless of algorithm configuration.




# PL_I SMOTE --------------------------------------------------------------


PL_I_SMOTE <- PLI_SMOTE_Norm_I_Prebio

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_I_SMOTE[is.na(PL_I_SMOTE)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_I_SMOTE <- PL_I_SMOTE[PL_I_SMOTE$total_coeffs <= 10,] # COUNTING WITH INTERCEPT

#Encode and label
PL_I_SMOTE <- PL_I_SMOTE[,c(1,4,9,11:18)]
cols_to_factor <- colnames(PL_I_SMOTE[,c(1,2,3)])
cols_to_numeric <- colnames(PL_I_SMOTE[,-c(1,2,3)])
PL_I_SMOTE[,cols_to_factor] <- lapply(PL_I_SMOTE[,cols_to_factor],as.factor)
PL_I_SMOTE[,cols_to_numeric] <- lapply(PL_I_SMOTE[,cols_to_numeric],as.numeric)

# Group by Label

data_SMOTE_pl1_I <- PL_I_SMOTE[PL_I_SMOTE$varNames == "pl1",]
data_SMOTE_pl2_I <- PL_I_SMOTE[PL_I_SMOTE$varNames == "pl2",]

mycontrols <- tableby.control(test =TRUE, total = FALSE)

summary_stats_pl1_I <- tableby(pmaxs_value ~ ., data = data_SMOTE_pl1_I, control = mycontrols)
summary_stats_pl1_I
summary(summary_stats_pl1_I, text = TRUE, title = 'PLI w/o Oversampling, Norm I, pl1', pfootnote = TRUE)

summary_stats_pl2_I <- tableby(pmaxs_value ~ ., data = data_SMOTE_pl2_I, control = mycontrols)
summary_stats_pl2_I
summary(summary_stats_pl2_I, text = TRUE, title = 'PLI w/o Oversampling, Norm I, pl2', pfootnote = TRUE)


#Group by pmaxs

pl2_I_SMOTE_pmax5 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == 5,]
pl2_I_SMOTE_pmax6 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == 6,]
pl2_I_SMOTE_pmax7 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == 7,]
pl2_I_SMOTE_pmax8 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == 8,]
pl2_I_SMOTE_pmax9 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == 9,]
pl2_I_SMOTE_pmax9 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == 9,]
pl2_I_SMOTE_maxcoef15 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == "1,5",]
pl2_I_SMOTE_maxcoef17 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == "1,7",]
pl2_I_SMOTE_maxcoef18 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == "1,8",]
pl2_I_SMOTE_maxcoef27 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == "2,7",]
pl2_I_SMOTE_maxcoef36 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == "3,6",]
pl2_I_SMOTE_maxcoef54 <- data_SMOTE_pl2_I[data_SMOTE_pl2_I$pmaxs_value == "5,4",]


# Statistical tests -------------------------------------------------------

t.test(pl2_I_maxcoef18$max_F, pl2_I_SMOTE_maxcoef18$max_F, paired = TRUE, alternative = "two.sided")
t.test(pl2_I_maxcoef18$max_F, pl2_I_maxcoef27$max_F, paired = TRUE, alternative = "two.sided")

res.aov <- aov(max_F ~ pmaxs_value, data = data_SMOTE_pl2_I)
summary(res.aov) 

# CONCLUSION: results are sign. different accordingly to their configuration.


library("ggpubr")
ggboxplot(data_SMOTE_pl2_I, x = "pmaxs_value", y = "max_F", 
          color = "pmaxs_value", #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          #order = c("ctrl", "trt1", "trt2"),
          ylab = "Fmax", xlab = "max.coefficients")






# PL_II_SMOTE -------------------------------------------------------------

PL_II_SMOTE <- PLII_SMOTE_Norm_I_Prebio
PL_II_SMOTE <- PL_II_SMOTE[!( PL_II_SMOTE$pmaxs_value==4 | PL_II_SMOTE$pmaxs_value==5 |
                                PL_II_SMOTE$pmaxs_value==6 |  PL_II_SMOTE$pmaxs_value== 7 | PL_II_SMOTE$pmaxs_value==8 | 
                                PL_II_SMOTE$pmaxs_value== 9 | PL_II_SMOTE$pmaxs_value== "1,4,4,4" | PL_II_SMOTE$pmaxs_value==3 |
                                PL_II_SMOTE$pmaxs_value== "3,3,3,3" | PL_II_SMOTE$pmaxs_value== "3,3,3,2") ,] 

# CLEANING
#data_LASSO <- data_LASSO[rowSums(is.na(data_LASSO)) == 0,] #keep only patients without NA values
PL_II_SMOTE[is.na(PL_II_SMOTE)] <- 0
# APPLYING EXCLUSION CRITERIA, keeping only models with nzcoef < 9
PL_II_SMOTE <- PL_II_SMOTE[PL_II_SMOTE$total_coeffs < 10,] # COUNTING WITH INTERCEPT

#Encode and label
PL_II_SMOTE <- PL_II_SMOTE[,c(1,4,13,15:22)]
cols_to_factor <- colnames(PL_II_SMOTE[,c(1,2,3)])
cols_to_numeric <- colnames(PL_II_SMOTE[,-c(1,2,3)])
PL_II_SMOTE[,cols_to_factor] <- lapply(PL_II_SMOTE[,cols_to_factor],as.factor)
PL_II_SMOTE[,cols_to_numeric] <- lapply(PL_II_SMOTE[,cols_to_numeric],as.numeric)


summary_stats_data_I <- tableby(varNames ~ ., data = PL_II_SMOTE, control = mycontrols)
summary_stats_data_I
summary(summary_stats_data_I, text = TRUE, title = 'PLII w/ SMOTE, Norm I, dataI', pfootnote = TRUE) 

# CONCLUSION #pll9 and pl20 are the best

# Group by Label

data_pl19_I_SMOTE <- PL_II_SMOTE[PL_II_SMOTE$varNames == "pl19",]
data_pl20_I_SMOTE <- PL_II_SMOTE[PL_II_SMOTE$varNames == "pl20",]

summary_stats_data_pl19_I_SMOTE <- tableby(pmaxs_value ~ ., data = data_pl19_I_SMOTE, control = mycontrols)
summary(summary_stats_data_pl19_I_SMOTE, text = TRUE, title = 'PLII_pl19 w/o Oversampling, Norm I, data_I', pfootnote = TRUE)

summary_stats_data_pl20_I_SMOTE <- tableby(pmaxs_value ~ ., data = data_pl20_I_SMOTE, control = mycontrols)
summary(summary_stats_data_pl20_I_SMOTE, text = TRUE, title = 'PLII_pl20 w/o Oversampling, Norm I, data_I', pfootnote = TRUE)


# CONCLUSION: PL19 1,7,0,1 performs better than PL20 1,7,0,1

library("ggpubr")
ggboxplot(data_pl19_I_SMOTE, x = "pmaxs_value", y = "max_F", 
          color = "pmaxs_value", #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          #order = c("ctrl", "trt1", "trt2"),
          ylab = "Fmax", xlab = "max.coefficients")


res.aov <- aov(max_F ~ pmaxs_value, data = data_pl19_I_SMOTE)
summary(res.aov)

# CONCLUSION: performance results are sig. different for different algorithm configurations

pairwise.t.test(data_pl19_I$max_F, data_pl19_I$pmaxs_value,
                p.adjust.method = "BH", pool.sd = FALSE)

data_pl19_I_SMOTE_best <- data_pl19_I_SMOTE[data_pl19_I_SMOTE$pmaxs_value == "1,7,0,1",]
data_pl20_I_SMOTE_best <- data_pl20_I_SMOTE[data_pl20_I_SMOTE$pmaxs_value == "1,7,0,1",]

t.test(data_pl19_I_best$max_F, data_pl19_I_SMOTE_best$max_F, paired = FALSE, alternative = "two.sided") 
t.test(data_pl19_I_best$max_F, data_pl19_I_SMOTE_best$max_F, paired = FALSE, alternative = "less") 
# CONCLUSION: SMOTE results are sig. different from NO SMOTE results, for the best configuration.

t.test(data_pl19_I_best$max_F, data_pl20_I_best$max_F, paired = FALSE, alternative = "two.sided") # CONCLUSION: PL19 and PL20, SMOTE, performance, is not sig. different (Fmax 0.882 vs 0.863)



t.test(PL_II$max_F, PL_II_SMOTE$max_F, paired = FALSE, alternative = "two.sided") 
t.test(PL_II$max_F, PL_II_SMOTE$max_F, paired = FALSE, alternative = "less") 
# CONCLUSION: Fmax performance is sig. different for NON SMOTE and SMOTE conditions, where SMOTE has higher mean values for Fmax metric.





# OVERSAMPLING BEST CLASSIFIER

data_pl19_I_SMOTE_best <- data_pl19_I_SMOTE[data_pl19_I_SMOTE$pmaxs_value == "1,7,0,1",]

wilcox.test(data_pl19_I_best$max_F, data_pl19_I_SMOTE_best$max_F, paired = FALSE)
wilcox.test(data_pl19_I_best$max_F, data_pl19_I_SMOTE_best$max_F, paired = FALSE, alternative = "less")



# best_LASSO_I_SMOTE <- LASSO_SMOTE_pmax9
# best_PL_I_SMOTE <- pl2_I_SMOTE_maxcoef18
# best_PL_II_SMOTE <- data_pl19_I_SMOTE_best


summary(best_LASSO_I)
summary(best_LASSO_I_SMOTE)
summary(best_PL_I)
summary(best_PL_I_SMOTE)
summary(best_PL_II)
summary(best_PL_II_SMOTE)

