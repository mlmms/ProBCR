setwd("/Users/monicasilva/Documents/2020_RStudio_directory")
SEED <- c(42,67,88,709,92,796,516,281,883,942,730,685,465,242)
BASE <- "/Users/monicasilva/Documents/Thesis_Writing/4_Results/4_Results_Normalization_II/1_PREBIOPSY/PLII_seedingPAIRS_SMOTE_"

lasso1 <- "pmax6"
lasso2 <- "pmax7"
lasso1_smote <- "pmax9"
lasso2_smote <- "pmax10"

pl1_1 <- "maxcoef19"
pl1_2 <- "maxcoef46"
pl1_1_smote <- "pmax7"
pl1_2_smote <- "maxcoef28"

pl2_1 <- "pmax9"
pl2_2 <- "maxcoef1601"
pl2_1_smote <- "maxcoef1701"
pl2_2_smote <- "maxcoef17013330"

#pmaxs_value <- c(5,6,7,8,9) #
#pmaxs_value <- rbind(c(1,8), c(1,5), c(1,7), c(5,4), c(2,7), c(3,6))
#pmaxs_value <- c(1,2,3,4,5,6,7,8,9) #
pmaxs_value <- rbind(c(3,3,3,3), c(3,3,3,2), c(3,3,3,0), c(3,3,2,2), c(3,3,2,0), c(3,3,0,0),
                     c(3,2,2,0), c(3,2,1,0), c(2,4,4,0), c(4,4,2,0), c(6,4,0,0), c(1,6,0,1), 
                     c(1,4,4,4), c(1,2,3,4), c(1,3,3,3), c(1,1,4,4), c(2,2,2,2), c(2,1,2,1), c(1,7,0,1))

# F_max <- vector()
# P_max <- vector()
# R_max <- vector()
# F0_max <- vector()
# P0_max <- vector()
# R0_max <- vector()
# optimal_cutpoint <- vector()
# cvlo <- vector()
# cvup <- vector()
# measure_lmin <- vector()
# std_measure_lmin <- vector()
# nzcoef <- vector()
# roc_train_auc <- vector()
# results_LASSO <- 
#TO WRITE CSV
# varNames <- vector() #will store the objects names of the priority lasso's
# aucsTrain <- vector()
# n_coeffs1 <- vector()
# n_coeffs2 <- vector()
# total_coeffs <- vector()
# str <- vector()
# max_cvm <- vector() #
# max_F <- vector()
# max_P <- vector()
# max_R <- vector()
# max_F0 <- vector()
# max_P0 <- vector()
# max_R0 <- vector()
# auc_pl_roc_train <- vector()
# results_PLI <- NULL
# optimal_cutpoint <- vector()

#TO WRITE CSV
varNames <- vector() #will store the objects names of the priority lasso's
aucsTrain <- vector()
stdTrain <- vector()
aucsTest <- vector()
stdTest <- vector()
cisTrain_l <- vector()
cisTrain_u <- vector()
cisTest_l <- vector()
cisTest_u <- vector()
n_coeffs1 <- vector()
n_coeffs2 <- vector()
n_coeffs3 <- vector()
n_coeffs4 <- vector()
total_coeffs <- vector()
str <- vector()
max_cvm <- vector()
max_F <- vector()
max_P <- vector()
max_R <- vector()
max_F0 <- vector()
max_P0 <- vector()
max_R0 <- vector()
pl_roc_train_auc <- vector()
results_PLII <- NULL
optimal_cutpoint <- vector()


perm <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

for (i in 1:length(pmaxs_value)) {
  for (j in 1:length(SEED)) { #length(SEED)
    set.seed(SEED[j])
    
    library(caret)
    library(lattice)
    library(ggplot2)
    library(pROC)
    library(writexl)
    library(precrec)
    library(dplyr)
    library(plyr)
    library(ROCR)
    library(onehot)
    library(binda)
    library(unbalanced) #masks objects ROCR::performance and caret::train
    library(lasso2)
    library(cutpointr) #auc objects is masked from glmnet, precrec, pROC (auc and roc), precision, recall, sensitivity, specificity from caret.
    library(glmnet) # library including functions for finding the lasso
    # glmnet requires design matrix and response variable\
    library(DMwR) #masks object plyr::join
    library(prioritylasso)
    
    
    # Data Loading -------------------------------------------------------------
    clin_path <- "/Users/monicasilva/Documents/1_Original_data/Clinical_data/ProBCR_2020_preBIOPSY_clinical_variables_93p_simple.csv"
    t2_ii_path <- "/Users/monicasilva/Documents/Final_Radiomics_Extraction_Force2DTRUE/results_2020_t2_II_bin_method_93p.csv"
    b1000_ii_path <- "/Users/monicasilva/Documents/Final_Radiomics_Extraction_Force2DTRUE/results_2020_b1000_II_bin_method_93.csv"
    adc_ii_path <- "/Users/monicasilva/Documents/Final_Radiomics_Extraction_Force2DTRUE/results_2020_adc_II_bin_method_93p.csv"
    
    clin_data <- read.csv(clin_path, header=TRUE, sep=";")
    
    #Clinical data: encode and label
    cols_to_factor <- colnames(clin_data[,c(1:4)])
    cols_to_numeric <- colnames(clin_data[,5])
    clin_data[,cols_to_factor] <- lapply(clin_data[,cols_to_factor],as.factor)
    clin_data[,cols_to_numeric] <- lapply(clin_data[,cols_to_numeric],as.numeric)
    #sapply(clin_data, class)
    #sapply(clin_data,summary)
    
    # #Defining model's COVARIATES: ----------------------------------
    clin_continuous_variables <- clin_data[,5]
    UID <- clin_data[,1]
    Label <- clin_data[,2]
    UID_Label_df <- clin_data[,c(1,2)]
    
    Clin_preopPSA <- clin_data[,5]
    Clin_preopPSA <- as.data.frame(Clin_preopPSA)
    Clin_preopPSA <- cbind(UID, Clin_preopPSA)
    
    Clin_preopPSA.10 <- clin_data[,3]
    Clin_preopPSA.10 <- as.data.frame(Clin_preopPSA.10)
    
    Clin_preopPSA_grading_0lessthan10_1_between10and20_2_higher20 <- clin_data[,4]
    Clin_preopPSA_grading_0lessthan10_1_between10and20_2_higher20 <- as.data.frame(Clin_preopPSA_grading_0lessthan10_1_between10and20_2_higher20)
    
    t2_ii_df <- read.csv(t2_ii_path,header=TRUE,fill=FALSE)
    t2_ii_df <- t2_ii_df[,-c(1,3:39)]
    
    b1000_ii_df <- read.csv(b1000_ii_path,header=TRUE,fill=FALSE)
    b1000_ii_df <- b1000_ii_df[,-c(1,3:39)]
    
    adc_ii_df <- read.csv(adc_ii_path,header=TRUE,fill=FALSE)
    adc_ii_df <- adc_ii_df[,-c(1,3:39)]
    
    radiomics_df <- merge(UID_Label_df, t2_ii_df, by.x = 'UID', by.y = 'Case')
    radiomics_df <- merge(radiomics_df, b1000_ii_df, by.x = 'UID', by.y = 'Case')
    radiomics_df <- merge(radiomics_df, adc_ii_df, by.x = 'UID', by.y = 'Case')
    #radiomics_df <- radiomics_df[rowSums(is.na(radiomics_df)) == 0,] #keep only patients without NA values
    
    # Data One-Hot Encoding ---------------------------------------------------------------
    to_one_hot_df <- clin_data[,c(3,4)]
    
    encoder <- onehot(to_one_hot_df)
    encoded_matrix <- predict(encoder, to_one_hot_df)
    encoded_df <- as.data.frame(encoded_matrix)
    
    clin_categorical_variables_multiplelevels_onehotencoded <- cbind(UID, encoded_df)
    
    clin_data <- clin_data[,-c(3,4)]
    clin_df_ii <- merge(clin_data, clin_categorical_variables_multiplelevels_onehotencoded, by = 'UID')
    
    #Clinical data: encode
    cols_to_factor <- colnames(clin_df_ii[,c(2,4:8)])
    cols_to_numeric <- colnames(clin_df_ii[,3])
    clin_df_ii[,cols_to_factor] <- lapply(clin_df_ii[,cols_to_factor],as.factor)
    clin_df_ii[,cols_to_numeric] <- lapply(clin_df_ii[,cols_to_numeric],as.numeric)
    
    # DataMerging -------------------------------------------------------------
    ## Final df
    df_ii <- merge(clin_df_ii, radiomics_df, by = 'UID')
    df_ii <- df_ii[rowSums(is.na(df_ii)) == 0,] #keep only patients without NA values
    
    
    # StratifiedTrainTestSplitting ------------------------------------------------------
    train_perc <- 0.70
    test_perc <- 0.30
    
    #Data shuffling (by row)
    set.seed(SEED[j])
    shuffled_df_ii <- df_ii[sample(nrow(df_ii)),]
    
    trainIndex <- createDataPartition(y=df_ii$Label.x, p=train_perc, list=FALSE) # the random sampling is done within the levels of y when y is a factor in an attempt to balance the class distributions within the splits.
    
    #Pre-standardization
    index_Label_x <- which(names(df_ii) == "Label.x")
    index_Label_y <- which(names(df_ii) == "Label.y")
    
    #TRAIN
    #features to standardize (continuous variables)
    x_train_to_std <- shuffled_df_ii[trainIndex, c(3,(index_Label_y+1):length(shuffled_df_ii))]
    #features to not standardize (categorical variables)
    x_train_to_not_std <- shuffled_df_ii[trainIndex, c(1,2,4:(index_Label_y-1))]
    
    #TEST
    #features to standardize (continuous variables)
    x_test_to_std <- shuffled_df_ii[-trainIndex, c(3,(index_Label_y+1):length(shuffled_df_ii))]
    #features to not standardize (categorical variables)
    x_test_to_not_std <- shuffled_df_ii[-trainIndex, c(1,2,4:(index_Label_y-1))]
    
    y_train <- shuffled_df_ii[trainIndex, c(1,index_Label_x)] #UID and Label
    y_test <- shuffled_df_ii[-trainIndex, c(1,index_Label_x)] #UID and Label
    
    
    # DataStandardization -----------------------------------------------------
    x_train_num <- (apply(x_train_to_std,2,as.numeric)) #2= apply_over_columns
    #y_train_num <- shuffled_df_ii[trainIndex,2] # NEW
    y_train_num <- y_train[,2]
    
    standardize <- preProcess(x_train_num, method=c("center", "scale","nzv", "zv"))
    x_train_scaled <- predict(standardize, x_train_num)
    
    x_test_num <- (apply(x_test_to_std,2,as.numeric)) #2= apply_over_columns
    x_test_scaled <- predict(standardize, x_test_num)
    y_test_num <- y_test[,2]
    
    # Merge all variables to train and to test sets (standardized and not) #!
    x_train_f <- cbind(x_train_to_not_std, x_train_scaled)
    x_test_f <- cbind(x_test_to_not_std, x_test_scaled)
    
    index_Label_train_x <- which(names(x_train_f) == "Label.x")
    index_Label_train_y <- which(names(x_train_f) == "Label.y")
    
    index_Label_test_x <- which(names(x_test_f) == "Label.x")
    index_Label_test_y <- which(names(x_test_f) == "Label.y")
    
    
    #remove UID and Label x and y
    x_train_f <- x_train_f[,-c(1,index_Label_train_x, index_Label_train_y)] 
    x_test_f <- x_test_f[,-c(1,index_Label_test_x, index_Label_test_y)]
    
    
    # Cross Validation
    set.seed(SEED[j])
    CVFOLD = 4
    # ADDED: create folds AND 'foldid' cv.glmnet parameter: to use stratified k folds
    flds <- createFolds(y_train_num, k = CVFOLD, list = TRUE, returnTrain = FALSE)
    foldids = rep(1,length(y_train_num))
    foldids[flds$Fold2] = 2
    foldids[flds$Fold3] = 3
    foldids[flds$Fold4] = 4
    
    
    #define outcome and predictors
    df_train_out <- apply(as.matrix(y_train_num), 2, as.numeric)
    df_train_pred <- apply(as.matrix(x_train_f), 2, as.numeric) # matrix with the predictors without label
    
    df_test_out <- apply(as.matrix(y_test_num), 2, as.numeric)
    df_test_pred <- apply(as.matrix(x_test_f), 2, as.numeric) # matrix with the predictors without label
    
    #------------------------------------------------------------------- LASSO 1 ------------------------------------------------------------
    
    
    #     results_path <- 
    #       paste0(BASE,"Lasso_nooversampling_testing.xlsx")
    #     
    #     
    #     
    #     # ChoosingLambda: unweighted data
    #     # A ) unweighted data: Observation weights; defaults to 1 per observation
    #     set.seed(SEED[j])
    #     cvlasso<-cv.glmnet(df_train_pred,as.factor(df_train_out),
    #                        alpha=1, family="binomial",
    #                        type.measure="auc", nlambda = 100,
    #                        standardize=FALSE, nfolds = 4, foldid = foldids, pmax = pmaxs_value[i] ) #weights = 
    #     
    #     #print(cvlasso)
    #     ##plotcvlasso)
    #     cvlasso$lambda.min
    #     
    #     #coefficients
    #     myCoefs <- coef(cvlasso, s="lambda.min");
    #     myCoefs[which(myCoefs != 0 ) ]; #coefficients: intercept included
    #     myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ];  #feature names: intercept included
    #     ## Asseble into a data.frame
    #     myResults <- data.frame(
    #       features = myCoefs@Dimnames[[1]][ which(myCoefs != 0 ) ], #intercept included
    #       coefs    = myCoefs[ which(myCoefs != 0 ) ]  #intercept included
    #     )
    #     print(myResults)
    #     print(paste0("number_coefs = ", nrow(myResults)))
    #     
    #     lmin <- cvlasso[["lambda.min"]]
    #     index_lmin <- which(cvlasso[["lambda"]] == lmin)
    #     measure_lmin[i] <- cvlasso[["cvm"]][index_lmin]
    #     std_measure_lmin[i] <- cvlasso[["cvsd"]][index_lmin]
    #     cvup[i] <- cvlasso[["cvup"]][index_lmin]
    #     cvlo[i] <- cvlasso[["cvlo"]][index_lmin]
    #     print(paste0("auclmin = ", round(measure_lmin,3), " std = ", round(std_measure_lmin,3)))
    #     
    #     
    #     #Predictions after CV: use predict()
    #     #1) lambda.min
    #     #Predictions after CV: use predict() #1) lambda.min
    #     
    #     #------ With OCP -----#
    #     response_predictions_train_lmin <- data.matrix(predict(cvlasso$glmnet.fit, newx = df_train_pred, type = "response", s=cvlasso$lambda.min))
    #     # response_predictions_test_lmin <- data.matrix(predict(cvlasso$glmnet.fit, newx = df_test_pred, type = "response", s=cvlasso$lambda.min))
    #     
    #     roc_train <- pROC::roc(response=y_train_num, predictor=response_predictions_train_lmin[,1])
    #     ##plotroc_train)
    #     roc_train["auc"]
    #     
    #     # Finding OCP:
    #     metric_to_optimize <- "F1_score"  # metric_to_optimize <- "Recall" #metric_to_optimize <- "ppv" #metric_to_optimize <- "roc01" #metric_to_optimize <- "cohens_kappa"
    #     
    #     set.seed(SEED[j])
    #     cutpoint <- cutpointr(x = response_predictions_train_lmin, 
    #                           class = as.factor(df_train_out),
    #                           use_midpoints = TRUE, 
    #                           maximize_boot_metric = metric_to_optimize)
    #     
    #     ##plotcutpoint)
    #     optimal_cutpoint[i] <- cutpoint$optimal_cutpoint
    #     new_predictions_train_lmin <- dichotomize(response_predictions_train_lmin, cutpoint$optimal_cutpoint)
    #     OCP_confusion_matrix_train_lmin <- confusionMatrix(reference =y_train_num, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "1")
    #     OCP_confusion_matrix_train_lmin
    #     
    #     roc_train_auc[i] <- roc_train[["auc"]][1]
    #     F_max[i] <- OCP_confusion_matrix_train_lmin[["byClass"]][["F1"]]
    #     P_max[i] <- OCP_confusion_matrix_train_lmin[["byClass"]][["Precision"]]
    #     R_max[i] <- OCP_confusion_matrix_train_lmin[["byClass"]][["Recall"]]
    #     
    #     signames <- myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ]
    #     nzcoef[i] <- length(myCoefs[which(myCoefs != 0 ) ])
    #     
    #     #new: Negative
    #     OCP_confusion_matrix_train_lmin_negative <- confusionMatrix(reference =y_train_num, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "0")
    #     OCP_confusion_matrix_train_lmin_negative
    #     
    #     F0_max[i] <- OCP_confusion_matrix_train_lmin_negative[["byClass"]][["F1"]]
    #     P0_max[i] <- OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Precision"]]
    #     R0_max[i] <- OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Recall"]]
    #     
    #     # new_predictions_test_lmin <- dichotomize(response_predictions_test_lmin, cutpoint$optimal_cutpoint)
    #     # OCP_confusion_matrix_test_lmin <- confusionMatrix(reference = y_test_num, data = as.factor(new_predictions_test_lmin), mode = "everything", positive = "1")
    #     # OCP_confusion_matrix_test_lmin
    #     
    #     results_LASSO <- rbind(results_LASSO, data.frame(pmaxs_value[i], SEED[j], "LASSO", cvlo[i], cvup[i], measure_lmin[i], std_measure_lmin[i],
    #                                                      + roc_train_auc[i], F_max[i], P_max[i], R_max[i], F0_max[i], P0_max[i], R0_max[i], nzcoef[i], optimal_cutpoint[i])) #
    #     
    #     write_xlsx(results_LASSO, results_path, col_names = TRUE, format_headers = FALSE)
    #     
    #   }
    # }
    # 
    # 
    #   
    #    # OVERSAMPLING (SMOTE) ----------------------------------------------------------
    #   
    #     results_path <-  paste0(BASE,"Lasso_SMOTE.xlsx")
    #   
    #     set.seed(SEED[j])
    #     df_for_SMOTE <- as.data.frame(cbind(df_train_out,df_train_pred))
    #     df_for_SMOTE$V1 <- factor(df_for_SMOTE$V1)
    #     oversampling <- SMOTE(form = V1~., data=df_for_SMOTE, perc.over = 200,perc.under=150, k=5)
    #     
    #     df_train_out_SMOTE <- oversampling[,1]
    #     df_train_pred_SMOTE <- oversampling[,-c(1)]
    #     cols_to_factor <- colnames(df_train_pred_SMOTE[,c(1:4)])
    #     df_train_pred_SMOTE[,cols_to_factor] <- lapply(df_train_pred_SMOTE[,cols_to_factor],as.factor)
    #     df_train_pred_SMOTE <- apply(as.matrix(df_train_pred_SMOTE), 2, as.numeric)
    #     
    #     for(k in 1:5){
    #       df_train_pred_SMOTE[,k] <- round(df_train_pred_SMOTE[,k], 0)
    #     }
    #     
    #     # Cross Validation
    #     set.seed(SEED[j])
    #     CVFOLD = 4
    #     # ADDED: create folds AND 'foldid' cv.glmnet parameter: to use stratified k folds
    #     flds <- createFolds(df_train_out_SMOTE, k = CVFOLD, list = TRUE, returnTrain = FALSE)
    #     foldids = rep(1,length(df_train_out_SMOTE))
    #     foldids[flds$Fold2] = 2
    #     foldids[flds$Fold3] = 3
    #     foldids[flds$Fold4] = 4
    #     #foldids[flds$Fold5] = 5
    #     #df_train_out_SMOTE <- matrix(df_train_out_SMOTE)
    #     
    #     
    #     set.seed(SEED[j])
    #     cvlasso_SMOTE<-cv.glmnet(df_train_pred_SMOTE,df_train_out_SMOTE,
    #                              alpha=1, family="binomial",
    #                              type.measure="auc", nlambda = 100,
    #                              standardize=FALSE, foldid = foldids, pmax = pmaxs_value[i]) # , weights = pesos2
    #     
    #     #plotcvlasso_SMOTE)
    #     
    #     #coefficients, weighted data
    #     myCoefs <- coef(cvlasso_SMOTE, s="lambda.min");
    #     myCoefs[which(myCoefs != 0 ) ]; #coefficients: intercept included
    #     myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ];  #feature names: intercept included
    #     ## Asseble into a data.frame
    #     myResults <- data.frame(
    #       features = myCoefs@Dimnames[[1]][ which(myCoefs != 0 ) ], #intercept included
    #       coefs    = myCoefs[ which(myCoefs != 0 ) ]  #intercept included
    #     )
    #     print(myResults)
    #     print(paste0("number_coefs = ", nrow(myResults)))
    #     
    #     # NEW: check performance
    #     lmin <- cvlasso_SMOTE[["lambda.min"]]
    #     index_lmin <- which(cvlasso_SMOTE[["lambda"]] == lmin)
    #     measure_lmin[i] <- cvlasso_SMOTE[["cvm"]][index_lmin]
    #     std_measure_lmin[i] <- cvlasso_SMOTE[["cvsd"]][index_lmin]
    #     cvup[i] <- cvlasso_SMOTE[["cvup"]][index_lmin]
    #     cvlo[i] <- cvlasso_SMOTE[["cvlo"]][index_lmin]
    #     print(paste0("auclmin = ", round(measure_lmin,3), " std = ", round(std_measure_lmin,3)))
    #     
    #     #Predictions after CV: use predict() #1) lambda.min
    #     #------ With OCP -----#
    #     response_predictions_train_lmin_SMOTE <- data.matrix(predict(cvlasso_SMOTE$glmnet.fit, newx = df_train_pred_SMOTE, type = "response", s=cvlasso_SMOTE$lambda.min))
    #     response_predictions_test_lmin_SMOTE <- data.matrix(predict(cvlasso_SMOTE$glmnet.fit, newx = df_test_pred, type = "response", s=cvlasso_SMOTE$lambda.min))
    #     
    #     roc_train_SMOTE <- pROC::roc(response=df_train_out_SMOTE, predictor=response_predictions_train_lmin_SMOTE[,1])
    #     #plotroc_train_SMOTE)
    #     roc_train_SMOTE["auc"]
    #     
    #     # Finding OCP:
    #     metric_to_optimize <- "F1_score"  # metric_to_optimize <- "Recall" #metric_to_optimize <- "ppv" #metric_to_optimize <- "roc01" #metric_to_optimize <- "cohens_kappa"
    #     
    #     set.seed(SEED[j])
    #     cutpoint_SMOTE <- cutpointr(x = response_predictions_train_lmin_SMOTE, class = df_train_out_SMOTE, use_midpoints = TRUE,
    #                                 maximize_boot_metric = metric_to_optimize) #, boot_runs = 200
    #     
    #     ##plotcutpoint_SMOTE)
    #     
    #     optimal_cutpoint[i] <- cutpoint_SMOTE$optimal_cutpoint
    #     new_predictions_train_lmin_SMOTE <- dichotomize(response_predictions_train_lmin_SMOTE, cutpoint_SMOTE$optimal_cutpoint)
    #     OCP_confusion_matrix_train_lmin_SMOTE <- confusionMatrix(reference =df_train_out_SMOTE, data = as.factor(new_predictions_train_lmin_SMOTE), mode = "everything", positive = "1")
    #     OCP_confusion_matrix_train_lmin_SMOTE
    #     
    #     # new_predictions_test_lmin_SMOTE <- dichotomize(response_predictions_test_lmin_SMOTE, cutpoint_SMOTE$optimal_cutpoint)
    #     # OCP_confusion_matrix_test_lmin_SMOTE <- confusionMatrix(reference = factor(df_test_out), data = as.factor(new_predictions_test_lmin_SMOTE), mode = "everything", positive = "1")
    #     # OCP_confusion_matrix_test_lmin_SMOTE
    #     
    #     roc_train_auc[i] <- roc_train_SMOTE[["auc"]][1]
    #     F_max[i] <- OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["F1"]]
    #     P_max[i] <- OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["Precision"]]
    #     R_max[i] <- OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["Recall"]]
    #     
    #     OCP_confusion_matrix_train_lmin_SMOTE_negative <- confusionMatrix(reference =df_train_out_SMOTE, data = as.factor(new_predictions_train_lmin_SMOTE), mode = "everything", positive = "0")
    #     OCP_confusion_matrix_train_lmin_SMOTE_negative
    #     
    #     F0_max[i] <- OCP_confusion_matrix_train_lmin_SMOTE_negative[["byClass"]][["F1"]]
    #     P0_max[i] <- OCP_confusion_matrix_train_lmin_SMOTE_negative[["byClass"]][["Precision"]]
    #     R0_max[i] <- OCP_confusion_matrix_train_lmin_SMOTE_negative[["byClass"]][["Recall"]]
    #     
    #     signames <- myCoefs@Dimnames[[1]][which(myCoefs != 0 ) ]
    #     nzcoef[i] <- length(myCoefs[which(myCoefs != 0 ) ])
    #     
    #     
    #     results_LASSO <- rbind(results_LASSO, data.frame(pmaxs_value[i], SEED[j], "LASSO", cvlo[i], cvup[i], measure_lmin[i], std_measure_lmin[i],
    #                                                      + roc_train_auc[i], F_max[i], P_max[i], R_max[i], F0_max[i], P0_max[i], R0_max[i], nzcoef[i], optimal_cutpoint[i])) #
    #     
    #     write_xlsx(results_LASSO, results_path, col_names = TRUE, format_headers = FALSE)
    #   }
    # }
    
    
    
    #------------------------------------------------------------------ PRIORITY LASSO----------------------------------------------------------
    
    #     # AUX Priority LASSO
    #     results_path <-  paste0(BASE,"PreBiosy_PLI_nooversampling.xlsx")
    #     
    #     # get indexes
    #     index_radiomics <- which(grepl("i_" , colnames(df_train_pred)))
    #     index_clinical <- which(grepl("Clin_" , colnames(df_train_pred)))
    #     
    #     # define blocks
    #     bp1_radiomics <- index_radiomics[1]:index_radiomics[length(index_radiomics)]
    #     bp2_clinical <- index_clinical[1]:index_clinical[length(index_clinical)]
    #     
    #     permutations <- perm(list(bp1_radiomics,bp2_clinical))
    #     names_permutations <- perm(list("Radiomics", "PreBIOPSY"))
    #     
    
    #     # Cross Validation (no resampling)
    #     set.seed(SEED[j])
    #     CVFOLD = 4
    #     # ADDED: create folds AND 'foldid' cv.glmnet parameter: to use stratified k folds
    #     flds <- createFolds(y_train_num, k = CVFOLD, list = TRUE, returnTrain = FALSE)
    #     foldids = rep(1,length(y_train_num))
    #     
    #     foldids[flds$Fold2] = 2
    #     foldids[flds$Fold3] = 3
    #     foldids[flds$Fold4] = 4
    #     
    #     #define outcome and predictors
    #     df_train_out <- apply(as.matrix(y_train_num), 2, as.numeric)
    #     df_train_pred <- apply(as.matrix(x_train_f), 2, as.numeric)
    #     
    #     df_test_out <- apply(as.matrix(y_test_num), 2, as.numeric)
    #     df_test_pred <- apply(as.matrix(x_test_f), 2, as.numeric) # matrix with the predictors without label
    #     
    #     
    #     
    #     for (m in 1:dim(permutations)[1]) { #dim(permutations)[1]
    #       print(m)
    #       perm_i <- permutations[m,]
    #       #print(perm_i)
    #       
    #       set.seed(SEED[j])
    #       assign(paste0("pl", m), prioritylasso(X = df_train_pred, Y = df_train_out, family = "binomial",
    #                                             type.measure="auc", nfolds = 4, foldid = foldids,
    #                                             blocks = perm_i, standardize = FALSE, cvoffset=TRUE, cvoffsetnfolds = 4,
    #                                             lambda.type="lambda.min", max.coef = pmaxs_value[i,] )) #   max.coef = c(10,10) )
    #       
    #       varNames[m] <- paste0("pl", m) #adds priority lasso's object names to list
    #       
    #       #### AUC CALCULATION AND PLOTTING
    #       d <- get(varNames[m]) # copy the object to an object named 'd'
    #       
    #       d$nzero
    #       n_coeffs1[m] <- d$nzero[[1]]
    #       n_coeffs2[m] <- d$nzero[[2]] #
    #       total_coeffs[m] <- n_coeffs1[m] + n_coeffs2[m]
    #       
    #       sig <- d$coefficients[d$coefficients != 0]
    #       #save coefficients as one string
    #       aux <- data.frame(sig)
    #       df <- tibble::rownames_to_column(aux, "feat")
    #       aux <- as.vector(t(df))
    #       str[m] <- toString(aux)
    #       max_cvm[m] <- round(d[["min.cvm"]][[2]],4)
    #       
    #       
    #       #------ With OCP -----#
    #       # Finding OCP:
    #       metric_to_optimize <- "F1_score"   #metric_to_optimize <- "Recall"   #metric_to_optimize <- "ppv"   #metric_to_optimize <- "roc01"   #metric_to_optimize <- "cohens_kappa"
    #       PL_response_predictions_train_lmin <- predict(d, df_train_pred, type = "response")
    #       
    #       set.seed(SEED[j])
    #       PL_cutpoint <- cutpointr(x = PL_response_predictions_train_lmin, class = y_train_num,
    #                                use_midpoints = TRUE,
    #                                maximize_boot_metric = metric_to_optimize ) # boot_runs = 200
    #       
    #       optimal_cutpoint[m] <- PL_cutpoint$optimal_cutpoint
    #       
    #       new_predictions_train_lmin <- dichotomize(PL_response_predictions_train_lmin, PL_cutpoint$optimal_cutpoint)
    #       OCP_confusion_matrix_train_lmin <- confusionMatrix(reference =y_train_num, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "1")
    #       print(OCP_confusion_matrix_train_lmin)
    #       
    #       pl_roc_train <- pROC::roc(response=y_train_num, predictor=PL_response_predictions_train_lmin[,1])
    #       #plotpl_roc_train)
    #       pl_roc_train["auc"]
    #       
    #       max_F[m] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["F1"]],4)
    #       max_P[m] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["Precision"]],4)
    #       max_R[m] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["Recall"]],4)
    #       
    #       auc_pl_roc_train[m] <- pl_roc_train[["auc"]]
    #       
    #       OCP_confusion_matrix_train_lmin_negative <- confusionMatrix(reference =y_train_num, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "0")
    #       print(OCP_confusion_matrix_train_lmin_negative)
    #       
    #       max_F0[m] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["F1"]],4)
    #       max_P0[m] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Precision"]],4)
    #       max_R0[m] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Recall"]],4)
    #       
    #       
    #       permutations_matrix <- matrix(unlist(names_permutations), nrow = dim(names_permutations)[1], ncol = dim(names_permutations)[2])
    #       
    #       results_PLI <- rbind(results_PLI, data.frame(pmaxs_value[i,1], pmaxs_value[i,2], SEED[j], "PL-I", varNames[m], permutations_matrix[m], 
    #                                                    permutations_matrix[2], n_coeffs1[m],n_coeffs2[m], total_coeffs[m],
    #                                                    str[m], max_cvm[m], auc_pl_roc_train[m], max_F[m], max_P[m], max_R[m], 
    #                                                    max_F0[m], max_P0[m], max_R0[m], optimal_cutpoint[m]))
    #       
    #       
    #       write_xlsx(results_PLI, results_path, col_names = TRUE, format_headers = FALSE)
    #       
    #     }
    #     
    #   }
    #   
    # }
    
    
    
    # OVERSAMPLING (SMOTE) 2 ----------------------------------------------------------
    
    #     results_path <- 
    #       paste0(BASE,"PreBiopsy_PLI_oversamplingSMOTE.xlsx")
    #     
    #     # AUX Priority LASSO
    #     # get indexes
    #     index_radiomics <- which(grepl("i_" , colnames(df_train_pred)))
    #     index_clinical <- which(grepl("Clin_" , colnames(df_train_pred)))
    #     
    #     # define blocks
    #     bp1_radiomics <- index_radiomics[1]:index_radiomics[length(index_radiomics)]
    #     bp2_clinical <- index_clinical[1]:index_clinical[length(index_clinical)]
    #     
    #     permutations <- perm(list(bp1_radiomics,bp2_clinical))
    #     names_permutations <- perm(list("Radiomics", "PreBIOPSY"))
    #     
    #     # SMOTE
    #     set.seed(SEED[j])
    #     df_for_SMOTE <- as.data.frame(cbind(df_train_out,df_train_pred))
    #     df_for_SMOTE$V1 <- factor(df_for_SMOTE$V1)
    #     oversampling <- SMOTE(form = V1~., data=df_for_SMOTE, perc.over = 200,perc.under=150, k=5)
    #     
    #     df_train_out_SMOTE <- oversampling[,1]
    #     df_train_pred_SMOTE <- oversampling[,-c(1)]
    #     cols_to_factor <- colnames(df_train_pred_SMOTE[,c(1:4)])
    #     df_train_pred_SMOTE[,cols_to_factor] <- lapply(df_train_pred_SMOTE[,cols_to_factor],as.factor)
    #     df_train_pred_SMOTE <- apply(as.matrix(df_train_pred_SMOTE), 2, as.numeric)
    #     
    #     for(p in 1:5){
    #       df_train_pred_SMOTE[,p] <- round(df_train_pred_SMOTE[,p],0)
    #     }
    #     
    #     # Cross Validation SMOTE
    #     set.seed(SEED[j])
    #     CVFOLD = 4
    #     # ADDED: create folds AND 'foldid' cv.glmnet parameter: to use stratified k folds
    #     flds <- createFolds(df_train_out_SMOTE, k = CVFOLD, list = TRUE, returnTrain = FALSE)
    #     foldids = rep(1,length(df_train_out_SMOTE))
    #     foldids[flds$Fold2] = 2
    #     foldids[flds$Fold3] = 3
    #     foldids[flds$Fold4] = 4
    #     
    #     #SMOTE
    #     for (n in 1:dim(permutations)[1]) { #dim(permutations)[1]
    #       print(n)
    #       perm_i <- permutations[n,]
    #       #print(perm_i)
    #       
    #       set.seed(SEED[j])
    #       assign(paste0("pl", n), prioritylasso(X = df_train_pred_SMOTE, Y = df_train_out_SMOTE, family = "binomial",
    #                                             type.measure="auc", nfolds = 4, foldid = foldids,
    #                                             blocks = perm_i, standardize = FALSE, cvoffset=TRUE, cvoffsetnfolds = 4,
    #                                             lambda.type="lambda.min", max.coef = pmaxs_value[i,] )) #   max.coef = c(10,10) )
    #       
    #       varNames[n] <- paste0("pl", n) #adds priority lasso's object names to list
    #       
    #       #### AUC CALCULATION AND PLOTTING
    #       d <- get(varNames[n]) # copy the object to an object named 'd'
    #       
    #       d$nzero
    #       n_coeffs1[n] <- d$nzero[[1]]
    #       n_coeffs2[n] <- d$nzero[[2]] #
    #       total_coeffs[n] <- n_coeffs1[n] + n_coeffs2[n]
    #       sig <- d$coefficients[d$coefficients != 0]
    #       
    #       #save coefficients as one string
    #       aux <- data.frame(sig)
    #       df <- tibble::rownames_to_column(aux, "feat")
    #       aux <- as.vector(t(df))
    #       str[n] <- toString(aux)
    #       max_cvm[n] <- round(d[["min.cvm"]][[2]],4)
    #       
    #       
    #       #------ With OCP -----#
    #       # Finding OCP:
    #       metric_to_optimize <- "F1_score"   #metric_to_optimize <- "Recall"   #metric_to_optimize <- "ppv"   #metric_to_optimize <- "roc01"   #metric_to_optimize <- "cohens_kappa"
    #       PL_response_predictions_train_lmin_SMOTE <- predict(d, df_train_pred_SMOTE, type = "response")
    #       
    #       set.seed(SEED[j])
    #       PL_cutpoint_SMOTE <- cutpointr(x = PL_response_predictions_train_lmin_SMOTE,
    #                                      class = df_train_out_SMOTE,
    #                                      use_midpoints = TRUE,
    #                                      maximize_boot_metric = metric_to_optimize, #boot_runs = 200
    #                                      pos_class = 1)
    #       
    #       optimal_cutpoint[n] <- PL_cutpoint_SMOTE$optimal_cutpoint
    #       new_predictions_train_lmin_SMOTE <- dichotomize(PL_response_predictions_train_lmin_SMOTE, PL_cutpoint_SMOTE$optimal_cutpoint)
    #       OCP_confusion_matrix_train_lmin_SMOTE <- confusionMatrix(reference =df_train_out_SMOTE, data = as.factor(new_predictions_train_lmin_SMOTE), mode = "everything", positive = "1")
    #       #print(OCP_confusion_matrix_train_lmin_SMOTE)
    #       
    #       max_F[n] <- round(OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["F1"]],4)
    #       max_P[n] <- round(OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["Precision"]],4)
    #       max_R[n] <- round(OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["Recall"]],4)
    #       
    #       OCP_confusion_matrix_train_lmin_SMOTE_negative <- confusionMatrix(reference =df_train_out_SMOTE, data = as.factor(new_predictions_train_lmin_SMOTE), mode = "everything", positive = "0")
    #       #print(OCP_confusion_matrix_train_lmin_SMOTE_negative)
    #       
    #       max_F0[n] <- round(OCP_confusion_matrix_train_lmin_SMOTE_negative[["byClass"]][["F1"]],4)
    #       max_P0[n] <- round(OCP_confusion_matrix_train_lmin_SMOTE_negative[["byClass"]][["Precision"]],4)
    #       max_R0[n] <- round(OCP_confusion_matrix_train_lmin_SMOTE_negative[["byClass"]][["Recall"]],4)
    #       
    #       # PL_response_predictions_test_lmin_SMOTE <- predict(d, df_test_pred, type = "response")
    #       # 
    #       # new_predictions_test_lmin_SMOTE <- dichotomize(PL_response_predictions_test_lmin_SMOTE, PL_cutpoint_SMOTE$optimal_cutpoint)
    #       # OCP_confusion_matrix_test_lmin_SMOTE <- confusionMatrix(reference = as.factor(df_test_out), data = as.factor(new_predictions_test_lmin_SMOTE), mode = "everything", positive = "1")
    #       # print(OCP_confusion_matrix_test_lmin_SMOTE)
    #       
    #       
    #       #score_train_SMOTE <- as.vector(as.matrix(df_train_pred_SMOTE[ ,names(sig)]) %*% sig)
    #       # score_test <- as.vector(as.matrix(df_test_pred[ ,names(sig)]) %*% sig)
    #       # 
    #       # pl_roc_train_SMOTE <- pROC::roc(as.factor(df_train_out_SMOTE), score_train_SMOTE)
    #       pl_roc_train_SMOTE <- pROC::roc(df_train_out_SMOTE, PL_response_predictions_train_lmin_SMOTE[,1])
    #       #plotpl_roc_train_SMOTE)
    #       #print(pl_roc_train_SMOTE)
    #       auc_pl_roc_train[n] <- pl_roc_train_SMOTE[["auc"]][1]
    #       
    #       #   str2 <- cbind(max_cvm,
    #       #                 pl_roc_train_SMOTE[["auc"]],
    #       #                 OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["F1"]], 
    #       #                 OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["Precision"]], 
    #       #                 OCP_confusion_matrix_train_lmin_SMOTE[["byClass"]][["Recall"]])
    #       #   print(str2)
    #       #   
    #       
    #       permutations_matrix <- matrix(unlist(names_permutations), nrow = dim(names_permutations)[1], ncol = dim(names_permutations)[2])
    #       
    #       results_PLI <- rbind(results_PLI, data.frame(pmaxs_value[i,1], pmaxs_value[i,2], SEED[j], "PL-I", varNames[n], permutations_matrix[n,1],
    #                                                    permutations_matrix[n,2], n_coeffs1[n],n_coeffs2[n], total_coeffs[n],
    #                                                    str[n], max_cvm[n], auc_pl_roc_train[n], max_F[n], max_P[n], max_R[n],
    #                                                    max_F0[n], max_P0[n], max_R0[n], optimal_cutpoint[n]))
    #       
    #       
    #       write_xlsx(results_PLI, results_path, col_names = TRUE, format_headers = FALSE)
    #       
    #     }
    #     
    #   }
    #   
    # }
    
    
    
    #--------------------------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------ PRIORITY ----------------------------------------------------------
    #--------------------------------------------------------------------LASSO-------------------------------------------------------------
    #----------------------------------------------------------------------II - 1---------------------------------------------------------------
    
    #     #define outcome and predictors
    #     results_path <- 
    #       paste0(BASE,"PreBiopsy_PLII_nooversampling.xlsx")
    #     
    #     df_train_out <- apply(as.matrix(y_train_num), 2, as.numeric)
    #     df_train_pred <- apply(as.matrix(x_train_f), 2, as.numeric)
    #     
    #     df_test_out <- apply(as.matrix(y_test_num), 2, as.numeric)
    #     df_test_pred <- apply(as.matrix(x_test_f), 2, as.numeric) # matrix with the predictors without label
    #     
    #     # Cross Validation
    #     set.seed(SEED[j])
    #     CVFOLD = 4
    #     # ADDED: create folds AND 'foldid' cv.glmnet parameter: to use stratified k folds
    #     flds <- createFolds(y_train_num, k = CVFOLD, list = TRUE, returnTrain = FALSE)
    #     foldids = rep(1,length(y_train_num))
    #     foldids[flds$Fold2] = 2
    #     foldids[flds$Fold3] = 3
    #     foldids[flds$Fold4] = 4
    #     #foldids[flds$Fold5] = 5
    #     
    #     #AUX Priority-lasso
    #     # define blocks
    #     #df_train_pred <- df_train_pred[,-c(1:5)]
    #     
    #     index_T2 <- which(grepl("T2i" , colnames(df_train_pred)))
    #     index_b1000 <- which(grepl("b1000i" , colnames(df_train_pred)))
    #     index_ADC <- which(grepl("ADCi" , colnames(df_train_pred)))
    #     index_clin <- which(grepl("Clin_" , colnames(df_train_pred)))
    #     
    #     bp1<-index_T2[1]:index_T2[length(index_T2)]
    #     bp2<-index_b1000[1]:index_b1000[length(index_b1000)]
    #     bp3<-index_ADC[1]:index_ADC[length(index_ADC)]
    #     bp4<-index_clin[1]:index_clin[length(index_clin)]
    #     
    #     permutations <- perm(list(bp1,bp2,bp3,bp4))
    #     names_permutations <- perm(list("T2","b1000", "ADC", "Clinical")) 
    #     
    #     # #TO WRITE CSV
    #     # varNames <- vector() #will store the objects names of the priority lasso's
    #     # aucsTrain <- vector()
    #     # stdTrain <- vector()
    #     # aucsTest <- vector()
    #     # stdTest <- vector()
    #     # cisTrain_l <- vector()
    #     # cisTrain_u <- vector()
    #     # cisTest_l <- vector()
    #     # cisTest_u <- vector()
    #     # n_coeffs1 <- vector()
    #     # n_coeffs2 <- vector()
    #     # n_coeffs3 <- vector()
    #     # n_coeffs4 <- vector()
    #     # str <- vector()
    #     # max_cvm <- vector()
    #     # max_F <- vector()
    #     # max_P <- vector()
    #     # max_R <- vector()
    #     # max_F0 <- vector()
    #     # max_P0 <- vector()
    #     # max_R0 <- vector()
    #     # pl_roc_train_auc <- vector()
    #     
    #     
    #     
    #     
    #     for (p in 1:dim(permutations)[1]) { #dim(permutations)[1]
    #       print(p)
    #       perm_i <- permutations[p,]
    #       #print(perm_i)
    #       
    #       set.seed(SEED[j])
    #       assign(paste0("pl", p), prioritylasso(X = df_train_pred, Y = df_train_out, family = "binomial",
    #                                             type.measure="auc", nfolds = 4, foldid = foldids,
    #                                             blocks = perm_i, standardize = FALSE, cvoffset=TRUE, cvoffsetnfolds = 4,
    #                                             lambda.type="lambda.min", max.coef = pmaxs_value[i,] ))  #max.coef = c(3,3,3,3) ))
    #       
    #       varNames[p] <- paste0("pl", p) #adds priority lasso's object names to list
    #       
    #       #### AUC CALCULATION AND PLOTTING
    #       d <- get(varNames[p]) # copy the object to an object named 'd'
    #       
    #       d$nzero
    #       n_coeffs1[p] <- d$nzero[[1]]
    #       n_coeffs2[p] <- d$nzero[[2]]
    #       n_coeffs3[p] <- d$nzero[[3]]
    #       n_coeffs4[p] <- d$nzero[[4]]
    #       total_coeffs[p] <- n_coeffs1[p] + n_coeffs2[p] + n_coeffs3[p] + n_coeffs4[p]
    #       
    #       sig <- d$coefficients[d$coefficients != 0]
    #       #save coefficients as one string
    #       aux <- data.frame(sig)
    #       df <- tibble::rownames_to_column(aux, "feat")
    #       aux <- as.vector(t(df))
    #       str[p] <- toString(aux)
    #       max_cvm[p] <- round(d[["min.cvm"]][[4]],4)
    #       
    #       
    #       #------ With OCP -----#
    #       #Finding OCP:
    #       metric_to_optimize <- "F1_score"   #metric_to_optimize <- "Recall"   #metric_to_optimize <- "ppv"   #metric_to_optimize <- "roc01"   #metric_to_optimize <- "cohens_kappa"
    #       PL_response_predictions_train_lmin <- predict(d, df_train_pred, type = "response")
    #       
    #       pl_roc_train <- pROC::roc(df_train_out, PL_response_predictions_train_lmin[,1])
    #       ##plotpl_roc_train)
    #       #print(pl_roc_train)
    #       pl_roc_train_auc[p] <- pl_roc_train[["auc"]][1]
    #       
    #       set.seed(SEED[j])
    #       PL_cutpoint <- cutpointr(x = PL_response_predictions_train_lmin,
    #                                class = y_train_num,
    #                                use_midpoints = TRUE,
    #                                maximize_boot_metric = metric_to_optimize, 
    #                                pos_class = 1) #, boot_runs = 200
    #       
    #       optimal_cutpoint[p] <- PL_cutpoint$optimal_cutpoint
    #       new_predictions_train_lmin <- dichotomize(PL_response_predictions_train_lmin, PL_cutpoint$optimal_cutpoint)
    #       OCP_confusion_matrix_train_lmin <- confusionMatrix(reference =y_train_num, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "1")
    #       print(OCP_confusion_matrix_train_lmin)
    #       
    #       max_F[p] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["F1"]],4)
    #       max_P[p] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["Precision"]],4)
    #       max_R[p] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["Recall"]],4)
    #       
    #       OCP_confusion_matrix_train_lmin_negative <- confusionMatrix(reference =y_train_num, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "0")
    #       print(OCP_confusion_matrix_train_lmin_negative)
    #       
    #       max_F0[p] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["F1"]],4)
    #       max_P0[p] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Precision"]],4)
    #       max_R0[p] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Recall"]],4)
    #       
    #       
    #       
    #       
    #       
    #       
    #       
    #       
    #       
    #       permutations_matrix <- matrix(unlist(names_permutations), nrow = dim(names_permutations)[1], ncol = dim(names_permutations)[2])
    #       
    #       results_PLII <- rbind(results_PLII, data.frame(pmaxs_value[i,1], pmaxs_value[i,2], pmaxs_value[i,3], pmaxs_value[i,4], SEED[j],
    #                                                      "PL-II", varNames[p], permutations_matrix[p,1],
    #                                                     permutations_matrix[p,2], permutations_matrix[p,3], permutations_matrix[p,4],
    #                                                     n_coeffs1[p],n_coeffs2[p], n_coeffs3[p], n_coeffs4[p], total_coeffs[p],
    #                                                     str[p], max_cvm[p], pl_roc_train_auc[p], max_F[p], max_P[p], max_R[p],
    #                                                     max_F0[p], max_P0[p], max_R0[p], optimal_cutpoint[p]))
    #       
    #       
    #       write_xlsx(results_PLII, results_path, col_names = TRUE, format_headers = FALSE)
    #       
    #     }
    #     
    #   }
    #   
    # }
    
    
    
    ###############################################################################
    # OVERSAMPLING (SMOTE) 1 ----------------------------------------------------------
    
    results_path <- paste0(BASE,"PreBiopsy_PLII_SMOTEoversampling.xlsx")
    
    #SMOTEoversampling2121
    
    set.seed(SEED[j])
    df_for_SMOTE <- as.data.frame(cbind(df_train_out,df_train_pred))
    df_for_SMOTE$V1 <- factor(df_for_SMOTE$V1)
    oversampling <- SMOTE(form = V1~., data=df_for_SMOTE, perc.over = 200,perc.under=150, k=5)
    
    df_train_out_SMOTE <- oversampling[,1]
    df_train_pred_SMOTE <- oversampling[,-c(1)]
    
    cols_to_factor <- colnames(df_train_pred_SMOTE[,c(1:4)])
    df_train_pred_SMOTE[,cols_to_factor] <- lapply(df_train_pred_SMOTE[,cols_to_factor],as.factor)
    df_train_pred_SMOTE <- apply(as.matrix(df_train_pred_SMOTE), 2, as.numeric)
    
    for(u in 1:5){
      df_train_pred_SMOTE[,u]<- round(df_train_pred_SMOTE[,u],0)
    }
    
    #AUX Priority-lasso
    # define blocks
    #df_train_pred <- df_train_pred[,-c(1:5)]
    
    index_T2 <- which(grepl("T2ii" , colnames(df_train_pred_SMOTE)))
    index_b1000 <- which(grepl("b1000ii" , colnames(df_train_pred_SMOTE)))
    index_ADC <- which(grepl("ADCii" , colnames(df_train_pred_SMOTE)))
    index_clin <- which(grepl("Clin_" , colnames(df_train_pred_SMOTE)))
    
    bp1<-index_T2[1]:index_T2[length(index_T2)]
    bp2<-index_b1000[1]:index_b1000[length(index_b1000)]
    bp3<-index_ADC[1]:index_ADC[length(index_ADC)]
    bp4<-index_clin[1]:index_clin[length(index_clin)]
    
    permutations <- perm(list(bp1,bp2,bp3,bp4))
    names_permutations <- perm(list("T2","b1000", "ADC", "Clinical"))
    
    # Cross Validation SMOTE
    set.seed(SEED[j])
    CVFOLD = 4
    # ADDED: create folds AND 'foldid' cv.glmnet parameter: to use stratified k folds
    flds <- createFolds(df_train_out_SMOTE, k = CVFOLD, list = TRUE, returnTrain = FALSE)
    foldids = rep(1,length(df_train_out_SMOTE))
    foldids[flds$Fold2] = 2
    foldids[flds$Fold3] = 3
    foldids[flds$Fold4] = 4
    #foldids[flds$Fold5] = 5
    #df_train_out_SMOTE <- matrix(df_train_out_SMOTE)
    
    
    # #TO WRITE CSV
    # varNames <- vector() #will store the objects names of the priority lasso's
    # aucsTrain <- vector()
    # stdTrain <- vector()
    # aucsTest <- vector()
    # stdTest <- vector()
    # cisTrain_l <- vector()
    # cisTrain_u <- vector()
    # cisTest_l <- vector()
    # cisTest_u <- vector()
    # n_coeffs1 <- vector()
    # n_coeffs2 <- vector()
    # n_coeffs3 <- vector()
    # n_coeffs4 <- vector()
    # str <- vector()
    # max_cvm <- vector()
    # max_F <- vector()
    # max_P <- vector()
    # max_R <- vector()
    # max_F0 <- vector()
    # max_P0 <- vector()
    # max_R0 <- vector()
    # pl_roc_train_auc <- vector()
    
    
    for (p in 1:dim(permutations)[1]) { #dim(permutations)[1]
      print(p)
      perm_i <- permutations[p,]
      #print(perm_i)
      
      set.seed(SEED[j])
      assign(paste0("pl", p), prioritylasso(X = df_train_pred_SMOTE, Y = df_train_out_SMOTE, family = "binomial",
                                            type.measure="auc", nfolds = 4, foldid = foldids,
                                            blocks = perm_i, standardize = FALSE, cvoffset=TRUE, cvoffsetnfolds = 4,
                                            lambda.type="lambda.min", max.coef = pmaxs_value[i,] ))  #pmax = 10 ) ) #  max.coef = c(3,3,3,3)
      
      varNames[p] <- paste0("pl", p) #adds priority lasso's object names to list
      
      #### AUC CALCULATION AND PLOTTING
      d <- get(varNames[p]) # copy the object to an object named 'd'
      
      d$nzero
      n_coeffs1[p] <- d$nzero[[1]]
      n_coeffs2[p] <- d$nzero[[2]]
      n_coeffs3[p] <- d$nzero[[3]]
      n_coeffs4[p] <- d$nzero[[4]]
      total_coeffs[p] <- n_coeffs1[p] + n_coeffs2[p] + n_coeffs3[p] + n_coeffs4[p]
      
      sig <- d$coefficients[d$coefficients != 0]
      #save coefficients as one string
      aux <- data.frame(sig)
      df <- tibble::rownames_to_column(aux, "feat")
      aux <- as.vector(t(df))
      str[p] <- toString(aux)
      max_cvm[p] <- round(d[["min.cvm"]][[4]],4)
      
      
      #------ With OCP -----#
      #Finding OCP:
      metric_to_optimize <- "F1_score"   #metric_to_optimize <- "Recall"   #metric_to_optimize <- "ppv"   #metric_to_optimize <- "roc01"   #metric_to_optimize <- "cohens_kappa"
      PL_response_predictions_train_lmin <- predict(d, df_train_pred_SMOTE, type = "response")
      
      pl_roc_train <- pROC::roc(df_train_out_SMOTE, PL_response_predictions_train_lmin[,1])
      ##plotpl_roc_train)
      #print(pl_roc_train)
      pl_roc_train_auc[p] <- pl_roc_train[["auc"]][1]
      
      set.seed(SEED[j])
      PL_cutpoint <- cutpointr(x = PL_response_predictions_train_lmin, class = df_train_out_SMOTE,use_midpoints = TRUE, 
                               maximize_boot_metric = metric_to_optimize, pos_class = 1) #boot_runs = 200)
      optimal_cutpoint[p] <- PL_cutpoint$optimal_cutpoint
      
      new_predictions_train_lmin <- dichotomize(PL_response_predictions_train_lmin, PL_cutpoint$optimal_cutpoint)
      OCP_confusion_matrix_train_lmin <- confusionMatrix(reference =df_train_out_SMOTE, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "1")
      #print(OCP_confusion_matrix_train_lmin)
      
      max_F[p] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["F1"]],4)
      max_P[p] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["Precision"]],4)
      max_R[p] <- round(OCP_confusion_matrix_train_lmin[["byClass"]][["Recall"]],4)
      
      OCP_confusion_matrix_train_lmin_negative <- confusionMatrix(reference =df_train_out_SMOTE, data = as.factor(new_predictions_train_lmin), mode = "everything", positive = "0")
      #print(OCP_confusion_matrix_train_lmin_negative)
      
      max_F0[p] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["F1"]],4)
      max_P0[p] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Precision"]],4)
      max_R0[p] <- round(OCP_confusion_matrix_train_lmin_negative[["byClass"]][["Recall"]],4)
      
      #PL_response_predictions_test_lmin <- predict(d, x_test_f, type = "response")
      
      #new_predictions_test_lmin <- dichotomize(PL_response_predictions_test_lmin, PL_cutpoint$optimal_cutpoint)
      #OCP_confusion_matrix_test_lmin <- confusionMatrix(reference = y_test_num, data = as.factor(new_predictions_test_lmin), mode = "everything", positive = "1")
      #OCP_confusion_matrix_test_lmin
      # 
      
      
      
      
      
      permutations_matrix <- matrix(unlist(names_permutations), nrow = dim(names_permutations)[1], ncol = dim(names_permutations)[2])
      
      results_PLII <- rbind(results_PLII, data.frame(pmaxs_value[i,1], pmaxs_value[i,2], pmaxs_value[i,3], pmaxs_value[i,4],
                                                     SEED[j],
                                                     "PL-II", varNames[p], permutations_matrix[p,1],
                                                     permutations_matrix[p,2], permutations_matrix[p,3], 
                                                     permutations_matrix[p,4],
                                                     n_coeffs1[p],n_coeffs2[p], n_coeffs3[p], n_coeffs4[p], total_coeffs[p],
                                                     str[p], max_cvm[p], pl_roc_train_auc[p], max_F[p], max_P[p], max_R[p],
                                                     max_F0[p], max_P0[p], max_R0[p], optimal_cutpoint[p]))
      
      
      write_xlsx(results_PLII, results_path, col_names = TRUE, format_headers = FALSE)
    }
    
  }
  
}      

