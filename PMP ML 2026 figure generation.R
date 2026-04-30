library(ggplot2)
library(hrbrthemes)
library(lattice)
library("epitools")
library("lmerTest")
library("survival")
library("survminer")
library(data.table)
library(cutpointr)
library(pROC)
library(ggforce)
library(reshape2)
library(extrafont)
library(tidyverse)
library(devtools)
library(readxl)
library(rms)
library(dplyr)
library(readxl)
library(stringr)
library(Seurat)
library(limma)
library(purrr)
library(edgeR)
library(msigdbr)
library(AUCell)
library(speckle)
library(tableone)
library(ggtext)
library(scales)
library(glmnet)
library(logistf)
library(ggh4x)
library(ggsurvfit)
library(colorspace)

options(scipen=4)
fonts()
loadfonts()

.cm = 2.54

setwd("C:/filepath")

#Figure 1B
FData <- read.csv("data_path.csv")
Final <- FData

convert_categorical <- function(x) {
  if (is.character(x) || is.factor(x)) {
    return(as.numeric(as.factor(x)))
  } else {
    return(x)
  }
}

cor_frame <- data.frame(Final$Age, Final$SexMale2, Final$BD, Final$Ulceration, Final$MitoticRate, Final$SiteHN, Final$HistoTwoLevel)
cor_frame <- as.data.frame(lapply(cor_frame, convert_categorical))

cor_labels <- c(
  "Age", "Sex", "Breslow depth", "Ulceration", "Mitotic rate", "Site", "Histo. subtype"
)
cor_labels_additional <- c(
  "Years", "Male", "mm", "Present", "/mm²", "Head/neck", "Nodular"
)

for (label_index in 1:length(cor_labels)) {
  cor_labels[label_index] <- paste0(cor_labels[label_index], "<br><span style='font-family: Montserrat; font-size:6pt;'>*", 
                                    cor_labels_additional[label_index], "*</span>")
}

cor_labels_x <- c(
  "Age", "Sex", "Depth", "Ulc.", "Mit. rate", "Site", "Subtype"
)

cor_matrix <- cor(cor_frame, method = "spearman", use = 'pairwise.complete.obs')
cor_matrix[lower.tri(cor_matrix)] <- NA
melted_matrix <- reshape2::melt(cor_matrix, na.rm = TRUE)

cor_heatmap <- ggplot(data = melted_matrix, aes(Var1, Var2, fill = value)) + 
  geom_tile(width = 1) +
  theme_bw(base_size = 8, base_family = "Montserrat SemiBold") +
  scale_x_discrete(name = "", expand=c(0,0), labels = cor_labels_x) +
  scale_y_discrete(name = "", expand=c(0,0), labels = rev(cor_labels), limits = rev) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.825, 0.79),
        legend.title = element_markdown(family = "Montserrat SemiBold", size = 8, margin = margin(t = 0, r = 0, b = 8, l = 0)),
        legend.text = element_text(size = 6, family = "Montserrat"),
        axis.ticks = element_blank(),
        axis.text.x = element_markdown(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1, margin = margin(t = 3, r = 0, b = 0, l = 0)),
        axis.text.y = element_markdown(color = "black", size = 8, margin = margin(t = 0, r = 4.5, b = 0, l = 0)),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.key.height = unit(10, 'pt'),
        legend.key.width = unit(8, 'pt')) +
  coord_fixed() +
  scale_fill_gradient2(low = "#FFFFFF", mid = "#FF7A7A", high = "#EA6060", 
                       midpoint = 0.5, limit = c(-0.01, 1), space = "Lab", 
                       name = "Spearman's <span style='font-family: Arial; font-size:8pt;'>ρ</span>") +
  geom_text(aes(Var1, Var2, label = sprintf("%.2f", value)), color = "black", family = "Montserrat", size = 8/.pt)
ggsave("COR.png", plot=print(cor_heatmap), width = 9.5/.cm, height = 8.5/.cm, dpi = 900, grDevices::png)

#Full correlation map
Final$OccTradesClerical <- ifelse(Final$Occupation == "Trade/elem. clerical", 1, 0)
Final$OccAdminService <- ifelse(Final$Occupation == "Inter./adv. clerical", 1, 0)
Final$OccProfessional <- ifelse(Final$Occupation == "Professional", 1, 0)

cor_frame_full <- data.frame(Final$Age, Final$SexMale2, Final$BD, Final$Ulceration, Final$MitoticRate, Final$SiteHN, Final$HistoTwoLevel,
                             Final$Regression, Final$PreNaevus, Final$Birthplace, Final$SEIFADecile, Final$RemotenessCat,
                             Final$OccTradesClerical, Final$OccAdminService, Final$OccProfessional, Final$Education, Final$Partner, Final$NumPregnancies,
                             Final$EverSmoker, Final$Sunscreen, Final$SunOccup, Final$SelfSkinChecks, Final$DoctorSkinChecks,
                             Final$PrevMelanoma, Final$PrevSkinCancer, Final$FHxMelanoma,
                             Final$BMI, Final$HeartDisease, Final$RespDisease, Final$KidneyDisease,
                             Final$CorticosteroidLongTerm, Final$NSAIDLongTerm, Final$PsychotropicsLongTerm, Final$StatinsLongTerm
                             )
cor_frame_full$Final.Birthplace <- factor(cor_frame_full$Final.Birthplace, levels = c("Elsewhere", "Aus/NZ"))
cor_frame_full$Final.RemotenessCat <- ifelse(cor_frame_full$Final.RemotenessCat != "Major cities", "Regional/remote", cor_frame_full$Final.RemotenessCat)
cor_frame_full <- as.data.frame(lapply(cor_frame_full, convert_categorical))

cor_labels_full <- c(
  "Age", "Male sex", "Breslow depth", "Ulceration", "Mitotic rate", "Head/neck site", "Nodular subtype",
  "Regression", "Pre-existing naevus", "Born in Aus or NZ", "SES (SEIFA) decile", "Regional or remote", 
  "Trades/clerical occ.", "Admin/service occ.", "Professional occ.", "Finished HS", "Partner at diagnosis", "Pregnancies",
  "Ever smoker", "Daily sunscreen use", "Occ. sun exposure", "Annual self SEs", "Annual doctor SEs",
  "Prev. melanoma", "Prev. skin cancer", "FHx of melanoma",
  "BMI", "Heart disease", "Resp. disease", "Renal disease",
  "LT corticosteroid use", "LT NSAID use", "LT CNS agent use", "LT statin use"
)

cor_matrix_full <- cor(cor_frame_full, method = "spearman", use = 'pairwise.complete.obs')
cor_matrix_full <- cor_matrix_full
cor_matrix_full[lower.tri(cor_matrix_full)] <- NA
melted_matrix_full <- reshape2::melt(cor_matrix_full, na.rm = TRUE)

cor_heatmap_full <- ggplot(data = melted_matrix_full, aes(Var1, Var2, fill = value)) + 
  geom_tile(width = 1) +
  theme_bw(base_size = 6, base_family = "Montserrat") +
  scale_x_discrete(name = "", expand=c(0,0), labels = cor_labels_full) +
  scale_y_discrete(name = "", expand=c(0,0), labels = rev(cor_labels_full), limits = rev) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.825, 0.79),
        legend.title = element_markdown(family = "Montserrat", size = 10, margin = margin(t = 0, r = 0, b = 8, l = 0)),
        legend.text = element_text(size = 8, family = "Montserrat"),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_markdown(size = 6, color = "black", angle = 40, hjust = 1, vjust = 1, margin = margin(t = 2.5, r = 0, b = 0, l = 0)),
        axis.text.y.left = element_markdown(size = 6.4, color = "black", vjust = 0.55, margin = margin(t = 0, r = 3.5, b = 0, l = 0)),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.key.height = unit(12, 'pt'),
        legend.key.width = unit(8, 'pt')) +
  coord_fixed() +
  scale_fill_gradient2(low = "#5070FA", mid = "#FFFFFF", high = "#ED6060", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Spearman's <span style='font-family: Arial; font-size:10pt;'>ρ</span>") +
  geom_text(aes(Var1, Var2, label = sprintf("%.2f", value)), color = "black", family = "Montserrat", size = 4.5/.pt)
ggsave("COR_full.png", plot=print(cor_heatmap_full), width = 18.5/.cm, height = 18/.cm, dpi = 900, grDevices::png)

#Figure 1C, 1D
auc_outcome_vars <- c("Y_Actual", "Y_DM", "Y_MelDeath")
auc_outcome_names <- c("Recurrence", "Distant metastasis", "Melanoma mortality")
auc_outcome_pos <- c(0.85, 0.825, 0.8) # legend xpos
aggregateprobs <- read.csv("XGB.csv")
aggregateprobs_tumour_only <- read.csv("XGB_res.csv")

aggregateprobs$LR_T_Prob <- aggregateprobs_tumour_only$LR_Prob
aggregateprobs$XGB_T_Prob <- aggregateprobs_tumour_only$XGB_Prob

inc_tumour_in_graph <- TRUE
disp_points <- FALSE

for (auc_var_index in 1:length(auc_outcome_vars)) {
  auc_var <- auc_outcome_vars[auc_var_index]
  auc_name <- auc_outcome_names[auc_var_index]
  auc_pos <- auc_outcome_pos[auc_var_index]
  
  train_data <- aggregateprobs[aggregateprobs$Set == "Train",]
  test_data <- aggregateprobs[aggregateprobs$Set == "Test",]
  test_noSLNB <- aggregateprobs[aggregateprobs$SLNBStatus == "Not performed" & aggregateprobs$Set == "Test",]
  test_delayedrec <- aggregateprobs[aggregateprobs$DelayedRec == 0 & aggregateprobs$Set == "Test",]
  test_IBIIA <- aggregateprobs[aggregateprobs$Stage %in% c("IB", "IIA") & aggregateprobs$Set == "Test",]
  
  train_data_tumour_only <- aggregateprobs_tumour_only[aggregateprobs_tumour_only$Set == "Train",]
  test_data_tumour_only <- aggregateprobs_tumour_only[aggregateprobs_tumour_only$Set == "Test",]
  test_noSLNB_tumour_only <- aggregateprobs_tumour_only[aggregateprobs_tumour_only$SLNBStatus == "Not performed" & aggregateprobs_tumour_only$Set == "Test",]
  test_delayedrec_tumour_only <- aggregateprobs_tumour_only[aggregateprobs_tumour_only$DelayedRec == 0 & aggregateprobs_tumour_only$Set == "Test",]
  test_IBIIA_tumour_only <- aggregateprobs_tumour_only[aggregateprobs_tumour_only$Stage %in% c("IB", "IIA") & aggregateprobs_tumour_only$Set == "Test",]
  
  data_sets <- list(train_data, test_data, test_noSLNB, test_delayedrec, test_IBIIA)
  set_names <- c("Train", "Test", "TestNoSLNB", "TestDelayedRec", "TestIBIIA")
  
  data_sets_tumour_only <- list(train_data_tumour_only, test_data_tumour_only, test_noSLNB_tumour_only, test_delayedrec_tumour_only, test_IBIIA_tumour_only)
  
  AUC_vals <- list()
  AUC_frames <- list()
  AUC_CI1 <- list()
  AUC_CI3 <- list()
  data_names <- c("LR", "LRIso", "XGB")
  data_columns <- c("LR_Prob", "LRIso_Prob", "XGB_Prob")
  
  optimal_vals <- list()
  
  for (name_index in 1:length(data_names)) {
    cur_column <- data_columns[name_index]
    data_name <- data_names[name_index]
    
    for (set_index in 1:length(data_sets)) {
      data_set <- data_sets[[set_index]]
      set_name <- set_names[set_index]
      
      data_set_tumour_only <- data_sets_tumour_only[[set_index]]
      
      roc_data_set <- roc(data_set[[auc_var]], data_set[[cur_column]])
      roc_data_set_tumour_only <- roc(data_set_tumour_only[[auc_var]], data_set_tumour_only[[cur_column]])
      
      model_name <- paste0(set_name, data_name)
      
      if (set_name == "Train") {
        youden_index <- roc_data_set$sensitivities + roc_data_set$specificities - 1
        optimal_index <- which.max(youden_index)
        optimal_cutoff <- roc_data_set$thresholds[optimal_index]
        
        youden_index_tumour_only <- roc_data_set_tumour_only$sensitivities + roc_data_set_tumour_only$specificities - 1
        optimal_index_tumour_only <- which.max(youden_index_tumour_only)
        optimal_cutoff_tumour_only <- roc_data_set_tumour_only$thresholds[optimal_index_tumour_only]
      }
      
      all_positives <- data_set[[auc_var]] == 1
      optimal_sensitivity <- sum(data_set[[cur_column]][all_positives] > optimal_cutoff) / sum(all_positives)
      all_negatives <- data_set[[auc_var]] == 0
      optimal_specificity <- sum(data_set[[cur_column]][all_negatives] <= optimal_cutoff) / sum(all_negatives)
      
      all_positives_tumour_only <- data_set_tumour_only[[auc_var]] == 1
      optimal_sensitivity_tumour_only <- sum(data_set_tumour_only[[cur_column]][all_positives_tumour_only] > optimal_cutoff_tumour_only) / sum(all_positives_tumour_only)
      all_negatives_tumour_only <- data_set_tumour_only[[auc_var]] == 0
      optimal_specificity_tumour_only <- sum(data_set_tumour_only[[cur_column]][all_negatives_tumour_only] <= optimal_cutoff_tumour_only) / sum(all_negatives_tumour_only)
      
      AUC_frames[[model_name]] <- data.frame(
        Model = model_name,
        Sensitivity = roc_data_set$sensitivities,
        Specificity = roc_data_set$specificities,
        Threshold = roc_data_set$thresholds
      )
      
      AUC_frames[[paste0(model_name, "_Tumour")]] <- data.frame(
        Model = paste0(model_name, "_Tumour"),
        Sensitivity = roc_data_set_tumour_only$sensitivities,
        Specificity = roc_data_set_tumour_only$specificities,
        Threshold = roc_data_set_tumour_only$thresholds
      )
      
      optimal_vals[[model_name]] <- data.frame(
        Model = model_name,
        Cutoff = optimal_cutoff,
        Sensitivity = optimal_sensitivity,
        Specificity = optimal_specificity
      )
      
      optimal_vals[[paste0(model_name, "_Tumour")]] <- data.frame(
        Model = paste0(model_name, "_Tumour"),
        Cutoff = optimal_cutoff_tumour_only,
        Sensitivity = optimal_sensitivity_tumour_only,
        Specificity = optimal_specificity_tumour_only
      )
      
      AUC_CI1[[model_name]] <- ci(roc_data_set)[1]
      AUC_CI3[[model_name]] <- ci(roc_data_set)[3]
      AUC_vals[[model_name]] <- roc_data_set$auc
      
      AUC_CI1[[paste0(model_name, "_Tumour")]] <- ci(roc_data_set_tumour_only)[1]
      AUC_CI3[[paste0(model_name, "_Tumour")]] <- ci(roc_data_set_tumour_only)[3]
      AUC_vals[[paste0(model_name, "_Tumour")]] <- roc_data_set_tumour_only$auc
    }
  }
  
  aligned_Final <- merge(Final, aggregateprobs, by = "ID", all.x = TRUE)
  aligned_Final_SLNB <- aligned_Final[aligned_Final$SLNB == "Performed",]
  aligned_Final_SLNB$SLNB_Binary <- ifelse(aligned_Final_SLNB$SLNBpos == "Positive", 1, 0)
  SLNB_data_set <- roc(aligned_Final_SLNB$SLNB_Binary, aligned_Final_SLNB$rec7year)
  AUC_frames[["SLNB"]] <- data.frame(
    Model = "SLNB",
    Sensitivity = SLNB_data_set$sensitivities,
    Specificity = SLNB_data_set$specificities,
    Threshold = SLNB_data_set$thresholds
  )
  
  TestTrainSplit <- list(c("Train", "TrainXGB", "TrainLR"), c("Test", "TestXGB", "TestLR"), 
                         c("Reg", "TestLR", "TestLRIso"),
                         c("NoSLNB", "TestNoSLNBXGB", "TestNoSLNBLR"),
                         c("DelayedRec", "TestDelayedRecXGB", "TestDelayedRecLR"),
                         c("IBIIA", "TestIBIIAXGB", "TestIBIIALR"))
  
  for (TestTrainIndex in 1:length(TestTrainSplit)) {
    TTName <- TestTrainSplit[[TestTrainIndex]][1]
    split_one <- TestTrainSplit[[TestTrainIndex]][2]
    split_two <- TestTrainSplit[[TestTrainIndex]][3]
    split_one_tumour <- paste0(split_one, "_Tumour")
    split_two_tumour <- paste0(split_two, "_Tumour")
    
    if (inc_tumour_in_graph) {
      main_AUC_data <- rbind(AUC_frames[[split_one]], AUC_frames[[split_two]],
                             AUC_frames[[split_one_tumour]], AUC_frames[[split_two_tumour]])
      boundary_points <- data.frame(
        Model = rep(c(split_one, split_two, split_one_tumour, split_two_tumour), each = 2), 
        Sensitivity = c(0, 1, 0, 1, 0, 1, 0, 1), Specificity = c(1, 0, 1, 0, 1, 0, 1, 0), Threshold = c(1, 0, 1, 0, 1, 0, 1, 0)
      )
      optimal_points <- do.call(rbind, optimal_vals[c(split_one, split_two, split_one_tumour, split_two_tumour)])
    } else {
      main_AUC_data <- rbind(AUC_frames[[split_one]], AUC_frames[[split_two]])
      boundary_points <- data.frame(
        Model = rep(c(split_one, split_two), each = 2), Sensitivity = c(0, 1, 0, 1), Specificity = c(1, 0, 1, 0), Threshold = c(1, 0, 1, 0)
      )
      optimal_points <- do.call(rbind, optimal_vals[c(split_one, split_two)])
    }
    
    main_AUC_data <- main_AUC_data %>%
      arrange(Model, desc(Specificity), desc(Sensitivity)) %>% group_by(Model, Specificity) %>% slice(1) %>% ungroup() %>%
      arrange(Model, desc(Sensitivity), desc(Specificity)) %>% group_by(Model, Sensitivity) %>% slice(1) %>% ungroup()
    main_AUC_data <- rbind(main_AUC_data, boundary_points)
    
    XGBColour = "#611683"
    LRColour = "#C06A10"
    LRIsoColour = "#1D4AB0"
    
    XGBLightColour = "#BFA0DB"
    LRLightColour = "#F2C9A3"
    LRIsoLightColour = "#A3C1E8"
    
    if (inc_tumour_in_graph) {
      main_AUC_data$Model <- factor(main_AUC_data$Model, levels = c(split_one, split_two, split_one_tumour, split_two_tumour))
    } else {
      main_AUC_data$Model <- factor(main_AUC_data$Model, levels = c(split_one, split_two))
    }
    
    breaks_graph <- c(split_one, split_two)
    
    if (!inc_tumour_in_graph) {
      if (TTName == "Reg") {
        LR_main_string <- "Ridge logistic regression<br>"
        LR_AUC_string <- paste0("<span style='font-family: Montserrat; font-size:5.4pt;'>AUC: ", sprintf("%.3f", AUC_vals[[split_one]]), 
                                " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_one]]), " - ", sprintf("%.3f", AUC_CI3[[split_one]]), ")</span>")
        LRIso_main_string <- "Binomial logistic regression<br>"
        LRIso_AUC_string <- paste0("<span style='font-family: Montserrat; font-size:5.4pt;'>AUC: ", sprintf("%.3f", AUC_vals[[split_two]]), 
                                   " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_two]]), " - ", sprintf("%.3f", AUC_CI3[[split_two]]), ")</span>")
        main_AUC_labels <- c(paste0(LR_main_string, LR_AUC_string),
                             paste0(LRIso_main_string, LRIso_AUC_string))
        
        cols <- setNames(c(LRColour, LRIsoColour), c(split_one, split_two))
        labs <- setNames(main_AUC_labels, c(split_one, split_two))
      
      } else {
        XGB_main_string <- "Ensemble model (XGBoost)<br>"
        XGB_AUC_string <- paste0("<span style='font-family: Montserrat; font-size:5.4pt;'>AUC: ", sprintf("%.3f", AUC_vals[[split_one]]), 
                                 " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_one]]), " - ", sprintf("%.3f", AUC_CI3[[split_one]]), ")</span>")
        XGB_Sens_string <- paste0("<span style='font-family: Montserrat; font-size:5.4pt;'><br>Sensitivity: ", 
                                 sprintf("%.3f", optimal_vals[[split_one]]$Sensitivity), "; Specificity: ", 
                                 sprintf("%.3f", optimal_vals[[split_one]]$Specificity), "</span>")
        LR_main_string <- "Ridge/L2 logistic regression<br>"
        LR_AUC_string <- paste0("<span style='font-family: Montserrat; font-size:5.4pt;'>AUC: ", sprintf("%.3f", AUC_vals[[split_two]]), 
                                " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_two]]), " - ", sprintf("%.3f", AUC_CI3[[split_two]]), ")</span>")
        LR_Sens_string <- paste0("<span style='font-family: Montserrat; font-size:5.4pt;'><br>Sensitivity: ", 
                                    sprintf("%.3f", optimal_vals[[split_two]]$Sensitivity), "; Specificity: ", 
                                    sprintf("%.3f", optimal_vals[[split_two]]$Specificity), "</span>")
        main_AUC_labels <- c(paste0(XGB_main_string, XGB_AUC_string),
                             paste0(LR_main_string, LR_AUC_string))
        
        cols <- setNames(c(XGBColour, LRColour), c(split_one, split_two))
        labs <- setNames(main_AUC_labels, c(split_one, split_two))
      }
    } else {
      if (TTName == "Reg") {
        next
      } else {
        XGB_main_string <- "Ensemble model (XGBoost)<br><span style='font-family: Montserrat; font-size:4.5pt;'><i>All factors</i><br>"
        XGB_AUC_string <- paste0("AUC: ", sprintf("%.3f", AUC_vals[[split_one]]), 
                                 " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_one]]), " - ", sprintf("%.3f", AUC_CI3[[split_one]]), ")</span>")
        LR_main_string <- "Ridge/L2 logistic regression<br><span style='font-family: Montserrat; font-size:4.5pt;'><i>All factors</i><br>"
        LR_AUC_string <- paste0("AUC: ", sprintf("%.3f", AUC_vals[[split_two]]), 
                                " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_two]]), " - ", sprintf("%.3f", AUC_CI3[[split_two]]), ")</span>")
        XGB_main_string_tumour <- "Ensemble model (XGBoost)<br><span style='font-family: Montserrat; font-size:4.5pt;'><i>Tumour-specific factors only</i><br>"
        XGB_AUC_string_tumour <- paste0("AUC: ", sprintf("%.3f", AUC_vals[[split_one_tumour]]), 
                                        " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_one_tumour]]), " - ", sprintf("%.3f", AUC_CI3[[split_one_tumour]]), ")</span>")
        LR_main_string_tumour <- "Ridge/L2 logistic regression<br><span style='font-family: Montserrat; font-size:4.5pt;'><i>Tumour-specific factors only</i><br>"
        LR_AUC_string_tumour <- paste0("AUC: ", sprintf("%.3f", AUC_vals[[split_two_tumour]]), 
                                       " (95% CI: ", sprintf("%.3f", AUC_CI1[[split_two_tumour]]), " - ", sprintf("%.3f", AUC_CI3[[split_two_tumour]]), ")</span>")
        main_AUC_labels <- c(paste0(XGB_main_string, XGB_AUC_string),
                             paste0(LR_main_string, LR_AUC_string),
                             paste0(XGB_main_string_tumour, XGB_AUC_string_tumour),
                             paste0(LR_main_string_tumour, LR_AUC_string_tumour))
        
        cols <- setNames(c(XGBColour, LRColour, XGBLightColour, LRLightColour), c(split_one, split_two, split_one_tumour, split_two_tumour))
        labs <- setNames(main_AUC_labels, c(split_one, split_two, split_one_tumour, split_two_tumour))
        
        breaks_graph <- c(split_one, split_two, split_one_tumour, split_two_tumour)
      }
    }
    
    main_AUC_plot <- ggplot(main_AUC_data, aes(x = Specificity, y = Sensitivity)) +
      coord_cartesian(ylim = c(0, 1)) +
      #geom_point(data = optimal_points, aes(x = Specificity, y = Sensitivity, color = Model), size = 1, stroke = 0) +
      #geom_step(aes(color = Model), direction = "vh", size = 0.45) +
      geom_smooth(aes(color = Model), method = "gam", formula = y ~ s(x, bs = "cs"), method.args = list(family = quasibinomial(link = "logit")), se = FALSE, size = 0.45) +
      #geom_smooth(aes(color = Model), method = "loess", span = 0.5, se = FALSE, size = 0.45) +
      scale_x_reverse(name = "Specificity", expand = expansion(add = c(0, 0.02)), limits = c(1, 0), n.breaks = 6) + 
      scale_y_continuous(name = "Sensitivity", expand = expansion(add = c(0, 0.02)), n.breaks = 6) +
      geom_segment(x = -1, y = 0, xend = 0, yend = 1, linetype = "dashed") +
      theme_minimal(base_size = 8, base_family = "Montserrat") +
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(1, 0.025),
            legend.title = element_blank(),
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(0, 0, 0, 0),
            legend.text = element_markdown(size = 6, family = "Montserrat SemiBold", lineheight = 1.15),
            legend.key.spacing.y = unit(1.75, "pt"),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.text.x = element_text(color="black", size = 7.2),
            axis.text.y = element_text(color="black", size = 7.2),
            axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), family = "Montserrat SemiBold"),
            axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), family = "Montserrat SemiBold")) +
      scale_color_manual(values = cols, breaks = breaks_graph, labels = labs) +
      guides(color = guide_legend(byrow = TRUE)) +
      labs(color = "Key")
    
    plotfilename = paste0("AUCplot_", TTName, "_", auc_var, ".png")
    ggsave(plotfilename, plot = print(main_AUC_plot), width = 7.5/.cm, height = 7.5/.cm, dpi = 900, grDevices::png)
  }
  
  #Figure 1D
  Final_test <- aligned_Final[aligned_Final$Set == "Test",]
  Final_test$Y_Ground <- factor(Final_test[[auc_var]], levels = c(1, 0), labels = c("Positive", "Negative"))
  
  if (auc_var == "Y_Actual") { # benchmark using recurrence thresholds
    opt_XGB <- optimal_vals$TestXGB$Cutoff
    opt_LR <- optimal_vals$TestLR$Cutoff
    opt_LRIso <- optimal_vals$TestLRIso$Cutoff
    
    if (inc_tumour_in_graph) {
      opt_XGB_T <- optimal_vals$TestXGB_Tumour$Cutoff
      opt_LR_T <- optimal_vals$TestLR_Tumour$Cutoff
    }
  }
  
  colorvar = auc_name
  
  colorstring = paste0(colorvar, "<br><span style='font-family: Montserrat; font-size:5.4pt;'>*7-year follow-up*</span>")
  
  cor_var <- cor.test(Final_test$LR_Prob, Final_test$XGB_Prob, method = "spearman")
  if (cor_var$p.value < 0.001) {
    cor_string <- paste0("Spearman's ρ: ", sprintf("%.3f", cor_var$estimate), "; p < 0.001")
  } else {
    cor_string <- paste0("Spearman's ρ: ", sprintf("%.3f", cor_var$estimate), "; p: ", sprintf("%.3f", cor_var$p.value))
  }
  
  print(cor_string)
  
  cor_var_LR <- cor.test(Final_test$LR_Prob, Final_test$LRIso_Prob, method = "spearman")
  if (cor_var_LR$p.value < 0.001) {
    cor_string_LR <- paste0("Spearman's ρ: ", sprintf("%.3f", cor_var_LR$estimate), "; p < 0.001")
  } else {
    cor_string_LR <- paste0("Spearman's ρ: ", sprintf("%.3f", cor_var_LR$estimate), "; p: ", sprintf("%.3f", cor_var_LR$p.value))
  }
  
  FACS_plot <- ggplot(Final_test, aes(x = LR_Prob, y = XGB_Prob, color = Y_Ground)) +
    geom_point(size = 0.75, stroke = 0) +
    scale_x_continuous(name = "Predicted probability (LR)", 
                       expand = c(0, 0), limits = c(0, 1), breaks = c(0, opt_LR, 1),
                       labels = c("0", sprintf("%.2f", opt_LR), "1"),
                       trans = trans_new(name = "sqrt", transform = sqrt, inverse = function(x) x^2)) + 
    scale_y_continuous(name = "Predicted probability (XGBoost)", 
                       expand = c(0, 0), limits = c(0, 1), breaks = c(0, opt_XGB, 1),
                       labels = c("0", sprintf("%.2f", opt_XGB), "1"),
                       trans = trans_new(name = "sqrt", transform = sqrt, inverse = function(x) x^2)) +
    #geom_text(x = 0.01, y = 0.01, label = cor_string, color = "black", size = 6/.pt, hjust = 0, family = "Montserrat", inherit.aes = FALSE) +
    geom_segment(x = 0, y = opt_XGB^0.5, xend = 1, yend = opt_XGB^0.5, linetype = "solid", color = XGBColour, linewidth = 0.3) +
    geom_segment(x = opt_LR^0.5, y = 0, xend = opt_LR^0.5, yend = 1, linetype = "solid", color = LRColour, linewidth = 0.3) +
    theme_minimal(base_size = 7.2, base_family = "Montserrat") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = c(auc_pos, 0.1),
          legend.title = element_markdown(size = 6, family = "Montserrat SemiBold", lineheight = 1.2,
                                          margin = margin(t = 0, r = 0, b = -2, l = 0)),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.text = element_markdown(size = 5.4, family = "Montserrat"),
          legend.spacing.y = unit(-10, "pt"),
          legend.key.spacing.y = unit(-10, "pt"),
          legend.key.width = unit(4, "pt"),
          axis.ticks = element_line(color = "black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          axis.text.x.bottom = element_text(color="black", size = 7.2, margin = margin(t = 2, r = 0, b = 0, l = 0)),
          axis.text.y.left = element_text(color="black", size = 7.2, margin = margin(t = 0, r = 2, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), family = "Montserrat SemiBold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), family = "Montserrat SemiBold")) +
    scale_color_manual(values = c("Negative" = "#00005A", "Positive" = "#CC0000"),
                       labels = c("Negative" = "Negative", "Positive" = "Positive")) +
    guides(color = guide_legend(byrow = TRUE, override.aes = list(size = 0.5))) +
    labs(color = colorstring)
  
  plotfilename = paste0("FACSplot_", auc_var, ".png")
  ggsave(plotfilename, plot = print(FACS_plot), width = 7.5/.cm, height = 7.5/.cm, dpi = 900, grDevices::png)
  
  FACS_plot_LR <- ggplot(Final_test, aes(x = LR_Prob, y = LRIso_Prob, color = Y_Ground)) +
    geom_point(size = 0.75, stroke = 0) +
    scale_x_continuous(name = "Predicted probability (L2-regularized LR)", 
                       expand = c(0, 0), limits = c(0, 1), breaks = c(0, opt_LR, 1),
                       labels = c("0", sprintf("%.2f", opt_LR), "1"),
                       trans = trans_new(name = "sqrt", transform = sqrt, inverse = function(x) x^2)) + 
    scale_y_continuous(name = "Predicted probability (Non-regularized LR)", 
                       expand = c(0, 0), limits = c(0, 1), breaks = c(0, opt_LRIso, 1),
                       labels = c("0", sprintf("%.2f", opt_LRIso), "1"),
                       trans = trans_new(name = "sqrt", transform = sqrt, inverse = function(x) x^2)) +
    #geom_text(x = 0.01, y = 0.01, label = cor_string_LR, color = "black", size = 6/.pt, hjust = 0, family = "Montserrat", inherit.aes = FALSE) +
    geom_segment(x = 0, y = opt_LRIso^0.5, xend = 1, yend = opt_LRIso^0.5, linetype = "solid", color = LRIsoColour, linewidth = 0.3) +
    geom_segment(x = opt_LR^0.5, y = 0, xend = opt_LR^0.5, yend = 1, linetype = "solid", color = LRColour, linewidth = 0.3) +
    theme_minimal(base_size = 7.2, base_family = "Montserrat") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = c(auc_pos, 0.1),
          legend.title = element_markdown(size = 6, family = "Montserrat SemiBold", lineheight = 1.2,
                                          margin = margin(t = 0, r = 0, b = -2, l = 0)),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.text = element_markdown(size = 5.4, family = "Montserrat"),
          legend.spacing.y = unit(-10, "pt"),
          legend.key.spacing.y = unit(-10, "pt"),
          legend.key.width = unit(4, "pt"),
          axis.ticks = element_line(color = "black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          axis.text.x.bottom = element_text(color="black", size = 7.2, margin = margin(t = 2, r = 0, b = 0, l = 0)),
          axis.text.y.left = element_text(color="black", size = 7.2, margin = margin(t = 0, r = 2, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), family = "Montserrat SemiBold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), family = "Montserrat SemiBold")) +
    scale_color_manual(values = c("Negative" = "#00005A", "Positive" = "#CC0000"),
                       labels = c("Negative" = "Negative", "Positive" = "Positive")) +
    guides(color = guide_legend(byrow = TRUE, override.aes = list(size = 0.5))) +
    labs(color = colorstring)
  
  plotfilename = paste0("FACSplotLR_", auc_var, ".png")
  ggsave(plotfilename, plot = print(FACS_plot_LR), width = 7.5/.cm, height = 7.5/.cm, dpi = 900, grDevices::png)
}

#Table 1
selected_data <- read.csv("selected_data.csv")
feature_names <- colnames(selected_data)[!colnames(selected_data) %in% c("rec7year", "ID")]

remove_columns <- c("Site", "HistoComType")
add_columns <- c("SiteThreeLevel", "HistoSubtype")

feature_names <- feature_names[!feature_names %in% remove_columns]
feature_names <- c(feature_names, add_columns)

FinalLR <- aligned_Final[, c(feature_names, "ID", "rec2year", "rec7year", "dm7year", "meldeath7year", "TStageFull", "ClinicalStage", "Set")]
FinalLR$NumPregnancies[FinalLR$Sex == "Male"] <- 0

# FinalLR$BDCat <- cut(FinalLR$BD, breaks = c(0, 1, 2, 4, Inf), labels = c("T1", "T2", "T3", "T4"), right = FALSE)
# FinalLR$BDCat <- factor(FinalLR$BDCat, levels = c("T1", "T2", "T3", "T4"))
# feature_names <- c(feature_names, "BDCat")

categorical_columns <- sapply(FinalLR[, feature_names], is.factor) | 
  sapply(FinalLR[, feature_names], is.character)

continuous_columns <- sapply(FinalLR[, feature_names], is.numeric)

subset_LR <- FinalLR[, c(feature_names, "rec7year", "dm7year", "meldeath7year")]
subset_LR$SiteThreeLevel <- relevel(factor(subset_LR$SiteThreeLevel), ref = "Trunk")
subset_LR$HistoSubtype <- relevel(factor(subset_LR$HistoSubtype), ref = "SSM")
#subset_LR$SLNBTotal <- relevel(factor(subset_LR$SLNBTotal), ref = "Not performed")

outcome_list <- c("rec7year", "dm7year", "meldeath7year")
for (outcome in outcome_list) {
  print(paste0("OUTCOME: ", outcome))
  model_list <- list()
  
  for (feature_index in 1:length(feature_names)) {
    feature <- feature_names[feature_index]
    formula <- as.formula(paste0(outcome, " ~ ", feature))
    model <- glm(formula, data = subset_LR, family = binomial)
    
    model_list[[feature]] <- model
  }
  
  print_models <- c("Age", "Sex", "BD", "Ulceration", "MitoticRate", "SiteThreeLevel", "HistoSubtype") #, "SLNBTotal"
  for (model_name in print_models) {
    cur_model <- model_list[[model_name]]
    
    print(model_name)
    print(exp(coef(summary(cur_model))[,"Estimate"]))
    print(exp(confint(profile(cur_model))))
    print(coef(summary(cur_model))[, "Pr(>|z|)"])
  }
}

multivariate_model <- glm(as.formula(paste("rec7year ~", paste(print_models, collapse = " + "))), data = subset_LR, family = binomial)
print(exp(coef(summary(multivariate_model))[,"Estimate"]))
print(exp(confint(multivariate_model)))
print(coef(summary(multivariate_model))[, "Pr(>|z|)"])

multivariate_model_dm <- glm(as.formula(paste("dm7year ~", paste(print_models, collapse = " + "))), data = subset_LR, family = binomial)
print(exp(coef(summary(multivariate_model_dm))[,"Estimate"]))
print(exp(confint(multivariate_model_dm)))
print(coef(summary(multivariate_model_dm))[, "Pr(>|z|)"])

multivariate_model_death <- glm(as.formula(paste("meldeath7year ~", paste(print_models, collapse = " + "))), data = subset_LR, family = binomial)
print(exp(coef(summary(multivariate_model_death))[,"Estimate"]))
print(exp(confint(multivariate_model_death)))
print(coef(summary(multivariate_model_death))[, "Pr(>|z|)"])

TableOneColumns <- c("Age", "Sex", "BD", "Ulceration", "MitoticRate", "SiteThreeLevel", "HistoSubtype", "SLNBTotal")
table1_stratified <- CreateTableOne(vars = TableOneColumns, strata = "rec2year", data = FinalLR, addOverall = TRUE); print(table1_stratified)
table1_stratified <- CreateTableOne(vars = TableOneColumns, strata = "rec7year", data = FinalLR, addOverall = TRUE); print(table1_stratified)
table1_stratified <- CreateTableOne(vars = TableOneColumns, strata = "dm7year", data = FinalLR, addOverall = TRUE); print(table1_stratified)
table1_stratified <- CreateTableOne(vars = TableOneColumns, strata = "meldeath7year", data = FinalLR, addOverall = TRUE); print(table1_stratified)
table1_train_test <- CreateTableOne(vars = TableOneColumns, strata = "Set", data = FinalLR, addOverall = TRUE, test = TRUE); print(table1_train_test)

TableTColumns <- c("TStageFull", "ClinicalStage", "SLNBTotal")
tableT_stratified <- CreateTableOne(vars = TableTColumns, strata = "rec2year", data = FinalLR, addOverall = TRUE); print(tableT_stratified)
tableT_stratified <- CreateTableOne(vars = TableTColumns, strata = "rec7year", data = FinalLR, addOverall = TRUE); print(tableT_stratified)
tableT_stratified <- CreateTableOne(vars = TableTColumns, strata = "dm7year", data = FinalLR, addOverall = TRUE); print(tableT_stratified)
tableT_stratified <- CreateTableOne(vars = TableTColumns, strata = "meldeath7year", data = FinalLR, addOverall = TRUE); print(tableT_stratified)

#Figure 2
feature_dict <- c(
  # Tumour
  "BD" = "Breslow depth",
  "Ulceration_Present" = "Ulceration",
  "MitoticRate" = "Mitotic rate",
  "Regression_Present" = "Regression",
  
  "HistoComType_Nodular" = "Histological subtype: Nodular",
  "HistoComType_Other" = "Histological subtype: Other",
  "PreNaevus_Present" = "Histologic association with naevus",
  
  "Site_HeadNeck" = "Melanoma site: Head/Neck",
  "Site_Arms" = "Melanoma site: Extremities",
  
  "SLNBTotal_Positive" = "SLNB status: Positive",
  "SLNBTotal_Negative" = "SLNB status: Negative",
  
  # Sociodemographic
  "Age" = "Age",
  "Sex_Male" = "Sex: Male",
  
  "RemotenessCat_Inner reg." = "Remoteness: Inner regional",
  "RemotenessCat_Outer reg./Remote" = "Remoteness: Outer regional/remote",
  
  "Birthplace_Elsewhere" = "Birthplace: Not in Aus/NZ",
  
  "SEIFADecile" = "Socioeconomic status (SEIFA)",
  
  "Occupation_Trade/elem. clerical" = "Occupation: Trades and clerical",
  "Occupation_Inter./adv. clerical" = "Occupation: Admin/service work",
  "Occupation_Professional" = "Occupation: Professional",
  
  "Education_Grade 12/Trade/Diploma" = "Education: High school/Trades",
  "Education_Uni/College" = "Education: University/College",
  
  "Partner_Yes" = "Established partner at diagnosis",
  "NumPregnancies" = "Number of pregnancies",
  
  # Health behaviour
  "EverSmoker_Yes" = "Ever smoker: Yes",
  
  "Sunscreen_Yes" = "Daily sunscreen use",
  "SunOccup_Mixed" = "Occupational sun exposure: Mixed",
  "SunOccup_Outdoors" = "Occupational sun exposure: Outdoors",
  
  "SelfSkinChecks_Once/yr or more" = "Self skin checks: Annually or more",
  "DoctorSkinChecks_Once/yr or more" = "Doctor skin checks: Annually or more",
  
  "FACTG_Emotional" = "FACT-G Emotional sub-score",
  "FACTG_Functional" = "FACT-G Functional sub-score",
  "FACTG_Physical" = "FACT-G Physical sub-score",
  "FACTG_Social" = "FACT-G Social sub-score",
  
  # Medical
  "PrevMelanoma_Yes" = "Previous diagnosis of melanoma",
  "PrevSkinCancer_Yes" = "Previous diagnosis of skin cancer",
  "FHxMelanoma_Yes" = "Family history of melanoma",
  
  "HeartDisease_Yes" = "Cardiovascular co-morbidity",
  "RespDisease_Yes" = "Respiratory co-morbidity",
  "KidneyDisease_Yes" = "Renal co-morbidity",
  
  "Arthritis_Yes" = "Diagnosis of arthritis",
  "Dementia_Yes" = "Diagnosis of dementia",
  "Depression_Yes" = "Diagnosis of depression",
  "Diabetes_Yes" = "Diagnosis of diabetes",
  "HIV_Yes" = "Diagnosis of HIV/AIDS",
  "Hypertension_Yes" = "Diagnosis of hypertension",
  "OtherMalignancy_Yes" = "Concomitant non-skin malignancy",
  
  "BMI" = "Body mass index",
  
  "CorticosteroidLongTerm_Yes" = "Long-term corticosteroid use",
  "NSAIDLongTerm_Yes" = "Long-term NSAID use",
  "PsychotropicsLongTerm_Yes" = "Long-term psychotropic use",
  "StatinsLongTerm_Yes" = "Long-term statin use",

  "Chemo5YR_Yes" = "Chemotherapy, 5y pre-diagnosis",
  "Rad5YR_Yes" = "Radiotherapy, 5y pre-diagnosis",
  "Hormone5YR_Yes" = "Hormone therapy, 5y pre-diagnosis",
  
  "HealthServiceAccess_Yes" = "Liaison with cancer services"
)

SHAP_annotations <- read.csv("SHAP.csv")

SHAP_annotations$ScaledSHAPXGB <- SHAP_annotations$MeanAbsSHAP_XGB / sum(SHAP_annotations$MeanAbsSHAP_XGB)
SHAP_annotations$ScaledSHAPLR <- SHAP_annotations$MeanAbsSHAP_LR / sum(SHAP_annotations$MeanAbsSHAP_LR)

feature_order <- SHAP_annotations %>%
  arrange(desc(ScaledSHAPXGB), desc(ScaledSHAPLR)) %>%
  pull(Feature)

SHAP_annotations_long <- SHAP_annotations %>%
  mutate(ScaledSHAPXGB = -ScaledSHAPXGB) %>%
  dplyr::select(Feature, ScaledSHAPXGB, ScaledSHAPLR, LR_Direction, XGB_Direction) %>%
  pivot_longer(cols = c("ScaledSHAPXGB", "ScaledSHAPLR"),
               names_to = "Model",
               values_to = "ScaledSHAP") %>%
  mutate(Model = factor(Model,
                        levels = c("ScaledSHAPXGB", "ScaledSHAPLR"),
                        labels = c("XGB", "LR"))) %>%
  mutate(Direction = if_else(Model == "XGB", XGB_Direction, LR_Direction)) %>%
  group_by(Feature) %>%
  mutate(
    xgb_val = if_else(Model == "XGB", ScaledSHAP, NA_real_),
    lr_val  = if_else(Model == "LR", ScaledSHAP, NA_real_)) %>%
  tidyr::fill(xgb_val, lr_val, .direction = "downup") %>%
  ungroup() %>%
  #arrange(desc(xgb_val), lr_val) %>%
  arrange((abs(xgb_val) + abs(lr_val)), lr_val) %>%
  mutate(Feature = recode(Feature, !!!feature_dict)) %>%
  mutate(Feature = factor(Feature, levels = unique(Feature))) %>%
  dplyr::select(Feature, Model, ScaledSHAP, Direction)

unique_features <- unique(SHAP_annotations_long$Feature)
ymax <- length(unique_features) + 0.35

if (TRUE) {
  SHAP_annotations_long <- SHAP_annotations_long[35:108,] # cutoff for SHAP sum > 0.01 i.e. shared 1% importance across both models 
  # or average of 0.5%, but most features below 10% are specific to only one model in any case
}

shapley_string <- paste0("<span style='font-size:7.2pt;'>Recurrence risk</span>", 
                         "<br><span style='font-family: Montserrat; font-size:5.4pt;'>*7-year follow-up*</span>")

SHAPLEY_plot <- ggplot(SHAP_annotations_long, aes(x = ScaledSHAP, y = Feature, fill = Direction)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  scale_x_continuous(name = "Scaled SHAP values", labels = c("0.15", "0.10", "0.05", "0.00", "0.05", "0.10"), 
                     limits = c(-0.15, 0.1), breaks = seq(-0.15, 0.1, 0.05),
                     expand = expansion(mult = 0.02)) +
  #facet_grid(~ Model, scales = "free_x", space = "free") +
  scale_fill_manual(name = shapley_string, values = c("Positive" = "#880000", "Negative" = "#000090"),
                    labels = c("Positive" = "Increased", "Negative" = "Decreased")) +
  theme_minimal(base_size = 7.2, base_family = "Montserrat") +
  theme(axis.text.y = element_text(size = 7.2), 
        axis.text.x = element_text(size = 7.2), 
        legend.title = element_markdown(size = 7.2, family = "Montserrat SemiBold"),
        legend.text = element_text(size = 7.2, margin = margin(t = 0.9, r = 0, b = 0.9, l = 3)),
        legend.key.spacing.y = unit(5, "pt"), 
        axis.ticks = element_line(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        legend.position = "right", 
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0), family = "Montserrat"),
        axis.title.y = element_blank()) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.3) +
  guides(fill = guide_legend(keyheight = unit(2, "pt"), keywidth = unit(9, "pt"), byrow = TRUE))

ggsave("SHAPLEYplot.png", plot = print(SHAPLEY_plot), width = 16.5/.cm, height = 14.5/.cm, dpi = 900, grDevices::png)

# Binomial all factors
subset_LR_b <- FinalLR
print_models <- c("Age", "Sex", "BD", "Ulceration", "MitoticRate", "SiteThreeLevel", "HistoSubtype") #, "SLNBTotal"

factor_map <- list(
  Sex = c("Female", "Male"),
  Ulceration = c("Absent", "Present"),
  Regression = c("Absent", "Present"),
  PreNaevus = c("Absent", "Present"),
  PrevMelanoma = c("No", "Yes"),
  PrevSkinCancer = c("No", "Yes"),
  FHxMelanoma = c("No", "Yes"),
  Diabetes = c("No", "Yes"),
  Hypertension = c("No", "Yes"),
  HeartDisease = c("No", "Yes"),
  RespDisease = c("No", "Yes"),
  KidneyDisease = c("No", "Yes"),
  Arthritis = c("No", "Yes"),
  Depression = c("No", "Yes"),
  HIV = c("No", "Yes"),
  Dementia = c("No", "Yes"),
  OtherMalignancy = c("No", "Yes"),
  Chemo5YR = c("No", "Yes"),
  Rad5YR = c("No", "Yes"),
  Hormone5YR = c("No", "Yes"),
  StatinsLongTerm = c("No", "Yes"),
  NSAIDLongTerm = c("No", "Yes"),
  CorticosteroidLongTerm = c("No", "Yes"),
  PsychotropicsLongTerm = c("No", "Yes"),
  EverSmoker = c("No", "Yes"),
  Sunscreen = c("No", "Yes"),
  SelfSkinChecks = c("Less than once/yr", "Once/yr or more"),
  DoctorSkinChecks = c("Less than once/yr", "Once/yr or more"),
  HealthServiceAccess = c("No", "Yes"),
  Partner = c("No", "Yes"),
  Birthplace = c("Aus/NZ", "Elsewhere"),
  RemotenessCat = c("Major cities", "Inner reg.", "Outer reg./Remote"),
  Education = c("Less than Grade 12", "Grade 12/Trade/Diploma", "Uni/College"),
  Occupation = c("Trade/elem. clerical", "Inter./adv. clerical", "Professional", "Unemployed"),
  SiteThreeLevel = c("Trunk", "HeadNeck", "Limbs"),
  HistoSubtype = c("SSM", "Nodular", "Other", "Unclassified"),
  SunOccup = c("Indoors", "Mixed", "Outdoors"),
  SLNBTotal = c("Not performed", "Negative", "Positive")
)

for (name in intersect(names(factor_map), names(subset_LR_b))) {
  subset_LR_b[[name]] <- factor(subset_LR_b[[name]], levels = factor_map[[name]])
}

unique_features <- setdiff(feature_names, print_models)
var_order <- c("Regression", "PreNaevus", 
               "SEIFADecile", "NumPregnancies", "RemotenessCat", "Occupation", "Education", "Partner", "Birthplace", 
               "PrevMelanoma", "PrevSkinCancer", "FHxMelanoma", 
               "HeartDisease", "RespDisease", "KidneyDisease", 
               "Arthritis", "Dementia", "Depression", "Diabetes", "HIV", "Hypertension", 
               "OtherMalignancy", "Chemo5YR", "Rad5YR", "Hormone5YR", 
               "CorticosteroidLongTerm", "NSAIDLongTerm", "PsychotropicsLongTerm", "StatinsLongTerm", 
               "FACTG_Emotional", "FACTG_Functional", "FACTG_Physical", "FACTG_Social", 
               "EverSmoker", "Sunscreen", "SunOccup", "SelfSkinChecks", "DoctorSkinChecks", "HealthServiceAccess")

unique_features <- unique_features[unique_features %in% var_order]
unique_features <- var_order[var_order %in% unique_features]

table1_extra <- CreateTableOne(vars = unique_features, data = subset_LR_b, addOverall = TRUE); print(table1_extra, showAllLevels = TRUE)
table1_extra_rec <- CreateTableOne(vars = unique_features, strata = "rec7year", data = subset_LR_b, addOverall = TRUE); print(table1_extra_rec, showAllLevels = TRUE)
table1_extra_dm <- CreateTableOne(vars = unique_features, strata = "dm7year", data = subset_LR_b, addOverall = TRUE); print(table1_extra_dm, showAllLevels = TRUE)
table1_extra_md <- CreateTableOne(vars = unique_features, strata = "meldeath7year", data = subset_LR_b, addOverall = TRUE); print(table1_extra_md, showAllLevels = TRUE)

outcome_vars <- c("rec7year", "dm7year", "meldeath7year")

for (outcome_var in outcome_vars) {
  cat("\n", outcome_var, "\n", sep = "")
  
  for (uf_index in 1:length(unique_features)) {
    uf <- unique_features[uf_index]
    
    uv_m <- glm(as.formula(paste(outcome_var, " ~ ", uf)), data = subset_LR_b, family = binomial)
    uv_s <- summary(uv_m)$coef
    uv_ci <- confint(profile(uv_m))
    
    uv_idx <- rownames(uv_s) != "(Intercept)"
    uv_strings <- sprintf("%s %s (%s-%s) p = %s",
      rownames(uv_s)[uv_idx],
      sprintf("%.2f", exp(uv_s[uv_idx, "Estimate"])),
      sprintf("%.2f", exp(uv_ci[uv_idx, 1])),
      sprintf("%.2f", exp(uv_ci[uv_idx, 2])),
      sprintf("%.3f", uv_s[uv_idx, "Pr(>|z|)"])
    )
    
    cat(paste(uv_strings, collapse = "\n"), "\n", sep = "")
  }
}

for (outcome_var in outcome_vars) {
  cat("\n", outcome_var, "\n", sep = "")
  
  for (uf_index in 1:length(unique_features)) {
    uf <- unique_features[uf_index]
    
    mv_m <- glm(as.formula(paste(outcome_var, "~", uf, " + ", paste(print_models, collapse = " + "))), 
                data = subset_LR_b, family = binomial)
    mv_s <- summary(mv_m)$coef
    mv_ci <- confint(profile(mv_m))
    
    mv_idx <- grepl(uf, rownames(mv_s), fixed = TRUE)
    mv_strings <- sprintf("%s %s (%s-%s) p = %s",
                          rownames(mv_s)[mv_idx],
                          sprintf("%.2f", exp(mv_s[mv_idx, "Estimate"])),
                          sprintf("%.2f", exp(mv_ci[mv_idx, 1])),
                          sprintf("%.2f", exp(mv_ci[mv_idx, 2])),
                          sprintf("%.3f", mv_s[mv_idx, "Pr(>|z|)"])
    )
    
    cat(paste(mv_strings, collapse = "\n"), "\n", sep = "")
  }
}

#Figure 3
train_data <- aligned_Final[aligned_Final$Set == "Train",]
test_data <- aligned_Final[aligned_Final$Set == "Test",]

StageIBIIA_data <- test_data[test_data$Stage %in% c("IB", "IIA"),]
SLNBnone_data <- test_data[test_data$SLNBTotal == "Not performed",]
SLNBneg_data <- test_data[test_data$SLNBTotal == "Negative",]

DelayedRec_data <- test_data[test_data$rec2year == 0,]
DistantRec_data <- test_data
MelanomaDeath_data <- test_data

XGBColours = list(lighten(XGBColour, 0.45), darken(XGBColour, 0.6))
LRColours = list(lighten(LRColour, 0.3), darken(LRColour, 0.375))

data_sources <- list(test_data, test_data, SLNBnone_data, SLNBnone_data, SLNBneg_data, SLNBneg_data, StageIBIIA_data, StageIBIIA_data, 
                     DelayedRec_data, DelayedRec_data, DistantRec_data, DistantRec_data, MelanomaDeath_data, MelanomaDeath_data, 
                     test_data, test_data, test_data, test_data, train_data, train_data, train_data, train_data)
matrix_names <- c(rep(c("XGB", "LR"), 7),
                  "XGB", "LR", "XGB_T", "LR_T", "XGB", "LR", "XGB_T", "LR_T")
matrix_titles_fore <- rep(c("XGBoost prediction", "LR prediction"), 11)
matrix_titles_back <- c("Recurrence", "Recurrence", "Subgroup: no SLNB", "Subgroup: no SLNB",  "Subgroup: SLNB-negative", "Subgroup: SLNB-negative", 
                        "Subgroup: Stage IB/IIA", "Subgroup: Stage IB/IIA", "Subgroup: delayed recurrence", "Subgroup: delayed recurrence",
                        "Distant metastasis", "Distant metastasis", 
                        "Melanoma-specific death", "Melanoma-specific death",
                        "Recurrence (all factors)", "Recurrence (all factors)", "Recurrence (tumour-specific)", "Recurrence (tumour-specific)", 
                        "Recurrence (all factors)", "Recurrence (all factors)", "Recurrence (tumour-specific)", "Recurrence (tumour-specific)")
matrix_filenames <- c("XGB Test", "LR Test", "XGB SLNBnone", "LR SLNBnone", "XGB SLNBneg", "LR SLNBneg", "XGB IBIIA", "LR IBIIA",
                      "XGB DelayedRec", "LR DelayedRec", "XGB DistantRec", "LR DistantRec", "XGB MelanomaDeath", "LR MelanomaDeath",
                      "XGB Test-S", "LR Test-S", "XGB Test-S TO", "LR Test-S TO", "XGB Train-S", "LR Train-S", "XGB Train-S TO", "LR Train-S TO")
KM_colours <- rep(list(XGBColours, LRColours), 11)
cutoffs <- c(rep(c(opt_XGB, opt_LR), 8), opt_XGB_T, opt_LR_T, opt_XGB, opt_LR, opt_XGB_T, opt_LR_T)

# note run 8 is a replicate of run 1 with a more detailed axis title for supplementary figure

for (matrix_index in 1:length(matrix_names)) {
  model_name <- matrix_names[matrix_index]
  matrix_title_fore <- matrix_titles_fore[matrix_index]
  matrix_title_back <- matrix_titles_back[matrix_index]
  cutoff_value <- cutoffs[matrix_index]
  data_source <- data_sources[[matrix_index]]
  matrix_filename <- matrix_filenames[matrix_index]
  
  matrix_title <- paste0(matrix_title_fore, "<br><span style='font-family: Montserrat; font-size:4.8pt;'><i>", matrix_title_back, "</i></span>")
  
  predicted_val <- ifelse(data_source[[paste0(model_name, "_Prob")]] >= cutoff_value, 1, 0)
  
  if (matrix_title_back == "Melanoma-specific death") {
    actual_val <- data_source$Y_MelDeath
    event_label <- c("Non-fatal", "Fatal")
  } else if (matrix_title_back == "Distant metastasis") {
    actual_val <- data_source$Y_DM
    event_label <- c("No DM", "DM")
  } else {
    actual_val <- data_source$Y_Actual
    event_label <- c("No rec.", "Rec.")
  }
  
  
  predicted_val <- factor(predicted_val, levels = c(0, 1), labels = event_label)
  actual_val <- factor(actual_val, levels = c(0, 1), labels = event_label)
  
  cm_table <- table(actual = actual_val, predicted = predicted_val)
  melted_matrix <- reshape2::melt(cm_table, na.rm = TRUE)
  melted_matrix <- melted_matrix %>%
    group_by(actual) %>%
    mutate(total = sum(value)) %>%
    ungroup() %>%
    mutate(percentage = ifelse(total == 0, 0, value / total)) %>%
    mutate(doubletext = ifelse(value == 0 & predicted != actual, TRUE, FALSE)) %>%
    mutate(valtext = ifelse(doubletext, "", value)) %>%
    mutate(val2text = ifelse(doubletext, 0, "")) %>%
    mutate(perctext = ifelse(doubletext, "", paste0(sprintf("%.1f", percentage*100), "%"))) %>%
    mutate(colourperc = 0.925 - (percentage^0.8)*0.5) %>%
    mutate(colour = ifelse(predicted == actual, rgb(colourperc, 1, colourperc), 
                           ifelse(percentage == 0, rgb(0.92, 0.92, 0.92), rgb(1, colourperc, colourperc)))) %>%
    mutate(border = ifelse(predicted == actual, "#FFFFFF", "#FFFFFF")) %>%
    dplyr::select(-total)
  
  hlabel <- matrix_title
  
  initial_labels_h <- c("No rec.", "Rec.")
  initial_labels_v <- c("No rec.", "Rec.")
  
  pred0 = sum(predicted_val == event_label[1])
  pred1 = sum(predicted_val == event_label[2])
  actual0 = sum(actual_val == event_label[1])
  actual1 = sum(actual_val == event_label[2])

  if (matrix_title_back == "Melanoma-specific death") {
    unique_labels_h <- c(paste0("Low (", pred0, ")"), 
                         paste0("High (", pred1, ")"))
    unique_labels_v <- c(paste0("Non-fatal (", actual0, ")"), 
                         paste0("Fatal (", actual1, ")"))
    conf_label <- "7-year melanoma mortality"
    KM_label <- "Melanoma-specific survival"
  } else if (matrix_title_back == "Distant metastasis") {
    unique_labels_h <- c(paste0("Low (", pred0, ")"), 
                         paste0("High (", pred1, ")"))
    unique_labels_v <- c(paste0("No DM (", actual0, ")"), 
                         paste0("DM (", actual1, ")"))
    conf_label <- "7-year distant met."
    KM_label <- "Distant metastasis-free survival"
  } else {
    unique_labels_h <- c(paste0("Low (", pred0, ")"), 
                         paste0("High (", pred1, ")"))
    unique_labels_v <- c(paste0("No rec. (", actual0, ")"), 
                         paste0("Rec. (", actual1, ")"))
    conf_label <- "7-year recurrence"
    KM_label <- "Recurrence-free survival"
  }
  
  conf_matrix <- ggplot(melted_matrix, aes(x = predicted, y = actual, fill = colour)) +
    geom_tile(width = 1, color = "#FFFFFF", linewidth = 0.75/.pt) + scale_fill_identity() +
    scale_x_discrete(name = hlabel, position = "top", expand = c(0, 0), labels = rev(unique_labels_h), limits = rev(event_label)) +
    scale_y_discrete(name = conf_label, expand = c(0, 0), labels = unique_labels_v) +
    theme_bw(base_size = 6, base_family = "Montserrat") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.75 * .pt/.stroke),
          legend.position = "none",
          plot.margin = margin(1, 1, 1, 1),
          axis.ticks = element_blank(),
          axis.text.x.top = element_markdown(color = "black", size = 6,
                                     margin = margin(t = 0, r = 0, b = 2, l = 0)),
          axis.text.y = element_markdown(color = "black", size = 6, hjust = 0.5, angle = 90, 
                                     margin = margin(t = 0, r = 2.5, b = 0, l = 0)),
          axis.title.x.top = element_markdown(margin = margin(t = 0, r = 0, b = 4, l = 0), size = 6, family = "Montserrat SemiBold", lineheight = 1),
          axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0), size = 6, family = "Montserrat SemiBold")) +
    coord_fixed() +
    geom_text(aes(predicted, actual, label = valtext), color = "black", family = "Montserrat", size = 8/.pt, vjust = -0.15) +
    geom_text(aes(predicted, actual, label = perctext), color = "black", family = "Montserrat SemiBold", size = 6/.pt, vjust = 1.45) +
    geom_text(aes(predicted, actual, label = val2text), color = "black", family = "Montserrat", size = 8/.pt, vjust = 0.5)
  filename = paste("ConfMatrix - ", matrix_filename, ".png", sep = "")
  ggsave(filename, plot = print(conf_matrix), width = 4/.cm, height = 4/.cm, dpi = 900, grDevices::png)
  
  # hook CoxPH into same loop
  KM_colour_palette <- KM_colours[[matrix_index]]
  data_source$predicted_val <- predicted_val
  
  if (matrix_title_back == "Melanoma-specific death") {
    surv_object <- Surv(time = data_source$survFU, event = data_source$meldeath7year)
    scale_values = c("predicted_val=Non-fatal" = KM_colour_palette[[1]], "predicted_val=Fatal" = KM_colour_palette[[2]])
    scale_labels = c("predicted_val=Non-fatal" = "Low", "predicted_val=Fatal" = "High")
  } else if (matrix_title_back == "Distant metastasis") {
    surv_object <- Surv(time = data_source$dmFU, event = data_source$dm7year)
    scale_values = c("predicted_val=No DM" = KM_colour_palette[[1]], "predicted_val=DM" = KM_colour_palette[[2]])
    scale_labels = c("predicted_val=No DM" = "Low", "predicted_val=DM" = "High")
  } else {
    surv_object <- Surv(time = data_source$recFU, event = data_source$rec7year)
    scale_values = c("predicted_val=No rec." = KM_colour_palette[[1]], "predicted_val=Rec." = KM_colour_palette[[2]])
    scale_labels = c("predicted_val=No rec." = "Low", "predicted_val=Rec." = "High")
  }
  
  rec_fit <- survfit(surv_object ~ predicted_val, data = data_source)

  cox_fit <- coxph(surv_object ~ predicted_val, data = data_source)
  
  pstring <- ifelse(summary(cox_fit)$coefficients[5] < 0.001, 
                    paste0("p < 0.001"), 
                    paste0("p: ", sprintf("%.3f", summary(cox_fit)$coefficients[5])))
  cox_text <- paste0("HR: ", sprintf("%.2f", exp(coef(cox_fit))), 
                     "; ", pstring,
                     "; 95% CI: ", sprintf("%.2f", exp(confint(cox_fit)[1])), "-", sprintf("%.2f", exp(confint(cox_fit)[2])))
  
  print(cox_text)
  
  matrix_small_title <- gsub("prediction", "model", matrix_title)
  matrix_small_title <- gsub("delayed recurrence", "delayed rec.", matrix_small_title)
  #matrix_small_title <- gsub("metastasis", "met.", matrix_small_title)
  matrix_small_title <- gsub("Melanoma-specific death", "Melanoma mortality", matrix_small_title)
  
  survplot <- ggsurvfit(rec_fit, size = 0.3) +
    labs(x = "Months", y = KM_label, color = matrix_small_title) +
    geom_text(x = 2.6, y = 0.042, label = cox_text, size = 5/.pt, hjust = 0, vjust = 0, family = "Montserrat", inherit.aes = FALSE, check_overlap = TRUE) +
    theme_minimal(base_size = 6, base_family = "Montserrat") +
    theme(
      legend.title = element_markdown(size = 6, family = "Montserrat SemiBold", hjust = 0,
                                      margin = margin(t = 0, r = 3, b = 0, l = -8.5, unit = "pt")),
      legend.text = element_text(size = 5, margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "pt"), vjust = 0.5),
      legend.position = "top",
      legend.justification = "left",
      legend.key.size = unit(2, "pt"),
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
      panel.grid = element_blank(),            
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"), 
      axis.text.x.bottom = element_text(size = 5, color = "black", margin = margin(t = 1, r = 0, b = 0, l = 0)),
      axis.text.y.left = element_text(size = 5, color = "black", margin = margin(t = 0, r = 1, b = 0, l = 0)),
      axis.title.x = element_text(size = 6, margin = margin(t = 2, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size = 6, margin = margin(t = 0, r = 2, b = 0, l = 0))) +
    scale_x_continuous(limits = c(0, 87), breaks = seq(0, 84, 12), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.0, 1), breaks = seq(0.0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
    scale_color_manual(values = scale_values,
                       labels = scale_labels) + 
    guides(color = guide_legend(keywidth = 0.4, keyheight = 0, nrow = 1, byrow = TRUE))
  
  filename = paste("KMPlot - ", matrix_filename, ".png", sep = "")
  ggsave(filename, plot = print(survplot), width = 4.1/.cm, height = 4.5/.cm, dpi = 1200, grDevices::png)
}
