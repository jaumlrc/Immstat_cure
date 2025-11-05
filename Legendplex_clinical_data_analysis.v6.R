setwd("/mnt/DATA/Joao_cunha/Paul_Kaye/EDCTP_citometry_markers/Combining_all_data/Script_to_paper")

##########################################################################
# This script does all the data cleaning and analysis for the manuscript #
# Systemic inflammatory markers of VL treatment response in East Africa  #
##########################################################################
# It receives as input the csv tables with:                              #
# 1-Patient clinical metadata                                            #
# 2-Legendplex results and Limit of Quantification                       #
##########################################################################
# Script tasks:                                                          #
# Generates several quality control cleaning                             #
# Separates the data in male and female patients                         #
# PCA and UMAPs                                                          #
# Compares the Pre and post treatment data                               #
# Identify potential makers linked tp clinical states                    #
##########################################################################

# Expected input data:
# Legendplex immune panel data, as in:
#well	experiment	sample	dilution	replicate	sample_type	TGF.β1..Free.Active...A4.	PAI.1..A5.	sTREM.1..A6.	PTX3..A7.	sCD40L..A8.	sCD25..IL.2Ra...A10.	CXCL12..B2.	sST2..B3.	sTNF.RI..B4.	sTNF.RII..B5.	sRAGE..B6.	CX3CL1..B7.	sCD130..gp130...B9.	Sample_name	group
#02-Well-A9.fcs	In_K_plate 3	02-Well-A9	1	1	Sample	19.23	20519.17	37.7	2129.05	18726.26	3469.37	819.34	1561	1182.25	438.21	298.15	9826.67	46440.25	VL024	V2
#02-Well-B10.fcs	In_K_plate 3	02-Well-B10	1	1	Sample	7.7	10330.28	32.56	1030.78	4446.65	3230.9	612.27	934.23	1216.53	970.39	738.5	10578.61	112821.11	VL028	V2
#02-Well-B11.fcs	In_K_plate 3	02-Well-B11	1	1	Sample	0	8971.49	30.73	1651.99	3133.82	12095.17	2308.97	3255.08	3182.56	1923.62	644.1	12922.09	96981.67	VL032	V2
#02-Well-B3.fcs	In_K_plate 3	02-Well-B3	1	1	Sample	0	8262.93	149.84	2482.91	10084.11	5303.83	679.34	2478.51	1571.75	615.69	376.21	23816.19	69932.39	VL034	V1

#Legendplex LOQ values data, as in:
#Plate	Cytokine	LOD
#plate1	CX3CL1	2255.26
#plate1	CXCL12	216.18
#plate1	PAI.1	749.46
#plate1	PTX3	173.92
#plate1	sCD130	335535.43
#plate1	sCD25	2892.87
#plate1	sCD40L	150.10

#Clinical metadata as in (note: This example does not represetn real data):
# Patient	Class	Sex	Age	Height	Weight	Hepatomegaly_York	Splenomegaly_York	Auxiliar_Lymphnodes	SittingSystolicBloodPressure	SittingDiastoliccBloodPressure	Temperature	Pulse	Hemoglobin	WBCells	Neutrophil	Lymphocyte	Platelets	ALT_Results	Creatinine	Albumin	Total_Bilirubin	Aspirate_grade	Liver_size	Spleen_size
# A00A	V1	1	30	1.7	77.5	0	0	0	NA		32	99	10	7	2.73	1.24	101	47.8	32.6	2.8	5.5	2	NA	NA
# A00B	V1	1	25	1.8	85.4	0	1	0	144	78	34.6	100	7.8	3.62	1.6	1.63	93	20.1	98.8	1.49	104.5	1	NA	10


# Disclaimers.
#This script was specifically tailored to work with the legendplex and clinical data from the "Systemic inflammatory markers of VL treatment response in East Africa" manuscript.
# It will not work with other data unless modifications are made
# It was tested and run on R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 20.04.6 LTS
# Bioconductor version 3.20

#Silencing warnings
oldw <- getOption("warn")
options(warn = -1)
#options(warn = oldw)
#options(warn = 0)

#Loading the required libraries:
library(ggplot2)
library(reshape2)
library(viridis)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(factoextra)
library(ggrepel)
library(plotly) 
library(umap) 
library(pheatmap)
library(Hmisc)
library(corrplot)
library(ggradar)
library(fmsb)
library(scales)
library(cowplot)
library(ggiraphExtra)
library(grid)
library(gridExtra)
library(ggprism)
library(ggfortify)
library(cluster)
library(caret)
library(mixOmics)
library(sva)
library(logistf)
library(ggpointdensity)
library(MASS)
library(ggVennDiagram)
library(tibble)
library(ggtext)


#Creating the directories to separate the male and female data
dir.create("immstat_Male_data.dir")
dir.create("immstat_female_data.dir")

##############################################
#Custom libraries to be used in the script ###
##############################################
## This section loads all the required custom functions for the data analysis #
## These are all custim-made functions to do/authomatize the data analysis

#Used Functions:
#These functions correct the limit of detection of the Legendplex data
correct_measurements_limit <- function(input_table, limit_table, dilution_facotr_to_multiply) {
  temp_raw_table <- input_table
  temp_limit_table <- limit_table
  
  #remove non-existent columns from here
  temp_raw_table_melt <- melt(temp_raw_table, id.vars = c("Sample_name", "group", "plate", "class", "experiment", "well", "sample", "dilution", "replicate", "sample_type"))
  temp_raw_table_melt$variable <- gsub("\\.\\..*","",temp_raw_table_melt$variable)
  
  #Removing NAs from the sample names:
  temp_raw_table_melt <- temp_raw_table_melt[!is.na(temp_raw_table_melt$Sample_name),]
  
  #Correcting by the LOD:
  colnames(temp_limit_table) <- c("plate", "variable","LOD")
  merged_with_LOD <- left_join(temp_raw_table_melt, temp_limit_table, by=c("plate","variable"))
  merged_with_LOD$Raw_values <- merged_with_LOD$value
  
  a=1
  for (a in 1:nrow(merged_with_LOD)) {
    if(merged_with_LOD[a,"value"] <= merged_with_LOD[a,"LOD"]) {merged_with_LOD[a,"value"] <- merged_with_LOD[a,"LOD"]}
  }
  
  merged_with_LOD$value <- merged_with_LOD$value * dilution_facotr_to_multiply
  merged_with_LOD$dilution <- dilution_facotr_to_multiply
  
  return(merged_with_LOD)
  
}
#Uganda has a different custom function as some samples were diluted and others not. This was accounted for here
correct_measurements_limit_UGANDA <- function(input_table, limit_table, dilution_facotr_to_multiply) {
  temp_raw_table <- input_table
  temp_limit_table <- limit_table
  
  #remove non-existent columns from here
  temp_raw_table_melt <- melt(temp_raw_table, id.vars = c("Sample_name", "group", "plate", "class", "experiment", "well", "sample", "dilution", "replicate", "sample_type"))
  temp_raw_table_melt$variable <- gsub("\\.\\..*","",temp_raw_table_melt$variable)
  
  #Removing NAs from the sample names:
  temp_raw_table_melt <- temp_raw_table_melt[!is.na(temp_raw_table_melt$Sample_name),]
  
  #Correcting by the LOD:
  colnames(temp_limit_table) <- c("plate", "variable","LOD")
  merged_with_LOD <- left_join(temp_raw_table_melt, temp_limit_table, by=c("plate","variable"))
  # nrow(merged_with_LOD)
  # nrow(temp_raw_table_melt)
  merged_with_LOD$Raw_values <- merged_with_LOD$value
  
  a=1
  for (a in 1:nrow(merged_with_LOD)) {
    if(merged_with_LOD[a,"value"] <= merged_with_LOD[a,"LOD"]) {merged_with_LOD[a,"value"] <- merged_with_LOD[a,"LOD"]}
  }
  
  merged_with_LOD$dilution[merged_with_LOD$Sample_name %in% c("VL001","VL002","VL003","VL004","VL005","VL006")] = 2
  merged_with_LOD$value <- merged_with_LOD$value * merged_with_LOD$dilution
  
  return(merged_with_LOD)
  
}

#Establish and plot the limit of detection of each plate
plot_limit_of_detection <- function(input_full_table, output_prefix) {
  #first 
  to_limit_of_detection <- input_full_table
  temp_name_LOD <- output_prefix
  
  to_limit_of_detection$Class <- "Lower"
  to_limit_of_detection$Class[to_limit_of_detection$Raw_values>=to_limit_of_detection$LOD] <- "Higher"
  to_limit_of_detection$sample <- gsub(".*-","",to_limit_of_detection$sample)
  to_limit_of_detection$ID <- paste(to_limit_of_detection$Sample_name, to_limit_of_detection$group, to_limit_of_detection$sample, sep = "_")
  
  #Separating in samples, HV and pools:
  to_limit_of_detection_samples <- to_limit_of_detection[to_limit_of_detection$group %in% c("V1","V2"),]
  to_limit_of_detection_control <- to_limit_of_detection[to_limit_of_detection$group %in% c("HV_Pool","VL_Pool"),]
  to_limit_of_detection_hv <- to_limit_of_detection[to_limit_of_detection$group =="HV",]
  
  #Selecting the relevant columns and changing to 0 and 1:
  to_limit_of_detection_samples2 <- to_limit_of_detection_samples[,c("ID", "Class", "variable")]
  to_limit_of_detection_control2 <- to_limit_of_detection_control[,c("ID", "Class", "variable")]
  to_limit_of_detection_hv2 <- to_limit_of_detection_hv[,c("ID", "Class", "variable")]
  
  to_limit_of_detection_samples2[to_limit_of_detection_samples2=="Lower"] <- 0
  to_limit_of_detection_samples2[to_limit_of_detection_samples2=="Higher"] <- 1
  
  to_limit_of_detection_control2[to_limit_of_detection_control2=="Lower"] <- 0
  to_limit_of_detection_control2[to_limit_of_detection_control2=="Higher"] <- 1
  
  to_limit_of_detection_hv2[to_limit_of_detection_hv2=="Lower"] <- 0
  to_limit_of_detection_hv2[to_limit_of_detection_hv2=="Higher"] <- 1
  
  #Transforming the data frame from long to wide:
  to_limit_of_detection_samples2_heatmap <- data.frame(to_limit_of_detection_samples2 %>% pivot_wider(names_from = variable, values_from = Class))
  to_limit_of_detection_hv2_heatmap <- data.frame(to_limit_of_detection_hv2 %>% pivot_wider(names_from = variable, values_from = Class))
  
  #Patients
  rownames(to_limit_of_detection_samples2_heatmap) <- to_limit_of_detection_samples2_heatmap$ID
  to_limit_of_detection_samples2_heatmap <- to_limit_of_detection_samples2_heatmap[order(to_limit_of_detection_samples2_heatmap$ID),]
  to_limit_of_detection_samples2_heatmap <- to_limit_of_detection_samples2_heatmap[,colnames(to_limit_of_detection_samples2_heatmap) != "ID"]
  to_limit_of_detection_samples2_heatmap <- mutate_all(to_limit_of_detection_samples2_heatmap, function(x) as.numeric(as.character(x)))
  
  #HV
  rownames(to_limit_of_detection_hv2_heatmap) <- to_limit_of_detection_hv2_heatmap$ID
  to_limit_of_detection_hv2_heatmap <- to_limit_of_detection_hv2_heatmap[order(to_limit_of_detection_hv2_heatmap$ID),]
  to_limit_of_detection_hv2_heatmap <- to_limit_of_detection_hv2_heatmap[,colnames(to_limit_of_detection_hv2_heatmap) != "ID"]
  to_limit_of_detection_hv2_heatmap <- mutate_all(to_limit_of_detection_hv2_heatmap, function(x) as.numeric(as.character(x)))
  
  #Setting the colorstrips
  tocol_column <- data.frame(rownames(to_limit_of_detection_samples2_heatmap))
  colnames(tocol_column) <- "Ids"
  tocol_column$class <- NA
  tocol_column$class[grepl("HV", tocol_column$Ids)] <- "HV"
  tocol_column$class[grepl("V1", tocol_column$Ids)] <- "V1"
  tocol_column$class[grepl("V2", tocol_column$Ids)] <- "V2"
  rownames(tocol_column) <- tocol_column$Ids
  tocol_column <- tocol_column[,colnames(tocol_column) !="Ids", drop=FALSE]
  my_colour = list(
    class = c(HV = "#CC79A7", V1 = "#E69F00", V2 ="#56B4E9")
  )
  
  png(paste(temp_name_LOD, "samples_heatmap_LOD.png", sep = "_"), width = 1500, height = 4000, res = 300)
  print(pheatmap(to_limit_of_detection_samples2_heatmap, fontsize_row =3 , color=c("#DC143C", "white"),
                 annotation_row =tocol_column, annotation_colors=my_colour, cluster_cols = FALSE, cluster_rows = FALSE ))
  dev.off()
  
  
  png(paste(temp_name_LOD, "HV_heatmap_LOD.png", sep = "_"), width = 1500, height = 4000, res = 300)
  print(pheatmap(to_limit_of_detection_hv2_heatmap, fontsize_row =5 , color=c("#DC143C", "white"),
                 cluster_cols = FALSE, cluster_rows = FALSE ))
  dev.off()
  
}

#Standardize the different plates based on the control Pool, represented by VL patient data
standardizing_plates <- function(input_full_table, output_prefix) {
  
  #Obtaining the mean of the replicates from the pools
  controls_mean <- data.frame(input_full_table[input_full_table$Sample_name %in% c("HV_Pool","VL_Pool"),] 
                              %>% group_by(variable, Sample_name, plate) %>% summarise(mean_replicates = mean(value)))
  
  #Organizing the table and obtaining the mediam value of all plates
  controls_mean_pivoted <- data.frame(pivot_wider(controls_mean, names_from = plate, values_from = mean_replicates))
  controls_mean_pivoted$median_norm <- apply(controls_mean_pivoted[,3:ncol(controls_mean_pivoted)], 1, FUN = median)
  
  #Generating the normalizer for each plate
  plates_ids <- unique(controls_mean$plate)
  
  #Setting the normalizer for each plate
  a=1
  for (a in 1:length(plates_ids)) {
    tempid <- plates_ids[a]
    temp_df_ratio <- data.frame(controls_mean_pivoted[, tempid]/controls_mean_pivoted$median_norm)
    temp_df_ratio$variable <- controls_mean_pivoted$variable
    temp_df_ratio$Sample_name <- controls_mean_pivoted$Sample_name
    
    #Using only the normalizer from VL:
    temp_df_ratio_wider <- data.frame(pivot_wider(temp_df_ratio, names_from = Sample_name, values_from = controls_mean_pivoted...tempid..controls_mean_pivoted.median_norm ))
    temp_df_ratio_wider2 <- temp_df_ratio_wider[,c("variable", "VL_Pool")]
    colnames(temp_df_ratio_wider2) <- c("variable",paste(tempid,"normalizer", sep = "_"))
    temp_df_ratio_wider2_VL <- temp_df_ratio_wider2
    temp_df_ratio_wider2_HV <- temp_df_ratio_wider2
    temp_df_ratio_wider2_VL$Sample_name <- "VL_Pool"
    temp_df_ratio_wider2_HV$Sample_name <- "HV_Pool"
    temp_df_ratio_wider3 <- rbind(temp_df_ratio_wider2_VL,temp_df_ratio_wider2_HV)
    
    temp_df_ratio <- temp_df_ratio_wider3
    
    controls_mean_pivoted <- merge(controls_mean_pivoted, temp_df_ratio, by=c("variable","Sample_name"))
    temp_normalized <- data.frame(controls_mean_pivoted[,tempid]/controls_mean_pivoted[,paste(tempid,"normalizer", sep = "_")])
    colnames(temp_normalized) <- paste(tempid,"normalized", sep = "_")
    temp_normalized$variable <- controls_mean_pivoted$variable
    temp_normalized$Sample_name <- controls_mean_pivoted$Sample_name
    
    
    
    controls_mean_pivoted <- merge(controls_mean_pivoted, temp_normalized, by=c("variable","Sample_name"))
    
  }
  
  #Plotting the before and after of each plate normalizer
  
  controls_mean_pivoted_to_plot <- controls_mean_pivoted
  
  colnames_to_select <- colnames(controls_mean_pivoted_to_plot[,grepl("plate", colnames(controls_mean_pivoted_to_plot))])
  colnames_to_select <- c("variable","Sample_name",colnames_to_select)
  
  controls_mean_pivoted_to_plot_bn <- melt(controls_mean_pivoted_to_plot[,colnames_to_select], id.vars = c("variable","Sample_name"))
  colnames(controls_mean_pivoted_to_plot_bn)[1] <- "marker"
  
  controls_mean_pivoted_to_plot_bn2 <- data.frame(controls_mean_pivoted_to_plot_bn %>% group_by(marker,Sample_name, variable) %>% summarise(mean_value =mean(value)))
  #Removing the normalizer data 
  controls_mean_pivoted_to_plot_bn2 <- controls_mean_pivoted_to_plot_bn2[!grepl("normalizer",controls_mean_pivoted_to_plot_bn2$variable),]
  
  png(paste(output_prefix,"barplot_normalizer_LOD.png", sep = "_"), width = 3500,height = 3500, res = 300)
  print(ggplot(controls_mean_pivoted_to_plot_bn2, aes(x=variable, y=mean_value, fill=Sample_name, group=Sample_name)) + geom_col(alpha=0.5, position = "dodge2") +  
          facet_wrap(~marker, scales = "free") + theme_bw() + ylab("pg/mL") + scale_fill_manual(values = c("#CC79A7","#E69F00"), breaks = c("HV_Pool", "VL_Pool")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_log10())
  dev.off()
  
  
  #Generating the heatmap of the normalizers
  columns_to_select <- colnames(controls_mean_pivoted)[grepl("_normalizer", colnames(controls_mean_pivoted))]
  columns_to_select <- c("variable", "Sample_name", columns_to_select)
  controls_mean_pivoted_to_heatmap <- controls_mean_pivoted[,columns_to_select]
  rownames(controls_mean_pivoted_to_heatmap) <- paste(controls_mean_pivoted_to_heatmap$variable, controls_mean_pivoted_to_heatmap$Sample_name, sep = "_")
  controls_mean_pivoted_to_heatmap2 <- controls_mean_pivoted_to_heatmap[,!colnames(controls_mean_pivoted_to_heatmap) %in% c("variable","Sample_name")]
  
  #Setting the colours for the heatmap
  tocol_column <- data.frame(rownames(controls_mean_pivoted_to_heatmap2))
  colnames(tocol_column) <- "Ids"
  tocol_column$class <- NA
  tocol_column$class[grepl("HV", tocol_column$Ids)] <- "HV"
  tocol_column$class[grepl("VL", tocol_column$Ids)] <- "VL"
  rownames(tocol_column) <- tocol_column$Ids
  tocol_column <- tocol_column[,colnames(tocol_column) !="Ids", drop=FALSE]
  my_colour = list(
    class = c(HV = "#CC79A7", VL = "#E69F00"))
  
  png(paste(output_prefix,"heatmap_Standard_normalizers.png", sep = "_"), width = 1500,height = 4000, res = 300)
  print(pheatmap(controls_mean_pivoted_to_heatmap2, fontsize_row =5, cluster_rows=FALSE, cluster_cols = FALSE, scale="none",
                 annotation_row =tocol_column, annotation_colors=my_colour, display_numbers=TRUE))
  dev.off()
  
  write.table(controls_mean_pivoted, paste(output_prefix,"normalizers.csv", sep = "_"), sep = ",", row.names = FALSE, quote = FALSE)
  
  #Now to add the normalizing information in the table
  colnames(controls_mean_pivoted_to_heatmap)[colnames(controls_mean_pivoted_to_heatmap)=="variable"] <- "marker"
  controls_mean_pivoted_to_heatmap_to_merge1 <- melt(controls_mean_pivoted_to_heatmap, id.vars = c("marker","Sample_name"))
  controls_mean_pivoted_to_heatmap_to_merge1$variable <- gsub("_normalizer", "", controls_mean_pivoted_to_heatmap_to_merge1$variable)
  
  colnames(controls_mean_pivoted_to_heatmap_to_merge1) <- c("variable","normalizer_origin","plate","normalizer")
  
  #Setting a marker to separate HV as controls and VL as samples:
  input_full_table$sample <- NA
  input_full_table$sample[grepl("VL",input_full_table$group)] <- "sample"
  input_full_table$sample[grepl("V1",input_full_table$group)] <- "sample"
  input_full_table$sample[grepl("V2",input_full_table$group)] <- "sample"
  input_full_table$sample[grepl("HV",input_full_table$group)] <- "control"
  
  
  controls_mean_pivoted_to_heatmap_to_merge1$sample <- NA
  controls_mean_pivoted_to_heatmap_to_merge1$sample[grepl("HV_Pool",controls_mean_pivoted_to_heatmap_to_merge1$normalizer_origin)] <- "control"
  controls_mean_pivoted_to_heatmap_to_merge1$sample[grepl("VL_Pool",controls_mean_pivoted_to_heatmap_to_merge1$normalizer_origin)] <- "sample"
  
  input_full_table2 <- left_join(input_full_table, controls_mean_pivoted_to_heatmap_to_merge1, by=c("variable","sample","plate"))
  nrow(input_full_table)
  nrow(input_full_table2)
  
  input_full_table2$norm_value <- input_full_table2$value/input_full_table2$normalizer
  
  options(scipen=999)
  input_full_table2_no_controls <- input_full_table2[!input_full_table2$group %in% c("HV_Pool", "VL_Pool"),]
  png(paste(output_prefix,"Boxplot_plates_control.png", sep = "_"), width = 3000,height = 2500, res = 300)
  print(ggplot(input_full_table2_no_controls, aes(x=group, y=norm_value, fill=plate)) + geom_boxplot() + facet_wrap(~variable, scales = "free") +
          theme_bw() + scale_y_continuous(trans='log10')  + xlab("") + ylab("Normalized_values"))
  dev.off()
  
  png(paste(output_prefix,"Boxplot_plates_control_before_scaling.png", sep = "_"), width = 3000,height = 2500, res = 300)
  print(ggplot(input_full_table2_no_controls, aes(x=group, y=value, fill=plate)) + geom_boxplot() + facet_wrap(~variable, scales = "free") +
          theme_bw() + scale_y_continuous(trans='log10')  + xlab("") + ylab("Normalized_values"))
  dev.off()
  
  #Printing the normalizers:
  input_full_table2_no_controls_unique <- unique(input_full_table2_no_controls[,c("plate", "normalizer", "class", "variable")])
  
  
  png(paste(output_prefix,"boxplot_normalizer_LOD.png", sep = "_"), width = 3500,height = 3500, res = 300)
  print(ggplot(input_full_table2_no_controls_unique, aes(x=plate, y=normalizer, fill=class, group=class)) + geom_col(alpha=0.5, position = "dodge2") +  
          theme_bw() + ylab("Normalizer") + facet_wrap(~variable) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  return(input_full_table2)
  
  
}

#Assessing the duplicates variance in each plate
comparing_replicates <- function(input_full_table3, outprefix2) {
  
  summarized_table <- data.frame(input_full_table3 %>% group_by(Sample_name, group, variable) 
                                 %>% summarise(Mean_rep=mean(norm_value), SD_rep=sd(norm_value)))
  
  summarized_table$cv <- (summarized_table$SD_rep/summarized_table$Mean_rep)*100
  
  table(summarized_table$group)
  
  summarized_table_samples <- summarized_table[!summarized_table$group %in% c("HV_Pool", "VL_Pool"),]
  summarized_table_controls <-summarized_table[!summarized_table$group %in% c("HV_Pool", "VL_Pool"),]
  
  png(paste(outprefix2,"CV_error_samples.png", sep = "_"), width = 5000, height = 1300, res = 300)
  print(ggplot(summarized_table_samples, aes(x=Sample_name, y= cv, color = variable, shape=variable)) + geom_point() + 
          facet_wrap(~group, scales = "free_x") + theme_bw() + scale_shape_manual(values = rep(c(0,1,2,3,4,5),10)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  
}

# Function to estimate the missing data both in the clinical and Legendplex data
estimate_missing_data <- function(input_table, name){
  #Estimating the missing data:
  
  #Iterating by country:
  countries <- unique(input_table$Country)
  
  d=1
  for (d in 1:length(countries)) {
    country <- countries[d]
    
    metadata_table <-  input_table[which(input_table$Country==country),]
    
    #Estimating the proportion of missing data for each character:
    metadata_table_V1 <- metadata_table[which(metadata_table$Class == "V1"),]
    na_counts_base_V1 <- data.frame(lapply(metadata_table_V1, function(x) sum(is.na(x))), Group="V1")
    metadata_table_V2 <- metadata_table[which(metadata_table$Class == "V2"),]
    na_counts_base_V2 <- data.frame(lapply(metadata_table_V2, function(x) sum(is.na(x))), Group="V2")
    metadata_table_HV <- metadata_table[which(metadata_table$Class == "HV"),]
    na_counts_base_HV <- data.frame(lapply(metadata_table_HV, function(x) sum(is.na(x))), Group="HV")
    
    #Combining and transposing the DF:
    na_counts_all <- rbind(na_counts_base_V1, na_counts_base_V2, na_counts_base_HV)
    na_counts_allmelt <- melt(na_counts_all, id.vars = "Group")
    
    na_counts_allmelt$sample_number <- NA
    na_counts_allmelt$sample_number[which(na_counts_allmelt$Group=="V1")] <- nrow(metadata_table_V1)
    na_counts_allmelt$sample_number[which(na_counts_allmelt$Group=="V2")] <- nrow(metadata_table_V2)
    na_counts_allmelt$sample_number[which(na_counts_allmelt$Group=="HV")] <- nrow(metadata_table_HV)
    
    na_counts_allmelt$MissingPercentage <- (na_counts_allmelt$value/na_counts_allmelt$sample_number)*100
    
    png(paste(name, country, "missing_data.png", sep = "_"), width = 3800, height = 1500, res = 300)
    print(ggplot(na_counts_allmelt, aes(x=variable,y=MissingPercentage)) + geom_col() + facet_wrap(~Group, scales = "free_x") +
            theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,100))
    dev.off()
    
  }
  
  metadata_table_2 <- input_table
  #Separating in V1 and V2
  metadata_table_2_V1 <- metadata_table_2[which(metadata_table_2$Class=="V1"),]
  metadata_table_2_V2 <- metadata_table_2[which(metadata_table_2$Class=="V2"),]
  #Including a suffix for the colnames:
  colnames(metadata_table_2_V1)[3:ncol(metadata_table_2_V1)] <- paste(colnames(metadata_table_2_V1[3:ncol(metadata_table_2_V1)]), "_V1", sep = "")
  colnames(metadata_table_2_V2)[3:ncol(metadata_table_2_V2)] <- paste(colnames(metadata_table_2_V2[3:ncol(metadata_table_2_V2)]), "_V2", sep = "")
  metadata_table_2_V1$UID <- gsub("V1_","", metadata_table_2_V1$UID)
  metadata_table_2_V2$UID <- gsub("V2_","", metadata_table_2_V2$UID)
  
  metadata_table_2_V1_V2 <- full_join(metadata_table_2_V1, metadata_table_2_V2, by=c("UID", "Patient"))
  metadata_table_2_V1_V2 <- metadata_table_2_V1_V2[order(metadata_table_2_V1_V2$Country_V1),]
  metadata_table_2_V1_V2$Country_V1 <- as.factor(metadata_table_2_V1_V2$Country_V1) 
  
  rownames(metadata_table_2_V1_V2) <- metadata_table_2_V1_V2$UID
  metadata_table_2_V1_V2 <-  subset(metadata_table_2_V1_V2, select = -c(Patient, Class_V1, Country_V1, Class_V2, Country_V2, UID))
  metadata_table_2_V1_V2[!is.na(metadata_table_2_V1_V2)] <- 1
  metadata_table_2_V1_V2[is.na(metadata_table_2_V1_V2)] <- 0
  
  
  #Setting the colorstrips
  tocol_column <- data.frame(rownames(metadata_table_2_V1_V2))
  colnames(tocol_column) <- "Ids"
  tocol_column$class <- NA
  tocol_column$class[grepl("Ethiopia", tocol_column$Ids)] <- "Ethiopia"
  tocol_column$class[grepl("Kenya", tocol_column$Ids)] <- "Kenya"
  tocol_column$class[grepl("Uganda", tocol_column$Ids)] <- "Uganda"
  tocol_column$class[grepl("Sudan", tocol_column$Ids)] <- "Sudan"
  
  rownames(tocol_column) <- tocol_column$Ids
  tocol_column <- tocol_column[,colnames(tocol_column) !="Ids", drop=FALSE]
  my_colour = list(
    class = c(Ethiopia = "#F8766D", Kenya = "#00BFC4", Uganda ="#C77CFF", Sudan= "#7CAE00")
  )
  
  
  pdf(paste(name, "missing_data_heatmap.pdf", sep = "_"), width = 18, height = 12)
  print(pheatmap(metadata_table_2_V1_V2, fontsize_row =3 , color=c("#DC143C", "white"),
                 cluster_cols = TRUE, cluster_rows = FALSE, annotation_row =tocol_column, annotation_colors=my_colour))
  
  dev.off()
  
}

# Plotting the evaluation of the continuous and discrete data
plot_continuous_and_disctreta_data <- function(clinical_data_traits, output_prefix){
  
  output_prefix_name <- output_prefix
  clinical_data_traits_melt <- melt(clinical_data_traits, id.vars = c("UID", "Patient", "Class",  "Country"))
  
  relevant_data_continuous <- c("Patient", "Class", "Age", "Height", "Weight", "SSystolicBP", "Spleen_size", "Liver_size",
                                "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin",  "WBCells", "Neutrophil", "Aspirate_grade",
                                "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin", "Country", "Spleen_size")
  
  to_continuous_data <- clinical_data_traits_melt[clinical_data_traits_melt$variable %in% relevant_data_continuous,]
  
  to_continuous_data_summary <- data.frame(to_continuous_data %>% group_by(Class, Country, variable) %>% 
                                             filter(!is.na(value)) %>%
                                             summarise(Median_value = median(value), Q25=quantile(value,probs = 0.25),
                                                       Q75=quantile(value,probs = 0.75), SD=sd(value), Count=n()))
  
  
  to_continuous_data_summary$IQR <- to_continuous_data_summary$Q75-to_continuous_data_summary$Q25
  to_continuous_data_summary$IQR_median <- to_continuous_data_summary$IQR/to_continuous_data_summary$Median_value
  
  to_continuous_data_summary$variable <- factor(to_continuous_data_summary$variable, levels = c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                                                                                                "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                                                                                                "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                                                "Lymphocyte"))
  png(paste(output_prefix_name, "Overall_continuous_clinical_data.png", sep = "_"), width = 3000, height = 3000, res=300)
  print(ggplot(to_continuous_data_summary, aes(y=Median_value, x=Class, color=Country, group=Country)) + 
          geom_line() + facet_wrap(~variable, scales = "free") + theme_bw() + expand_limits(y=0) + ylab("Median values") + xlab("Patient stage") +
          scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), breaks = c("Ethiopia", "Kenya", "Sudan", "Uganda"))+
          geom_pointrange(aes(ymin=Q25, ymax=Q75, alpha=0.5)) + scale_alpha(guide = 'none'))
  dev.off()
  
  #Removing the line from the first point
  significant_plot_all_melt_summary2 <- to_continuous_data_summary
  significant_plot_all_melt_summary2$Countries <- NA
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Ethiopia" &
                                                 significant_plot_all_melt_summary2$Class=="HV" ] <- "EthiopiaHV"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Ethiopia" &
                                                 significant_plot_all_melt_summary2$Class=="V1" ] <- "EthiopiaP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Ethiopia" &
                                                 significant_plot_all_melt_summary2$Class=="V2" ] <- "EthiopiaP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Kenya" &
                                                 significant_plot_all_melt_summary2$Class=="HV" ] <- "KenyaHV"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Kenya" &
                                                 significant_plot_all_melt_summary2$Class=="V1" ] <- "KenyaP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Kenya" &
                                                 significant_plot_all_melt_summary2$Class=="V2" ] <- "KenyaP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Uganda" &
                                                 significant_plot_all_melt_summary2$Class=="HV" ] <- "UgandaHV"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Uganda" &
                                                 significant_plot_all_melt_summary2$Class=="V1" ] <- "UgandaP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Uganda" &
                                                 significant_plot_all_melt_summary2$Class=="V2" ] <- "UgandaP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Sudan" &
                                                 significant_plot_all_melt_summary2$Class=="HV" ] <- "SudanHV"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Sudan" &
                                                 significant_plot_all_melt_summary2$Class=="V1" ] <- "SudanP"
  significant_plot_all_melt_summary2$Countries[significant_plot_all_melt_summary2$Country =="Sudan" &
                                                 significant_plot_all_melt_summary2$Class=="V2" ] <- "SudanP"
  
  
  significant_plot_all_melt_summary2$variable <- factor(significant_plot_all_melt_summary2$variable, levels = c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                                                                                                                "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                                                                                                                "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                                                                "Lymphocyte"))
  
  
  png(paste(output_prefix_name,"Overall_continuous_clinical_data.dots_nolink.png", sep = "_"), width = 3000, height = 3000, res=300)
  print(ggplot(significant_plot_all_melt_summary2, aes(y=Median_value, x=Class, color=Country, group=Countries)) +
          geom_line() +
          facet_wrap(~variable, scales = "free") + theme_bw() + expand_limits(y=0) + ylab("Median values") + xlab("Patient stage") +
          scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), breaks = c("Ethiopia", "Kenya", "Sudan", "Uganda"))+
          geom_pointrange(aes(ymin=Q25, ymax=Q75, alpha=0.5), size=0.8) + scale_alpha(guide = 'none') +
          theme(strip.text = element_text(size = 10), legend.key.size = unit(2, 'cm'),
                axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
          geom_text_repel(aes(label=Count), size=3))
  dev.off()
  
  
  png(paste(output_prefix_name,"Overall_continuous_clinical_data.dots_nolink_no_number.png", sep = "_"), width = 3000, height = 3000, res=300)
  print(ggplot(significant_plot_all_melt_summary2, aes(y=Median_value, x=Class, color=Country, group=Countries)) +
          geom_line() +
          facet_wrap(~variable, scales = "free") + theme_bw() + expand_limits(y=0) + ylab("Median values") + xlab("Patient stage") +
          scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), breaks = c("Ethiopia", "Kenya", "Sudan", "Uganda"))+
          geom_pointrange(aes(ymin=Q25, ymax=Q75, alpha=0.5), size=0.8) + scale_alpha(guide = 'none') +
          theme(strip.text = element_text(size = 10), legend.key.size = unit(2, 'cm'),
                axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)))
  dev.off()
  
  
  #Now  for the discrete data
    relevant_data_discrete <- c("Patient", "Country", "Class", "Males","Hepatomegaly","Splenomegaly","Auxiliar_Lymphnodes")
  
  
  significant_plot_categorical_all <- clinical_data_traits_melt[clinical_data_traits_melt$variable %in% relevant_data_discrete,]
  significant_plot_categorical_all2 <- significant_plot_categorical_all[!is.na(significant_plot_categorical_all$value),]
  significant_plot_categorical_all2$value <- as.numeric(as.character(significant_plot_categorical_all2$value))
  significant_plot_categorical_all2_Binary2 <- data.frame(significant_plot_categorical_all2 %>% 
                                                            group_by(Class, Country, variable) %>% summarise(Sum_values =sum(value), Count=n()))
  significant_plot_categorical_all2_Binary2$Proportion <- significant_plot_categorical_all2_Binary2$Sum_values/significant_plot_categorical_all2_Binary2$Count
  significant_plot_categorical_all2_Binary2$Ratio <- paste(significant_plot_categorical_all2_Binary2$Sum_values, significant_plot_categorical_all2_Binary2$Count, sep = "/")
  
  
  png(paste(output_prefix_name,"Overall_discrete_clinical_data.png", sep = "_"), width = 3000, height = 3000, res=300)
  print(ggplot(significant_plot_categorical_all2_Binary2, aes(x=Class, y=Proportion, fill=Country)) + geom_col(alpha=0.5, position = "dodge", stat = "identity") + 
          facet_grid(variable ~ Country) + ylim(0,1) + theme_bw()  + 
          scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), breaks = c("Ethiopia", "Kenya", "Sudan", "Uganda")) +
          theme(strip.text = element_text(size = 15)) +
          geom_text(aes(label=Ratio), position = position_dodge(width=1),  size=3))
  dev.off()
  
  
  significant_plot_categorical_all2_Binary2_noMales <- significant_plot_categorical_all2_Binary2[which(significant_plot_categorical_all2_Binary2$variable != "Males"),]
  
  
  significant_plot_categorical_all2_Binary2_noMales <- significant_plot_categorical_all2_Binary2[significant_plot_categorical_all2_Binary2$variable != "Males",]
  png(paste(output_prefix_name,"Overall_discrete_clinical_data_noMales.png", sep = "_"), width = 3000, height = 3000, res=300)
  print(ggplot(significant_plot_categorical_all2_Binary2_noMales, aes(x=Class, y=Proportion, fill=Country)) + geom_col(alpha=0.5, position = "dodge", stat = "identity") + 
          facet_grid(variable ~ Country) + ylim(0,1) + theme_bw()  + 
          scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"), breaks = c("Ethiopia", "Kenya", "Sudan", "Uganda")) +
          theme(strip.text = element_text(size = 15)) +
          geom_text(aes(label=Ratio), position = position_dodge(width=1),  size=3))
  dev.off()
  
  
}

# Correlation plots. Including both the clinical metadata and the Legendplex data:
correlation_plots_with_metadata <- function(input_table, sample_order_to_plot, country, output_plots) {
  
  #Correlation plots:
  #http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
  #I canceled the rcorr. It was OK, but had not multiple test correction. I implemented the correlation analysis myself with base R.
  
  input_table[,5:ncol(input_table)] <- lapply(input_table[,5:ncol(input_table)], as.numeric)
  
  #Exporting the dataframes with the correlation values and p-values
  
  # Creating a directory and plotting the significant correlations:
  dir.create(file.path(paste(output_plots,"corplot_individual.dir", sep = ".")), showWarnings = FALSE)
  dir.create(file.path(paste(paste(output_plots,"corplot_individual.dir", sep = ".")), "Significant", sep=""), showWarnings = FALSE)
  dir.create(file.path(paste(paste(output_plots,"corplot_individual.dir", sep = ".")), "Not_significant", sep=""), showWarnings = FALSE)
  
  #Generating a dataframe with all combinations:
  # Generate all unique unordered pairs (no self-pairs, no duplicates)
  input_table_ids <- colnames(input_table)
  input_table_ids <- input_table_ids[!input_table_ids %in% c("UID", "Patient", "Class", "Country", "Males")]
  
  #Generating all the combinations of comparisons
  input_table_ids_allcomb <- as.data.frame(t(combn(input_table_ids, 2)))
  colnames(input_table_ids_allcomb) <- c("Var1", "Var2")
  
  #DF to save each iteration
  df_to_return_spearman <- data.frame()
  
  for (i in 1:nrow(input_table_ids_allcomb)) {
    paramter1 <- as.character(input_table_ids_allcomb[i,"Var1"])
    paramter2 <- as.character(input_table_ids_allcomb[i,"Var2"])
    
    input_table3 <- input_table[,c("UID", paramter1,paramter2)]
    input_table3 <- input_table3[!is.na(input_table3[,paramter1]), ]
    input_table3 <- input_table3[!is.na(input_table3[,paramter2]), ]
    
    #If check if the values are not NA and if there are at least three values
    if (!all(is.na(input_table3[,paramter1])) && !all(is.na(input_table3[,paramter2]))
        && nrow(input_table3) >3) {
      
      temp_correlation <- cor.test(as.numeric(as.character(input_table3[,paramter1])), as.numeric(as.character(input_table3[,paramter2])), method = "spearman")
      temp_correlation_p <- format(round(temp_correlation$p.value, 4), nsmall = 2)
      temp_correlation_cor <- format(round(temp_correlation$estimate, 4), nsmall = 2)
      
      temp_iteration_sperman <- data.frame(P1 = paramter1, P2 =paramter2, Pvalue=temp_correlation_p, Cor=temp_correlation_cor)
      rownames(temp_iteration_sperman) <- ""
      
      df_to_return_spearman <- rbind(df_to_return_spearman, temp_iteration_sperman)
      
      
      if (temp_correlation_p <= 0.05) {
        pdf(paste(paste(output_plots,"corplot_individual.dir", sep = "."),"/Significant/",paste(paramter1,paramter2,"corplot.pdf", sep = "."), sep = ""), height = 6, width = 6)
        print(ggplot(input_table3, aes_string(x=paramter1,y=paramter2)) + geom_point(alpha=0.7, shape=1) + theme_bw() + theme(plot.title = element_text(size = 8)) +
                geom_smooth(method = "lm", se = TRUE, color = "blue") +
                ggtitle(paste(paramter1,paramter2,"Pvalue:",temp_correlation_p,"rho:",temp_correlation_cor, sep = "_")))
        dev.off()
      }
      
      else {
        pdf(paste(paste(output_plots,"corplot_individual.dir", sep = "."),"/Not_significant/",paste(paramter1,paramter2,"corplot.pdf", sep = "."), sep = ""), height = 6, width = 6)
        print(ggplot(input_table3, aes_string(x=paramter1,y=paramter2)) + geom_point(alpha=0.7, shape=1) + theme_bw() + theme(plot.title = element_text(size = 8)) +
                geom_smooth(method = "lm", se = TRUE, color = "blue")+
                ggtitle(paste(paramter1,paramter2,"Pvalue:",temp_correlation_p,"rho:",temp_correlation_cor, sep = "_")))
        dev.off()
      }
    } 
    
    else { temp_iteration_sperman <- data.frame(P1 = paramter1, P2 = paramter2, Pvalue= NA, Cor= NA)
    rownames(temp_iteration_sperman) <- ""
    
    df_to_return_spearman <- rbind(df_to_return_spearman, temp_iteration_sperman)
    }
    
  }  
  
  
  #Generating the plots corrected by multiple testing
  df_to_return_spearman_2 <- df_to_return_spearman
  df_to_return_spearman_2_temp <- df_to_return_spearman_2[,c("P2", "P1", "Pvalue","Cor")]
  colnames(df_to_return_spearman_2_temp) <- c("P1", "P2", "Pvalue","Cor")
  
  df_to_return_spearman_2$P_adjust <- p.adjust(df_to_return_spearman_2$Pvalue, method = "BH")
  df_to_return_spearman_2_temp$P_adjust <- p.adjust(df_to_return_spearman_2_temp$Pvalue, method = "BH")
  
  #Binding both to do the plot
  df_to_return_spearman_3 <- rbind(df_to_return_spearman_2, df_to_return_spearman_2_temp)
  df_to_return_spearman_3[which(df_to_return_spearman_3$P_adjust <0.05),]
  
  #Corr_data
  df_to_return_spearman_2_corr_matrix <- df_to_return_spearman_3[,c("P1", "P2", "Cor")] %>%
    pivot_wider(names_from = P2, values_from = Cor) %>%
    column_to_rownames("P1")
  
  #Replacing NA for 0 to plot
  df_to_return_spearman_2_corr_matrix[df_to_return_spearman_2_corr_matrix =="NA"] <- NA
  df_to_return_spearman_2_corr_matrix[is.na(df_to_return_spearman_2_corr_matrix)] <- 0
  
  nrow(df_to_return_spearman_2_corr_matrix)
  ncol(df_to_return_spearman_2_corr_matrix)
  
  # Convert values to numeric and build matrix
  df_to_return_spearman_2_corr_matrix2 <- as.matrix(sapply(df_to_return_spearman_2_corr_matrix, as.numeric))
  rownames(df_to_return_spearman_2_corr_matrix2) <- rownames(df_to_return_spearman_2_corr_matrix)
  
  # p-value adjusted
  df_to_return_spearman_2_p_matrix <- df_to_return_spearman_3[,c("P1", "P2", "P_adjust")] %>%
    pivot_wider(names_from = P2, values_from = P_adjust) %>%
    column_to_rownames("P1") 
  df_to_return_spearman_2_p_matrix[df_to_return_spearman_2_p_matrix =="NA"] <- NA
  df_to_return_spearman_2_p_matrix[is.na(df_to_return_spearman_2_p_matrix)] <- 1000
  
  # Convert values to numeric and build matrix
  df_to_return_spearman_2_p_matrix2 <- as.matrix(sapply(df_to_return_spearman_2_p_matrix, as.numeric))
  rownames(df_to_return_spearman_2_p_matrix2) <- rownames(df_to_return_spearman_2_p_matrix)
  
  #Pvalue non-adjusted:
  df_to_return_spearman_2_p_matrix_nonadjust <- df_to_return_spearman_3[,c("P1", "P2", "Pvalue")] %>%
    pivot_wider(names_from = P2, values_from = Pvalue) %>%
    column_to_rownames("P1") 
  df_to_return_spearman_2_p_matrix_nonadjust[df_to_return_spearman_2_p_matrix_nonadjust =="NA"] <- NA
  df_to_return_spearman_2_p_matrix_nonadjust[is.na(df_to_return_spearman_2_p_matrix_nonadjust)] <- 1000
  
  # Convert values to numeric and build matrix
  df_to_return_spearman_2_p_matrix_nonadjust2 <- as.matrix(sapply(df_to_return_spearman_2_p_matrix_nonadjust, as.numeric))
  rownames(df_to_return_spearman_2_p_matrix_nonadjust2) <- rownames(df_to_return_spearman_2_p_matrix_nonadjust)
  
  #Selecting the order of the samples:
  sample_order_to_plot <- sample_order_to_plot
  
  df_to_return_spearman_2_corr_matrix2_ordered <- df_to_return_spearman_2_corr_matrix2[sample_order_to_plot,sample_order_to_plot]
  df_to_return_spearman_2_p_matrix2_ordered <- df_to_return_spearman_2_p_matrix2[sample_order_to_plot,sample_order_to_plot]
  df_to_return_spearman_2_p_matrix_nonadjust2_oredered <- df_to_return_spearman_2_p_matrix_nonadjust2[sample_order_to_plot,sample_order_to_plot]
  
  
  pdf(paste(output_plots,"corplot_FDR_correct.pdf", sep = "."), height = 9, width = 18)
  corrplot(df_to_return_spearman_2_corr_matrix2_ordered, type = "upper", 
           diag = FALSE, col=colorRampPalette(c("#191970","white","#B22222"))(500),
           p.mat = df_to_return_spearman_2_p_matrix2_ordered, sig.level = 0.05, insig = "blank", addrect = 2, tl.col="black",
           tl.cex = 0.8)
  dev.off()
  
  pdf(paste(output_plots,"corplot_Not_correct.pdf", sep = "."), height = 9, width = 18)
  corrplot(df_to_return_spearman_2_corr_matrix2_ordered, type = "upper", 
           diag = FALSE, col=colorRampPalette(c("#191970","white","#B22222"))(500),
           p.mat = df_to_return_spearman_2_p_matrix_nonadjust2_oredered, sig.level = 0.05, insig = "blank", addrect = 2, tl.col="black",
           tl.cex = 0.8)
  dev.off()
  
  pdf(paste(output_plots,"corplot_No_pvalue.pdf", sep = "."), height = 9, width = 18)
  corrplot(df_to_return_spearman_2_corr_matrix2_ordered, type = "upper", 
           diag = FALSE, col=colorRampPalette(c("#191970","white","#B22222"))(500),
           p.mat = df_to_return_spearman_2_p_matrix2_ordered, sig.level = 1, insig = "blank", addrect = 2, tl.col="black",
           tl.cex = 0.8)
  dev.off()
  
  df_to_return_spearman_2$country <- country
  
  #Writing the output table:
  write.csv(df_to_return_spearman_2, paste(output_plots,"pearson_lognorm_stats_pvalue.csv", sep = "."), row.names = F, quote = F )
  
  names(df_to_return_spearman_2)[names(df_to_return_spearman_2) %in% c("Pvalue", "Cor", "P_adjust")] <- 
    c(paste(country, "Pvalue", sep = "_"), paste(country, "Cor", sep = "_"), paste(country, "P_adjust", sep = "_"))  
  
  return(df_to_return_spearman_2)
  
}

# Comparing clinical traits with reference values, to estimate the number and proportion of patients with values outside the reference range:
compare_clinical_results_to_referece <- function(input_table, output_prefix_dir) {
  #According to the project protocol, the reference values were:
  #Albumin 3 to 5g
  # Haemoglobin in men 11.0 – 16.0
  # Platelets adults 150–450
  # Creatinine <1.4 mg/dL
  # Total_Billirubin <1 mg/dL
  # ALT_results  < 40 U/L
  # WBCells 2.5 – 10 × 109/L
  # Neutrophil 1 – 6.6 × 10⁹/L 
  # Lymphocyte 0.66 – 4.4 × 10⁹/L
  
  
  input_table2 <- input_table[ 
    , c("Patient", "Class", "Country", "Hemoglobin", "WBCells",
        "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine" ,
        "Albumin", "Total_Bilirubin" )]
  
  
  #quantifying the number of patients from each country that are above each metric:
  input_table2_to_count <- input_table2
  #Albumin
  input_table2_to_count$Albumin_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Albumin_levels = ifelse(Albumin >= 3 & Albumin <= 5, "Normal", Albumin_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Albumin_levels = ifelse(Albumin < 3, "Low", Albumin_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Albumin_levels = ifelse(Albumin > 5, "High", Albumin_levels))
  
  #Haemoglobin
  input_table2_to_count$Haemoglobin_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Haemoglobin_levels = ifelse(Hemoglobin >= 11 & Hemoglobin <= 16, "Normal", Haemoglobin_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Haemoglobin_levels = ifelse(Hemoglobin < 11, "Low", Haemoglobin_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Haemoglobin_levels = ifelse(Hemoglobin > 16, "High", Haemoglobin_levels))
  
  #Platelets
  input_table2_to_count$Platelets_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Platelets_levels = ifelse(Platelets >= 150 & Platelets <= 450, "Normal", Platelets_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Platelets_levels = ifelse(Platelets < 150, "Low", Platelets_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Platelets_levels = ifelse(Platelets > 450, "High", Platelets_levels))
  
  #Creatinine
  input_table2_to_count$Creatinine_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Creatinine_levels = ifelse(Creatinine <=  1.4, "Normal", Creatinine_levels))
  #input_table2_to_count <- input_table2_to_count %>% mutate(Creatinine_levels = ifelse(Creatinine < 0.7, "Low", Creatinine_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Creatinine_levels = ifelse(Creatinine > 1.4, "High", Creatinine_levels))
  
  #Total_Billirubin
  input_table2_to_count$Total_Billirubin_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Total_Billirubin_levels = ifelse(Total_Bilirubin  <= 1, "Normal", Total_Billirubin_levels))
  #input_table2_to_count <- input_table2_to_count %>% mutate(Total_Billirubin_levels = ifelse(Total_Bilirubin < 0.1, "Low", Total_Billirubin_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Total_Billirubin_levels = ifelse(Total_Bilirubin > 1, "High", Total_Billirubin_levels))
  
  #ALT_results
  input_table2_to_count$ALT_results_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(ALT_results_levels = ifelse(ALT_Results <= 40, "Normal", ALT_results_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(ALT_results_levels = ifelse(ALT_Results > 40, "High", ALT_results_levels))
  
  #WBCells
  input_table2_to_count$WBCells_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(WBCells_levels = ifelse(WBCells >= 2.5 & WBCells <= 10, "Normal", WBCells_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(WBCells_levels = ifelse(WBCells < 2.5, "Low", WBCells_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(WBCells_levels = ifelse(WBCells > 10, "High", WBCells_levels))
  
  
  #Neutrophil
  input_table2_to_count$Neutrophil_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Neutrophil_levels = ifelse(Neutrophil >= 1 & Neutrophil <= 6.6, "Normal", Neutrophil_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Neutrophil_levels = ifelse(Neutrophil < 1, "Low", Neutrophil_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Neutrophil_levels = ifelse(Neutrophil > 6.6, "High", Neutrophil_levels))
  
  
  #Lymphocyte
  input_table2_to_count$Lymphocyte_levels <- NA
  input_table2_to_count <- input_table2_to_count %>% mutate(Lymphocyte_levels = ifelse(Lymphocyte >= 0.66 & Lymphocyte <= 4.4, "Normal", Lymphocyte_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Lymphocyte_levels = ifelse(Lymphocyte < 0.66, "Low", Lymphocyte_levels))
  input_table2_to_count <- input_table2_to_count %>% mutate(Lymphocyte_levels = ifelse(Lymphocyte > 4, "High", Lymphocyte_levels))
  
  
  input_table2_to_count_classes <-
    input_table2_to_count[,c("Patient", "Class",  "Country", "Albumin_levels", "Haemoglobin_levels", "Platelets_levels", "Creatinine_levels",
                             "Total_Billirubin_levels", "ALT_results_levels", "WBCells_levels", "Neutrophil_levels", "Lymphocyte_levels")]
  input_table2_to_count_classes <-
    input_table2_to_count_classes[order(input_table2_to_count_classes$Country),]
  
  
  input_table2_to_count_classes_melt <- melt(input_table2_to_count_classes, id.vars = c("Patient", "Class",  "Country"))
  
  
  input_table2_to_count_classes_melt$value[input_table2_to_count_classes_melt$value=="Low"] <- -1
  input_table2_to_count_classes_melt$value[input_table2_to_count_classes_melt$value=="Normal"] <- 0
  input_table2_to_count_classes_melt$value[input_table2_to_count_classes_melt$value=="High"] <- 1
  input_table2_to_count_classes_melt$UID <- paste(input_table2_to_count_classes_melt$Patient,
                                                  input_table2_to_count_classes_melt$Class,
                                                  input_table2_to_count_classes_melt$Country, sep = "_")
  
  normal_levels_recovery <- data.frame(input_table2_to_count_classes_melt %>% group_by(Class,  Country, variable, value) %>% summarise(Count=n()))
  normal_levels_recovery_wide <- data.frame(pivot_wider(normal_levels_recovery, names_from = "Class", values_from = "Count"))
  
  #Separating the column with HV:
  normal_levels_recovery_wide_HV <- normal_levels_recovery_wide[,c("Country", "variable", "value", "HV")]
  normal_levels_recovery_wide_HV2 <- data.frame(pivot_wider(normal_levels_recovery_wide_HV, names_from = "value", values_from = c("HV")))
  colnames(normal_levels_recovery_wide_HV2) <- gsub("X.1","HV_Low",colnames(normal_levels_recovery_wide_HV2))
  colnames(normal_levels_recovery_wide_HV2) <- gsub("X1","HV_High",colnames(normal_levels_recovery_wide_HV2))
  colnames(normal_levels_recovery_wide_HV2) <- gsub("X0","HV_Normal",colnames(normal_levels_recovery_wide_HV2))
  colnames(normal_levels_recovery_wide_HV2) <- gsub("NA.","HV_NA",colnames(normal_levels_recovery_wide_HV2))
  
  normal_levels_recovery_wide_HV2[is.na(normal_levels_recovery_wide_HV2)] <- 0
  
  normal_levels_recovery_wide_HV2$TotalHV <- normal_levels_recovery_wide_HV2$HV_Low + normal_levels_recovery_wide_HV2$HV_NA + normal_levels_recovery_wide_HV2$HV_Normal +
    normal_levels_recovery_wide_HV2$HV_High
  
  #Generating the percenntage WITHOUT NAs:
  
  normal_levels_recovery_wide_HV2$HV_low_Prop <- normal_levels_recovery_wide_HV2$HV_Low/normal_levels_recovery_wide_HV2$TotalHV 
  normal_levels_recovery_wide_HV2$HV_Normal_Prop <- normal_levels_recovery_wide_HV2$HV_Normal/normal_levels_recovery_wide_HV2$TotalHV 
  normal_levels_recovery_wide_HV2$HV_high_Prop <- normal_levels_recovery_wide_HV2$HV_High/normal_levels_recovery_wide_HV2$TotalHV 
  normal_levels_recovery_wide_HV2$HV_NA_Prop <- normal_levels_recovery_wide_HV2$HV_NA/normal_levels_recovery_wide_HV2$TotalHV 
  
  
  normal_levels_recovery_wide_HV2_to_plot <- normal_levels_recovery_wide_HV2[,c("Country","variable", "HV_low_Prop", "HV_Normal_Prop", "HV_high_Prop", "HV_NA_Prop")]
  normal_levels_recovery_wide_HV2_to_plot_melt <- melt(normal_levels_recovery_wide_HV2_to_plot, id.vars = c("Country","variable"))
  colnames(normal_levels_recovery_wide_HV2_to_plot_melt) <- c("Country","Marker", "Class", "Perecentage")
  normal_levels_recovery_wide_HV2_to_plot_melt$Marker <- gsub("_levels","",normal_levels_recovery_wide_HV2_to_plot_melt$Marker)
  normal_levels_recovery_wide_HV2_to_plot_melt$Class <- gsub("_Prop","",normal_levels_recovery_wide_HV2_to_plot_melt$Class)
  
  normal_levels_recovery_wide_HV2_to_plot_melt_split <- separate(normal_levels_recovery_wide_HV2_to_plot_melt, Class, into = c("Class", "Result"), sep = "_")
  normal_levels_recovery_wide_HV2_to_plot_melt_split$Country_Class <- paste(normal_levels_recovery_wide_HV2_to_plot_melt_split$Country,
                                                                            normal_levels_recovery_wide_HV2_to_plot_melt_split$Class, sep = "_")
  
  normal_levels_recovery_wide_HV2_to_plot_melt_split$Result <- factor(normal_levels_recovery_wide_HV2_to_plot_melt_split$Result, levels = c("high", "Normal", "low" ))
  
  
  write.table(normal_levels_recovery_wide_HV2_to_plot_melt_split, paste(output_prefix_dir, "Cdata_HV_NormLev_table.csv", sep = "."), row.names = F, quote = F, sep = ",")
  
  
  png(paste(output_prefix_dir, "percentages_with_NA_HV.png", sep = "."), width = 2000, height = 2000, res = 300)
  print(ggplot(normal_levels_recovery_wide_HV2_to_plot_melt_split, aes(x=Country_Class, y=Perecentage, fill=Result)) + geom_col() + facet_wrap(~Marker) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00",  "gray"), breaks = c("low","Normal","high", "NA")))
  dev.off()
  
  
  png(paste(output_prefix_dir, "percentages_with_NA_HV_1row.png", sep = "."), width = 4500, height = 1000, res = 300)
  print(ggplot(normal_levels_recovery_wide_HV2_to_plot_melt_split, aes(x=Country_Class, y=Perecentage, fill=Result)) + geom_col() + facet_wrap(~Marker, nrow = 1) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00",  "gray"), breaks = c("low","Normal","high", "NA")))
  dev.off()
  
  
  #Following without the HV for the downstream analysis:
  normal_levels_recovery_wide <- normal_levels_recovery_wide[,colnames(normal_levels_recovery_wide) != "HV"]
  
  normal_levels_recovery_wide2 <- data.frame(pivot_wider(normal_levels_recovery_wide, names_from = "value", values_from = c("V1", "V2")))
  colnames(normal_levels_recovery_wide2) <- gsub("_.1","_Low",colnames(normal_levels_recovery_wide2))
  colnames(normal_levels_recovery_wide2) <- gsub("_1","_High",colnames(normal_levels_recovery_wide2))
  colnames(normal_levels_recovery_wide2) <- gsub("_0","_Normal",colnames(normal_levels_recovery_wide2))
  normal_levels_recovery_wide2[is.na(normal_levels_recovery_wide2)] <- 0
  
  #Total Samples:
  normal_levels_recovery_wide2$TotalV1 <- normal_levels_recovery_wide2$V1_Low + normal_levels_recovery_wide2$V1_NA +normal_levels_recovery_wide2$V1_Normal +
    normal_levels_recovery_wide2$V1_High
  normal_levels_recovery_wide2$TotalV2 <- normal_levels_recovery_wide2$V2_Low + normal_levels_recovery_wide2$V2_NA +normal_levels_recovery_wide2$V2_Normal +
    normal_levels_recovery_wide2$V2_High
  normal_levels_recovery_wide2$TotalV1_noNA <- normal_levels_recovery_wide2$V1_Low + normal_levels_recovery_wide2$V1_Normal +
    normal_levels_recovery_wide2$V1_High
  normal_levels_recovery_wide2$TotalV2_noNA <- normal_levels_recovery_wide2$V2_Low + normal_levels_recovery_wide2$V2_Normal +
    normal_levels_recovery_wide2$V2_High
  
  #Generating the percenntage WITHOUT NAs:
  
  normal_levels_recovery_wide2$V1_low_Prop <- normal_levels_recovery_wide2$V1_Low/normal_levels_recovery_wide2$TotalV1_noNA 
  normal_levels_recovery_wide2$V1_Normal_Prop <- normal_levels_recovery_wide2$V1_Normal/normal_levels_recovery_wide2$TotalV1_noNA 
  normal_levels_recovery_wide2$V1_high_Prop <- normal_levels_recovery_wide2$V1_High/normal_levels_recovery_wide2$TotalV1_noNA 
  
  normal_levels_recovery_wide2$V2_low_Prop <- normal_levels_recovery_wide2$V2_Low/normal_levels_recovery_wide2$TotalV2_noNA 
  normal_levels_recovery_wide2$V2_Normal_Prop <- normal_levels_recovery_wide2$V2_Normal/normal_levels_recovery_wide2$TotalV2_noNA 
  normal_levels_recovery_wide2$V2_high_Prop <- normal_levels_recovery_wide2$V2_High/normal_levels_recovery_wide2$TotalV2_noNA 
  
  #And with NAs:
  normal_levels_recovery_wide2$V1_low_Prop_na <- normal_levels_recovery_wide2$V1_Low/normal_levels_recovery_wide2$TotalV1 
  normal_levels_recovery_wide2$V1_Normal_Prop_na <- normal_levels_recovery_wide2$V1_Normal/normal_levels_recovery_wide2$TotalV1 
  normal_levels_recovery_wide2$V1_high_Prop_na <- normal_levels_recovery_wide2$V1_High/normal_levels_recovery_wide2$TotalV1 
  normal_levels_recovery_wide2$V1_NA_Prop_na <- normal_levels_recovery_wide2$V1_NA/normal_levels_recovery_wide2$TotalV1 
  
  normal_levels_recovery_wide2$V2_low_Prop_na <- normal_levels_recovery_wide2$V2_Low/normal_levels_recovery_wide2$TotalV2 
  normal_levels_recovery_wide2$V2_Normal_Prop_na <- normal_levels_recovery_wide2$V2_Normal/normal_levels_recovery_wide2$TotalV2 
  normal_levels_recovery_wide2$V2_high_Prop_na <- normal_levels_recovery_wide2$V2_High/normal_levels_recovery_wide2$TotalV2 
  normal_levels_recovery_wide2$V2_NA_Prop_na <- normal_levels_recovery_wide2$V2_NA/normal_levels_recovery_wide2$TotalV2 
  
  write.table(normal_levels_recovery_wide2, paste(output_prefix_dir, "Cdata_nomlev_table.csv", sep = "."), row.names = F, quote = F, sep = ",")
  
  
  normal_levels_recovery_wide2_to_plot <- normal_levels_recovery_wide2[,c("Country","variable", "V1_low_Prop", "V1_Normal_Prop", "V1_high_Prop", "V2_low_Prop",
                                                                          "V2_Normal_Prop", "V2_high_Prop")]
  normal_levels_recovery_wide2_to_plot_melt <- melt(normal_levels_recovery_wide2_to_plot, id.vars = c("Country","variable"))
  colnames(normal_levels_recovery_wide2_to_plot_melt) <- c("Country","Marker", "Class", "Perecentage")
  normal_levels_recovery_wide2_to_plot_melt$Marker <- gsub("_levels","",normal_levels_recovery_wide2_to_plot_melt$Marker)
  normal_levels_recovery_wide2_to_plot_melt$Class <- gsub("_Prop","",normal_levels_recovery_wide2_to_plot_melt$Class)
  
  normal_levels_recovery_wide2_to_plot_melt_split <- separate(normal_levels_recovery_wide2_to_plot_melt, Class, into = c("Class", "Result"), sep = "_")
  normal_levels_recovery_wide2_to_plot_melt_split$Country_Class <- paste(normal_levels_recovery_wide2_to_plot_melt_split$Country,
                                                                         normal_levels_recovery_wide2_to_plot_melt_split$Class, sep = "_")
  
  normal_levels_recovery_wide2_to_plot_melt_split$Result <- factor(normal_levels_recovery_wide2_to_plot_melt_split$Result, levels = c("high", "Normal", "low" ))
  
  
  png(paste(output_prefix_dir, "percentages.png", sep = "."), width = 2000, height = 2000, res = 300)
  print(ggplot(normal_levels_recovery_wide2_to_plot_melt_split, aes(x=Country_Class, y=Perecentage, fill=Result)) + geom_col() + facet_wrap(~Marker) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00"), breaks = c("low","Normal","high")))
  dev.off()
  
  
  normal_levels_recovery_wide2_to_plot_na <- normal_levels_recovery_wide2[,c("Country","variable", "V1_low_Prop_na", "V1_Normal_Prop_na", "V1_high_Prop_na", "V1_NA_Prop_na",
                                                                             "V2_low_Prop_na", "V2_Normal_Prop_na", "V2_high_Prop_na", "V2_NA_Prop_na")]
  
  normal_levels_recovery_wide2_to_plot_na_melt <- melt(normal_levels_recovery_wide2_to_plot_na, id.vars = c("Country","variable"))
  colnames(normal_levels_recovery_wide2_to_plot_na_melt) <- c("Country","Marker", "Class", "Perecentage")
  normal_levels_recovery_wide2_to_plot_na_melt$Marker <- gsub("_levels","",normal_levels_recovery_wide2_to_plot_na_melt$Marker)
  normal_levels_recovery_wide2_to_plot_na_melt$Class <- gsub("_Prop","",normal_levels_recovery_wide2_to_plot_na_melt$Class)
  
  normal_levels_recovery_wide2_to_plot_na_meltt_split <- separate(normal_levels_recovery_wide2_to_plot_na_melt, Class, into = c("Class", "Result"), sep = "_")
  normal_levels_recovery_wide2_to_plot_na_meltt_split$Country_Class <- paste(normal_levels_recovery_wide2_to_plot_na_meltt_split$Country,
                                                                             normal_levels_recovery_wide2_to_plot_na_meltt_split$Class, sep = "_")
  
  normal_levels_recovery_wide2_to_plot_na_meltt_split$Result <- factor(normal_levels_recovery_wide2_to_plot_na_meltt_split$Result, levels = c("NA", "high", "Normal", "low" ))
  
  png(paste(output_prefix_dir, "percentages_with_NA.png", sep = "."), width = 2000, height = 2000, res = 300)
  print(ggplot(normal_levels_recovery_wide2_to_plot_na_meltt_split, aes(x=Country_Class, y=Perecentage, fill=Result)) + geom_col() + facet_wrap(~Marker) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00",  "gray"), breaks = c("low","Normal","high", "NA")))
  dev.off()
  
  png(paste(output_prefix_dir, "percentages_with_NA_1row.png", sep = "."), width = 4500, height = 1000, res = 300)
  print(ggplot(normal_levels_recovery_wide2_to_plot_na_meltt_split, aes(x=Country_Class, y=Perecentage, fill=Result)) + geom_col() + facet_wrap(~Marker, nrow = 1) +
          theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#56B4E9","#009E73","#E69F00",  "gray"), breaks = c("low","Normal","high", "NA")))
  dev.off()
  
  #to pheatmap:
  #Now, move the data to show V1 and V2 side by side
  input_table2_to_count_classes_temp_V1 <- 
    input_table2_to_count_classes[which(input_table2_to_count_classes$Class=="V1"),]
  input_table2_to_count_classes_temp_V2 <- 
    input_table2_to_count_classes[which(input_table2_to_count_classes$Class=="V2"),]
  
  
  colnames(input_table2_to_count_classes_temp_V1) <- paste(colnames(input_table2_to_count_classes_temp_V1), "_V1", sep = "")
  input_table2_to_count_classes_temp_V1$UID1 <- paste(input_table2_to_count_classes_temp_V1$Patient_V1,
                                                      input_table2_to_count_classes_temp_V1$Country_V1, sep = "_")
  
  
  colnames(input_table2_to_count_classes_temp_V2) <- paste(colnames(input_table2_to_count_classes_temp_V2), "_V2", sep = "")
  input_table2_to_count_classes_temp_V2$UID1 <- paste(input_table2_to_count_classes_temp_V2$Patient_V2,
                                                      input_table2_to_count_classes_temp_V2$Country_V2, sep = "_")
  
  input_table2_to_count_classes_temp_V1_V2 <- full_join(input_table2_to_count_classes_temp_V1,
                                                        input_table2_to_count_classes_temp_V2, by = "UID1")
  
  
  
  #to pheatmap:
  rownames(input_table2_to_count_classes_temp_V1_V2) <- input_table2_to_count_classes_temp_V1_V2$UID1
  
  input_table2_to_count_classes_temp_V1_V2[input_table2_to_count_classes_temp_V1_V2=="Low"] <- -1
  input_table2_to_count_classes_temp_V1_V2[input_table2_to_count_classes_temp_V1_V2=="Normal"] <- 0
  input_table2_to_count_classes_temp_V1_V2[input_table2_to_count_classes_temp_V1_V2=="High"] <- 1
  input_table2_to_count_classes_temp_V1_V2 <- 
    input_table2_to_count_classes_temp_V1_V2[,!colnames(input_table2_to_count_classes_temp_V1_V2) %in%
                                               c("Patient_V1", "Class_V1", "Country_V1",
                                                 "Patient_V2", "Class_V2", "Country_V2", "UID1")]
  
  input_table2_to_count_classes_temp_V1_V2_2 <- sapply( input_table2_to_count_classes_temp_V1_V2, as.numeric)
  rownames(input_table2_to_count_classes_temp_V1_V2_2) <- rownames(input_table2_to_count_classes_temp_V1_V2)
  
  
  input_table2_to_count_classes_temp_V1_V2_2_temp <- data.frame(input_table2_to_count_classes_temp_V1_V2_2)
  input_table2_to_count_classes_temp_V1_V2_2_temp_ordered <- input_table2_to_count_classes_temp_V1_V2_2_temp[, order(names(input_table2_to_count_classes_temp_V1_V2_2_temp))]
  
  #Setting the colorstrips
  tocol_column <- data.frame(rownames(input_table2_to_count_classes_temp_V1_V2))
  colnames(tocol_column) <- "Ids"
  tocol_column$class <- NA
  tocol_column$class[grepl("Ethiopia", tocol_column$Ids)] <- "Ethiopia"
  tocol_column$class[grepl("Kenya", tocol_column$Ids)] <- "Kenya"
  tocol_column$class[grepl("Uganda", tocol_column$Ids)] <- "Uganda"
  tocol_column$class[grepl("Sudan", tocol_column$Ids)] <- "Sudan"
  
  rownames(tocol_column) <- tocol_column$Ids
  tocol_column <- tocol_column[,colnames(tocol_column) !="Ids", drop=FALSE]
  my_colour = list(
    class = c(Ethiopia = "#F8766D", Kenya = "#00BFC4", Uganda ="#C77CFF", Sudan= "#7CAE00")
  )
  
  png(paste(output_prefix_dir, "normal_levels_heatmap.png", sep = "."), width = 2000, height = 2000, res = 300)
  pheatmap(input_table2_to_count_classes_temp_V1_V2_2_temp_ordered, na_col = "grey", color=c("#56B4E9","#009E73","#E69F00"), cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row=2,
           annotation_row =tocol_column, annotation_colors=my_colour, border_color=NA)
  dev.off()
}

# Performs the statistical comparison and plots, both from the Clinical data as well as the legendplex data
# Comparisons between V1 and V2 are done using wilcoxon paired, while between V1 or V2 and HV are done using wilcoxon non-paired
statistical_comparison_all <- function(input_table, output_plots) {
  #https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
  #Statisticall analysis of the differences:
  
  output_prefix_name <- output_plots
  input_table <- melt(input_table, id.vars = c("UID", "Patient", "Class",  "Country"))
  
  input_table$UID <- paste(input_table$Patient, input_table$Country, sep = "_")
  
  #Removing the discrete traits:
  #Note: removed "Spleen_size"
  relevant_data_discrete <- c("Males","Hepatomegaly","Splenomegaly","Auxiliar_Lymphnodes", "Aspirate_grade", "Spleen_size")
  input_table <- input_table[which(!input_table$variable %in% relevant_data_discrete),]
  
  #Selection of the Clinical data:
  
  relevant_data_continuous <- c("Patient", "Class", "Age", "Height", "Weight", "SSystolicBP", "Spleen_size", "Liver_size",
                                "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin",  "WBCells", "Neutrophil", "Aspirate_grade",
                                "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin", "Country")
  
  relevant_data_citokines<- c("CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", "sRAGE", "sST2",
                              "sTNF.RI", "sTNF.RII", "sTREM.1")
  
  analytes_to_test <- unique(as.character(input_table$variable))
  
  dunn_test_table <- data.frame()
  paired_wilcoxon <- data.frame()
  paired_ttest <- data.frame()
  
  samples_to_use_t_plot <- data.frame()
  
  for (c in 1:length(analytes_to_test)) {
    citokyne <- analytes_to_test[c]
    citokine_results <- input_table[input_table$variable==citokyne,]
    
    #Cleaning the data to remove samples missing in V1 or V2
    citokine_results_2 <- citokine_results
    
    #removing the samples that are not present in V1 and V2 for a given citokine:
    citokine_results_2_no_HV <- citokine_results_2[citokine_results_2$Class != "HV",]
    citokine_results_2_no_HV_2_temp <- data.frame(pivot_wider(citokine_results_2_no_HV, values_from = value, names_from = Class))
    
    #Temp testing the medians:
    citokine_results_2_no_HV_2_temp
    citokine_results_2_no_HV_2_temp <- citokine_results_2_no_HV_2_temp[!is.na(citokine_results_2_no_HV_2_temp$V1),]
    citokine_results_2_no_HV_2_temp <- citokine_results_2_no_HV_2_temp[!is.na(citokine_results_2_no_HV_2_temp$V2),]
    ids_to_use <- citokine_results_2_no_HV_2_temp$UID
    ids_to_use2 <- c(ids_to_use, citokine_results_2$UID[citokine_results_2$Class=="HV"])
    
    citokine_results_2 <- citokine_results_2[citokine_results_2$UID %in% ids_to_use2,]
    citokine_results_2[is.na(citokine_results_2$value),]
    
    citokine_results_2 <- citokine_results_2[which(!is.na(citokine_results_2$value)),]
    
    citokine_results_2_no_HV <- citokine_results_2[citokine_results_2$Class != "HV",]
    
    samples_to_use_t_plot <- rbind(samples_to_use_t_plot, citokine_results_2)
    
    
    temp_dunn <- data.frame(dunn_test(samples_to_use_t_plot, value ~ Class, p.adjust.method = "bonferroni") )
    temp_dunn$cytokine <- citokyne
    dunn_test_table <- rbind(dunn_test_table, temp_dunn)
    
    #Getting the max value
    Max_value <- max(citokine_results_2$value, na.rm=TRUE)
    
    #Emptying the DFs:
    temp_df_V1V2paired <- data.frame()
    temp_df_HVV1 <- data.frame()
    temp_df_HVV2 <- data.frame()
    temp_ttest_df_V1V2paired <- data.frame()
    temp_df_ttest_HVV1 <- data.frame()
    temp_df_ttest_HVV2 <- data.frame()
    
    #Comparing V1 and V2 paired 
    citokine_results_3 <- data.frame(citokine_results_2_no_HV %>% pivot_wider(names_from = "Class", values_from = "value"))
    
    #Running the wilcoxon test
    if(length(citokine_results_3$V1) >0 & length(citokine_results_3$V2) > 0) {
      
      temp_wilcoxon_paired <- wilcox.test(citokine_results_3$V1,citokine_results_3$V2, paired = TRUE)
      temp_df_V1V2paired <- data.frame(group1="V1", group2="V2", n1=length(citokine_results_3$V1), n2=length(citokine_results_3$V2), statistic=temp_wilcoxon_paired$statistic, 
                                       n1_median=median(citokine_results_3$V1, na.rm = TRUE), 
                                       n2_median=median(citokine_results_3$V2, na.rm = TRUE), 
                                       n1_mean=mean(citokine_results_3$V1, na.rm = TRUE),
                                       n2_mean=mean(citokine_results_3$V2, na.rm = TRUE),
                                       p=temp_wilcoxon_paired$p.value, cytokine=citokyne, Max.value=Max_value)
      
      
    } else if(length(citokine_results_3$V1)==0 | length(citokine_results_3$V2) ==0) {
      temp_df_V1V2paired <- data.frame(group1="V1", group2="V2", n1=length(HV_values), n2=length(V1_values), 
                                       statistic=NA, 
                                       n1_median=NA, 
                                       n2_median=NA, 
                                       n1_mean=NA, 
                                       n2_mean=NA,  
                                       p=NA, cytokine=citokyne, Max.value=Max_value)
      
      
    }
    
    #Compiring HV and V1 and HV and V2 not paired
    V1_values = citokine_results_2[citokine_results_2$Class=="V1","value"]
    V2_values = citokine_results_2[citokine_results_2$Class=="V2","value"]
    HV_values = citokine_results_2[citokine_results_2$Class=="HV","value"]
    
    
    
    if(length(V1_values)>0 & length(V2_values) >0 & length(HV_values) >0) {
      HV_V1_wilcoxon_not_paired <- wilcox.test(HV_values,V1_values, paired = FALSE)
      temp_df_HVV1 <- data.frame(group1="HV", group2="V1", n1=length(HV_values), n2=length(V1_values), statistic=HV_V1_wilcoxon_not_paired$statistic, 
                                 n1_median=median(HV_values, na.rm = TRUE), 
                                 n2_median=median(V1_values, na.rm = TRUE), 
                                 n1_mean=mean(HV_values, na.rm = TRUE), 
                                 n2_mean=mean(V1_values, na.rm = TRUE),  
                                 p=HV_V1_wilcoxon_not_paired$p.value, cytokine=citokyne, Max.value=Max_value)
      
      HV_V2_wilcoxon_not_paired <- wilcox.test(HV_values,V2_values, paired = FALSE)
      temp_df_HVV2 <- data.frame(group1="HV", group2="V2", n1=length(HV_values), n2=length(V2_values), statistic=HV_V2_wilcoxon_not_paired$statistic, 
                                 n1_median=median(HV_values, na.rm = TRUE), 
                                 n2_median=median(V2_values, na.rm = TRUE),
                                 n1_mean=mean(HV_values, na.rm = TRUE), 
                                 n2_mean=mean(V2_values, na.rm = TRUE),  
                                 p=HV_V2_wilcoxon_not_paired$p.value, cytokine=citokyne, Max.value=Max_value)
      
      
    } else if(length(V1_values)==0 | length(V2_values) ==0  | length(HV_values) == 0) {
      temp_df_HVV1 <- data.frame(group1="HV", group2="V1", n1=length(HV_values), n2=length(V1_values), 
                                 statistic=NA, 
                                 n1_median=NA, 
                                 n2_median=NA, 
                                 n1_mean=NA, 
                                 n2_mean=NA,  
                                 p=NA, cytokine=citokyne, Max.value=Max_value)
      temp_df_HVV2 <- data.frame(group1="HV", group2="V2", n1=length(HV_values), n2=length(V2_values), 
                                 statistic=NA, 
                                 n1_median=NA, 
                                 n2_median=NA,
                                 n1_mean=NA, 
                                 n2_mean=NA,  
                                 p=NA, cytokine=citokyne, Max.value=Max_value)
      
      
      
    }
    
    #Merging the Wilcox data:
    paired_wilcoxon <- rbind(paired_wilcoxon, temp_df_V1V2paired, temp_df_HVV1, temp_df_HVV2)
    
  }
  
  #I will consider separating the multiple testing correction later between cytokines and clinical traits if needed
  
  paired_wilcoxon_temp <- paired_wilcoxon
  paired_wilcoxon$p.adj = p.adjust(paired_wilcoxon$p, method = "BH")
  paired_wilcoxon$p.adj.signif <- "ns"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.05] <- "*"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.01] <- "**"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.001] <- "***"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.0001] <- "****"
  
  
  
  #Exporting the tables
  write.csv(dunn_test_table, paste(output_plots,"dunns.csv", sep = "."), row.names = FALSE, quote = FALSE)
  write.csv(paired_wilcoxon, paste(output_plots,"wilcoxon.csv", sep = "."), row.names = FALSE, quote = FALSE)

  #Generating the plots
    #Formatting the statistical results
  statistical_All4 <- paired_wilcoxon[,c("cytokine", "group1", "group2", "p", "p.adj", "p.adj", "p.adj.signif", "statistic", "Max.value")]
  colnames(statistical_All4) <- c(".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method", "Max.value")
  
  #Cleaning cases there the max value is -inf or inf due to missing data:
  statistical_All4 <- statistical_All4[statistical_All4$Max.value !=-Inf,]
  statistical_All4 <- statistical_All4[statistical_All4$Max.value !=Inf,]
  
  df_size2 <- nrow(statistical_All4)/3 
  statistical_All4$Norm <- rep(c(1.2,1.8,3),df_size2)
  statistical_All4$Max.value2 <- statistical_All4$Max.value*statistical_All4$Norm
  statistical_All4_2 <- statistical_All4[,c("group1", "group2", "p.signif", "Max.value2", ".y.")]
  colnames(statistical_All4_2) <- c("group1", "group2", "p.signif", "y.position", "variable")
  
  #Separating the data for the clinical results or immune markers
  statistical_All4_2_file <- statistical_All4_2[statistical_All4_2$p.signif != "ns",]
  
  #Separating the data for the clinical results or immune markers
  
  statistical_All4_2_file_immune <- statistical_All4_2_file[which(statistical_All4_2_file$variable %in% relevant_data_citokines),] 
  input_table_immune <- samples_to_use_t_plot[which(samples_to_use_t_plot$variable %in% relevant_data_citokines),]
  
  pdf(paste(output_plots, "violin_plots_immunemarkers.pdf", sep = "_"), height = 10, width = 15)
  print(ggplot(input_table_immune, aes(x=Class, y=value)) + geom_violin(alpha=0.7, aes(fill=Class)) + geom_boxplot(width=0.1, fill="white", color="black") +
          theme_minimal() + scale_y_log10() + ylab("pg/mL") + facet_wrap(~variable, scales = "free_y") +
          scale_fill_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")) + 
          theme(strip.text.x = element_text(size = 30)) +
          add_pvalue(statistical_All4_2_file_immune,  tip.length = 0, label.size=5))
  dev.off()
  
  #To plot the clinical data
  df_size2 <- nrow(statistical_All4)/3 
  statistical_All4$Norm <- rep(c(1.1,1.2,1.3),df_size2)
  statistical_All4$Max.value2 <- statistical_All4$Max.value*statistical_All4$Norm
  statistical_All4_2 <- statistical_All4[,c("group1", "group2", "p.signif", "Max.value2", ".y.")]
  colnames(statistical_All4_2) <- c("group1", "group2", "p.signif", "y.position", "variable")
  
  statistical_All4_2_file_clinical <- statistical_All4_2[which(statistical_All4_2$variable %in% relevant_data_continuous),] 
  statistical_All4_2_file_clinical <- statistical_All4_2_file_clinical[statistical_All4_2_file_clinical$p.signif != "ns",]
  
  input_table_clinical <- samples_to_use_t_plot[which(samples_to_use_t_plot$variable %in% relevant_data_continuous),]
  
  input_table_clinical$variable <- factor(input_table_clinical$variable, levels = c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                                                                                    "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                                                                                    "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                                    "Lymphocyte"))
  
  statistical_All4_2_file_clinical$variable <- factor(statistical_All4_2_file_clinical$variable, levels = c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                                                                                                            "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                                                                                                            "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                                                            "Lymphocyte"))
  
  dashed_lines <- data.frame(
    variable = c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin"),
    yintercept = c(13, 4, 1.5, 1.0, 150, 10, 0.7, 3.5, 0.1 ) 
  )
  
  dashed_lines_max <- data.frame(
    variable = c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin"),
    yintercept = c(17, 11, 8, 3.5, 400, 40, 1.3, 5.5, 1.2 ) 
  )
  dashed_lines$variable <- factor(dashed_lines$variable, c( "Albumin", "Hemoglobin",
                                                            "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                            "Lymphocyte"))
  
  dashed_lines_max$variable <- factor(dashed_lines_max$variable, c( "Albumin", "Hemoglobin",
                                                                    "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                    "Lymphocyte"))
  
  pdf(paste(output_plots, "violin_plots_immunemarkers_clinical.pdf", sep = "_"), height = 10, width = 15)
  print(ggplot(input_table_clinical, aes(x=Class, y=value)) + geom_violin(alpha=0.7, aes(fill=Class)) + geom_boxplot(width=0.1, fill="white", color="black") +
          theme_minimal() + ylab("pg/mL") + facet_wrap(~variable, scales = "free_y") +
          scale_fill_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")) + 
          theme(strip.text.x = element_text(size = 20)) +
          add_pvalue(statistical_All4_2_file_clinical,  tip.length = 0, label.size=5) +
          geom_hline(data = dashed_lines, aes(yintercept = yintercept),  linetype = "dashed", color = "#DC143C", size = 0.7) +
          geom_hline(data = dashed_lines_max, aes(yintercept = yintercept),  linetype = "dashed", color = "#228B22", size = 0.7))
  dev.off()
  
  
  #Generaing the tendency plot:
  input_table_clinical_temp <- input_table_clinical[input_table_clinical$Class %in% c("V1","V2"),]
  
  png(paste(output_plots, "Clinical_tendency.png", sep = "_"), height = 2000, width = 4000, res = 300)
  print(ggplot(input_table_clinical_temp, aes(x=Class, y=value, color=UID, group=UID, shape=UID)) + geom_point(size=3) +
          geom_line(size=0.1) +  facet_wrap(~variable, scales = "free_y") + scale_shape_manual(values = rep(c(0,1,2,3,4,5),40)) + theme_bw() +
          geom_hline(data = dashed_lines, aes(yintercept = yintercept),  linetype = "dashed", color = "#DC143C", size = 0.7) +
          geom_hline(data = dashed_lines_max, aes(yintercept = yintercept),  linetype = "dashed", color = "#228B22", size = 0.7))
  dev.off()
  
  
  
  #Now for the radial plots:
  all_plates_melt3_to_radar1_temp_meltmedian <- data.frame(input_table_immune %>%
                                                             group_by(Class, variable) %>% summarise(Median_count = median(value, na.rm = TRUE)))
  
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted <- data.frame(all_plates_melt3_to_radar1_temp_meltmedian %>% pivot_wider(names_from = variable, values_from = Median_count))
  rownames(all_plates_melt3_to_radar1_temp_meltmedian_pivoted) <- all_plates_melt3_to_radar1_temp_meltmedian_pivoted$Class
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted2 <- all_plates_melt3_to_radar1_temp_meltmedian_pivoted[, colnames(all_plates_melt3_to_radar1_temp_meltmedian_pivoted) != "Class"]
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled <- data.frame(apply(all_plates_melt3_to_radar1_temp_meltmedian_pivoted2, 2,  function(x){((x/max(x))*100)}))
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled$visit <- rownames(all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled)
  
  
  #To colour dataframe:
  #I will now only focus on V1 and V2 comparisons for the colour. Colouring blue if it is higher in V2 and orange if V1
  
  paired_wilcoxon_v1v2 <- paired_wilcoxon[which(paired_wilcoxon$group1=="V1" & paired_wilcoxon$group2=="V2"),]
  paired_wilcoxon_v1v2$colour <- NA
  paired_wilcoxon_v1v2$colour[which(paired_wilcoxon_v1v2$n1_mean > paired_wilcoxon_v1v2$n2_mean)] <- "#E69F00"
  paired_wilcoxon_v1v2$colour[which(paired_wilcoxon_v1v2$n2_mean > paired_wilcoxon_v1v2$n1_mean)] <- "#56B4E9"
  paired_wilcoxon_v1v2$colour[which(paired_wilcoxon_v1v2$p.adj.signif=="ns")] <- "black"
  paired_wilcoxon_v1v2 <- paired_wilcoxon_v1v2[which(paired_wilcoxon_v1v2$cytokine %in% relevant_data_citokines),]
  
  
  pdf(paste(output_plots, "Radar_plot_all_patients_combined_scaled_colored.pdf", sep = "."), width = 6, height = 6)
  print(ggRadar(data=all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled,aes(color=visit), legend.position="none", alpha = 0.1, size=0.5, rescale=FALSE) + theme_bw() +
          theme(axis.text=element_text(size=8), theme(legend.title=element_blank()), axis.text.x = element_text(face="bold", color=paired_wilcoxon_v1v2$colour)) + 
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")) +
          scale_fill_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
  dev.off()
  
  
  #Generaing the tendency plot:
  all_data_samples_melt2_temp <- input_table_immune[input_table_immune$Class %in% c("V1","V2"),]
  
  png(paste(output_plots, "tendency.png", sep = "_"), height = 2000, width = 4000, res = 300)
  print(ggplot(all_data_samples_melt2_temp, aes(x=Class, y=value, color=UID, group=UID, shape=UID)) + geom_point(size=3) +
          geom_line(size=0.3) +  facet_wrap(~variable, scales = "free_y") + scale_shape_manual(values = rep(c(0,1,2,3,4,5),40)) + theme_bw())
  dev.off()
  
  png(paste(output_plots, "tendency_log.png", sep = "_"), height = 2000, width = 4000, res = 300)
  print(ggplot(all_data_samples_melt2_temp, aes(x=Class, y=value, color=UID, group=UID, shape=UID)) + geom_point(size=3) +
          geom_line(size=0.3) +  facet_wrap(~variable, scales = "free_y") + scale_shape_manual(values = rep(c(0,1,2,3,4,5),40)) + theme_bw()  + ylab("pg/mL") + scale_y_log10())
  dev.off()
  
  
  
  return(paired_wilcoxon)
  
}
#This function is specific for Female patient data in Kenya:
statistical_comparison_all_noHV_nostats_inPlot <- function(input_table, output_plots) {
  #https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
  #Statisticall analysis of the differences:
  
  output_prefix_name <- output_plots
  input_table <- melt(input_table, id.vars = c("UID", "Patient", "Class",  "Country"))
  
  input_table$UID <- paste(input_table$Patient, input_table$Country, sep = "_")
  
  #Removing the discrete traits:
  #Note: removed "Spleen_size"
  relevant_data_discrete <- c("Males","Hepatomegaly","Splenomegaly","Auxiliar_Lymphnodes", "Aspirate_grade", "Spleen_size")
  input_table <- input_table[which(!input_table$variable %in% relevant_data_discrete),]
  
  #Selection of the Clinical data:
  
  relevant_data_continuous <- c("Patient", "Class", "Age", "Height", "Weight", "SSystolicBP", "Spleen_size", "Liver_size",
                                "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin",  "WBCells", "Neutrophil", "Aspirate_grade",
                                "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin", "Country")
  
  relevant_data_citokines<- c("CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", "sRAGE", "sST2",
                              "sTNF.RI", "sTNF.RII", "sTREM.1")
  
  analytes_to_test <- unique(as.character(input_table$variable))
  
  dunn_test_table <- data.frame()
  paired_wilcoxon <- data.frame()
  paired_ttest <- data.frame()
  
  samples_to_use_t_plot <- data.frame()

  for (c in 1:length(analytes_to_test)) {
    citokyne <- analytes_to_test[c]
    citokine_results <- input_table[input_table$variable==citokyne,]
    
    #Cleaning the data to remove samples missing in V1 or V2
    citokine_results_2 <- citokine_results
    
    #removing the samples that are not present in V1 and V2 for a given citokine:
    citokine_results_2_no_HV <- citokine_results_2[citokine_results_2$Class != "HV",]
    citokine_results_2_no_HV_2_temp <- data.frame(pivot_wider(citokine_results_2_no_HV, values_from = value, names_from = Class))
    
    #Temp testing the medians:
    citokine_results_2_no_HV_2_temp
    citokine_results_2_no_HV_2_temp <- citokine_results_2_no_HV_2_temp[!is.na(citokine_results_2_no_HV_2_temp$V1),]
    citokine_results_2_no_HV_2_temp <- citokine_results_2_no_HV_2_temp[!is.na(citokine_results_2_no_HV_2_temp$V2),]
    ids_to_use <- citokine_results_2_no_HV_2_temp$UID
    ids_to_use2 <- c(ids_to_use, citokine_results_2$UID[citokine_results_2$Class=="HV"])
    
    #This part removes the IDs that are missing in V1 or V2
    citokine_results_2 <- citokine_results_2[citokine_results_2$UID %in% ids_to_use2,]
    citokine_results_2[is.na(citokine_results_2$value),]
    
    citokine_results_2 <- citokine_results_2[which(!is.na(citokine_results_2$value)),]
    
    citokine_results_2_no_HV <- citokine_results_2[citokine_results_2$Class != "HV",]
    
    samples_to_use_t_plot <- rbind(samples_to_use_t_plot, citokine_results_2)
    
    
    temp_dunn <- data.frame(dunn_test(samples_to_use_t_plot, value ~ Class, p.adjust.method = "bonferroni") )
    temp_dunn$cytokine <- citokyne
    dunn_test_table <- rbind(dunn_test_table, temp_dunn)
    
    
    #Getting the max value
    Max_value <- max(citokine_results_2$value, na.rm=TRUE)
    
    #Emptying the DFs:
    temp_df_V1V2paired <- data.frame()
    temp_df_HVV1 <- data.frame()
    temp_df_HVV2 <- data.frame()
    temp_ttest_df_V1V2paired <- data.frame()
    temp_df_ttest_HVV1 <- data.frame()
    temp_df_ttest_HVV2 <- data.frame()
    
    #Comparing V1 and V2 paired 
    citokine_results_3 <- data.frame(citokine_results_2_no_HV %>% pivot_wider(names_from = "Class", values_from = "value"))
    
    # Instead, remove NAs manually
    citokine_results_3 <- citokine_results_3[!is.na(citokine_results_3$V1)]
    citokine_results_3 <- citokine_results_3[!is.na(citokine_results_3$V2)]
    
    #Running the wilcoxon test
    if(length(citokine_results_3$V1) >0 & length(citokine_results_3$V2) > 0) {
      
      temp_wilcoxon_paired <- wilcox.test(citokine_results_3$V1,citokine_results_3$V2, paired = TRUE)
      temp_df_V1V2paired <- data.frame(group1="V1", group2="V2", n1=length(citokine_results_3$V1), n2=length(citokine_results_3$V2), statistic=temp_wilcoxon_paired$statistic, 
                                       n1_median=median(citokine_results_3$V1, na.rm = TRUE), 
                                       n2_median=median(citokine_results_3$V2, na.rm = TRUE), 
                                       n1_mean=mean(citokine_results_3$V1, na.rm = TRUE),
                                       n2_mean=mean(citokine_results_3$V2, na.rm = TRUE),
                                       p=temp_wilcoxon_paired$p.value, cytokine=citokyne, Max.value=Max_value)
      

    } else if(length(citokine_results_3$V1)==0 | length(citokine_results_3$V2) ==0) {
      temp_df_V1V2paired <- data.frame(group1="V1", group2="V2", n1=length(HV_values), n2=length(V1_values), 
                                       statistic=NA, 
                                       n1_median=NA, 
                                       n2_median=NA, 
                                       n1_mean=NA, 
                                       n2_mean=NA,  
                                       p=NA, cytokine=citokyne, Max.value=Max_value)

      
    }
    
    #Compiring HV and V1 and HV and V2 not paired
    V1_values = citokine_results_2[citokine_results_2$Class=="V1","value"]
    V2_values = citokine_results_2[citokine_results_2$Class=="V2","value"]
    HV_values = citokine_results_2[citokine_results_2$Class=="HV","value"]
    
    #Mergint the wilcoxon data:
    # paired_wilcoxon <- rbind(paired_wilcoxon, temp_df_V1V2paired, temp_df_HVV1, temp_df_HVV2)
    paired_wilcoxon <- rbind(paired_wilcoxon, temp_df_V1V2paired)
    
    
  }
  
  #I will consider separating the multiple testing correction later between cytokines and clinical traits if needed
  
  paired_wilcoxon_temp <- paired_wilcoxon
  paired_wilcoxon$p.adj = p.adjust(paired_wilcoxon$p, method = "BH")
  paired_wilcoxon$p.adj.signif <- "ns"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.05] <- "*"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.01] <- "**"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.001] <- "***"
  paired_wilcoxon$p.adj.signif[paired_wilcoxon$p.adj <=0.0001] <- "****"
  
  #Exporting the tables
  write.csv(dunn_test_table, paste(output_plots,"dunns.csv", sep = "."), row.names = FALSE, quote = FALSE)
  write.csv(paired_wilcoxon, paste(output_plots,"wilcoxon.csv", sep = "."), row.names = FALSE, quote = FALSE)

  #Generating the plots
  #Formating the statistical results
  statistical_All4 <- paired_wilcoxon[,c("cytokine", "group1", "group2", "p", "p.adj", "p.adj", "p.adj.signif", "statistic", "Max.value")]
  colnames(statistical_All4) <- c(".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method", "Max.value")
  
  #Cleaning cases there the max value is -inf or inf due to missing data:
  statistical_All4 <- statistical_All4[statistical_All4$Max.value !=-Inf,]
  statistical_All4 <- statistical_All4[statistical_All4$Max.value !=Inf,]
  
  df_size2 <- nrow(statistical_All4)/2
  statistical_All4$Norm <- rep(c(1.2,1.8),df_size2)
  statistical_All4$Max.value2 <- statistical_All4$Max.value*statistical_All4$Norm
  statistical_All4_2 <- statistical_All4[,c("group1", "group2", "p.signif", "Max.value2", ".y.")]
  colnames(statistical_All4_2) <- c("group1", "group2", "p.signif", "y.position", "variable")
  
  #Separating the data for the clinical results or immune markers
  statistical_All4_2_file <- statistical_All4_2[statistical_All4_2$p.signif != "ns",]
 
  
  #Separating the data for the clinical results or immune markers
  statistical_All4_2_file_immune <- statistical_All4_2_file[which(statistical_All4_2_file$variable %in% relevant_data_citokines),] 
  input_table_immune <- samples_to_use_t_plot[which(samples_to_use_t_plot$variable %in% relevant_data_citokines),]
  
  pdf(paste(output_plots, "violin_plots_immunemarkers.pdf", sep = "_"), height = 10, width = 15)
  print(ggplot(input_table_immune, aes(x=Class, y=value)) + geom_violin(alpha=0.7, aes(fill=Class)) + geom_boxplot(width=0.1, fill="white", color="black") +
          theme_minimal() + scale_y_log10() + ylab("pg/mL") + facet_wrap(~variable, scales = "free_y") +
          scale_fill_manual(values = c("#E69F00", "#56B4E9"), breaks = c("V1","V2")) + 
          theme(strip.text.x = element_text(size = 30)))
  dev.off()
  
  #To plot the clinical data
  df_size2 <- nrow(statistical_All4)/2 
  statistical_All4$Norm <- rep(c(1.1,1.2),df_size2)
  statistical_All4$Max.value2 <- statistical_All4$Max.value*statistical_All4$Norm
  statistical_All4_2 <- statistical_All4[,c("group1", "group2", "p.signif", "Max.value2", ".y.")]
  colnames(statistical_All4_2) <- c("group1", "group2", "p.signif", "y.position", "variable")
  
  statistical_All4_2_file_clinical <- statistical_All4_2[which(statistical_All4_2$variable %in% relevant_data_continuous),] 
  statistical_All4_2_file_clinical <- statistical_All4_2_file_clinical[statistical_All4_2_file_clinical$p.signif != "ns",]
  
  input_table_clinical <- samples_to_use_t_plot[which(samples_to_use_t_plot$variable %in% relevant_data_continuous),]
  
  input_table_clinical$variable <- factor(input_table_clinical$variable, levels = c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                                                                                    "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                                                                                    "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                                    "Lymphocyte"))
  
  statistical_All4_2_file_clinical$variable <- factor(statistical_All4_2_file_clinical$variable, levels = c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                                                                                                            "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                                                                                                            "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                                                            "Lymphocyte"))
  
  dashed_lines <- data.frame(
    variable = c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin"),
    yintercept = c(13, 4, 1.5, 1.0, 150, 10, 0.7, 3.5, 0.1 )  
  )
  
  dashed_lines_max <- data.frame(
    variable = c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin"),
    yintercept = c(17, 11, 8, 3.5, 400, 40, 1.3, 5.5, 1.2 ) 
  )
  dashed_lines$variable <- factor(dashed_lines$variable, c( "Albumin", "Hemoglobin",
                                                            "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                            "Lymphocyte"))
  
  dashed_lines_max$variable <- factor(dashed_lines_max$variable, c( "Albumin", "Hemoglobin",
                                                                    "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                    "Lymphocyte"))
  
  pdf(paste(output_plots, "violin_plots_immunemarkers_clinical.pdf", sep = "_"), height = 10, width = 15)
  print(ggplot(input_table_clinical, aes(x=Class, y=value)) + geom_violin(alpha=0.7, aes(fill=Class)) + geom_boxplot(width=0.1, fill="white", color="black") +
          theme_minimal() + ylab("pg/mL") + facet_wrap(~variable, scales = "free_y") +
          scale_fill_manual(values = c("#E69F00", "#56B4E9"), breaks = c("V1","V2")) + 
          theme(strip.text.x = element_text(size = 20)) +
          geom_hline(data = dashed_lines, aes(yintercept = yintercept),  linetype = "dashed", color = "#DC143C", size = 0.7) +
          geom_hline(data = dashed_lines_max, aes(yintercept = yintercept),  linetype = "dashed", color = "#228B22", size = 0.7))
  dev.off()
  
  
  #Generaing the tendency plot:
  input_table_clinical_temp <- input_table_clinical[input_table_clinical$Class %in% c("V1","V2"),]
  
  png(paste(output_plots, "Clinical_tendency.png", sep = "_"), height = 2000, width = 4000, res = 300)
  print(ggplot(input_table_clinical_temp, aes(x=Class, y=value, color=UID, group=UID, shape=UID)) + geom_point(size=3) +
          geom_line(size=0.1) +  facet_wrap(~variable, scales = "free_y") + scale_shape_manual(values = rep(c(0,1,2,3,4,5),40)) + theme_bw() +
          geom_hline(data = dashed_lines, aes(yintercept = yintercept),  linetype = "dashed", color = "#DC143C", size = 0.7) +
          geom_hline(data = dashed_lines_max, aes(yintercept = yintercept),  linetype = "dashed", color = "#228B22", size = 0.7))
  dev.off()
  
  
  
  #Now for the radial plots:
  all_plates_melt3_to_radar1_temp_meltmedian <- data.frame(input_table_immune %>%
                                                             group_by(Class, variable) %>% summarise(Median_count = median(value, na.rm = TRUE)))
  
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted <- data.frame(all_plates_melt3_to_radar1_temp_meltmedian %>% pivot_wider(names_from = variable, values_from = Median_count))
  rownames(all_plates_melt3_to_radar1_temp_meltmedian_pivoted) <- all_plates_melt3_to_radar1_temp_meltmedian_pivoted$Class
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted2 <- all_plates_melt3_to_radar1_temp_meltmedian_pivoted[, colnames(all_plates_melt3_to_radar1_temp_meltmedian_pivoted) != "Class"]
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled <- data.frame(apply(all_plates_melt3_to_radar1_temp_meltmedian_pivoted2, 2,  function(x){((x/max(x))*100)}))
  all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled$visit <- rownames(all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled)
  
  
  #To colour dataframe:
  #I will now only focus on V1 and V2 comparisons for the colour. Colouring blue if it is higher in V2 and orange if V1
  
  paired_wilcoxon_v1v2 <- paired_wilcoxon[which(paired_wilcoxon$group1=="V1" & paired_wilcoxon$group2=="V2"),]
  paired_wilcoxon_v1v2$colour <- NA
  paired_wilcoxon_v1v2$colour[which(paired_wilcoxon_v1v2$n1_median > paired_wilcoxon_v1v2$n2_median)] <- "#E69F00"
  paired_wilcoxon_v1v2$colour[which(paired_wilcoxon_v1v2$n2_median > paired_wilcoxon_v1v2$n1_median)] <- "#56B4E9"
  paired_wilcoxon_v1v2$colour[which(paired_wilcoxon_v1v2$p.adj.signif=="ns")] <- "black"
  paired_wilcoxon_v1v2 <- paired_wilcoxon_v1v2[which(paired_wilcoxon_v1v2$cytokine %in% relevant_data_citokines),]
  
  
  pdf(paste(output_plots, "Radar_plot_all_patients_combined_scaled_colored.pdf", sep = "."), width = 6, height = 6)
  print(ggRadar(data=all_plates_melt3_to_radar1_temp_meltmedian_pivoted2_scaled,aes(color=visit), legend.position="none", alpha = 0.1, size=0.5, rescale=FALSE) + theme_bw() +
          theme(axis.text=element_text(size=8), theme(legend.title=element_blank()), axis.text.x = element_text(face="bold", color=paired_wilcoxon_v1v2$colour)) + 
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2", "HV")) +
          scale_fill_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2", "HV")))
  dev.off()
  
  #Generaing the tendency plot:
  all_data_samples_melt2_temp <- input_table_immune[input_table_immune$Class %in% c("V1","V2"),]
  
  png(paste(output_plots, "tendency.png", sep = "_"), height = 2000, width = 4000, res = 300)
  print(ggplot(all_data_samples_melt2_temp, aes(x=Class, y=value, color=UID, group=UID, shape=UID)) + geom_point(size=3) +
          geom_line(size=0.3) +  facet_wrap(~variable, scales = "free_y") + scale_shape_manual(values = rep(c(0,1,2,3,4,5),40)) + theme_bw())
  dev.off()
  
  png(paste(output_plots, "tendency_log.png", sep = "_"), height = 2000, width = 4000, res = 300)
  print(ggplot(all_data_samples_melt2_temp, aes(x=Class, y=value, color=UID, group=UID, shape=UID)) + geom_point(size=3) +
          geom_line(size=0.3) +  facet_wrap(~variable, scales = "free_y") + scale_shape_manual(values = rep(c(0,1,2,3,4,5),40)) + theme_bw()  + ylab("pg/mL") + scale_y_log10())
  dev.off()
  
  
  
  return(paired_wilcoxon)
  
}

#Statistical comparison, between two given groups, using wilcoxon non-paired:
statsitical_compare_plot_group <- function(input_table, feature_to_evaluate, outname) {
  
  select_info_to_plot <- c( "UID", feature_to_evaluate,
                            "Age", "Height", "Weight",  "SSystolicBP", "SDiastolicBP", "Temperature",  "Pulse", "Hemoglobin",
                            "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                            "Aspirate_grade", "Spleen_size", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", 
                            "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")
  
  input_table2 <- input_table[,colnames(input_table) %in% select_info_to_plot]
  
  analytes_to_test <- unique(as.character(colnames(input_table2)))
  analytes_to_test <- analytes_to_test[!analytes_to_test %in% c("UID", feature_to_evaluate)]
  
  input_table_melt <- melt(input_table2, id.vars = c("UID", feature_to_evaluate))
  
  #Run the comparison for each:
  wilcoxon_results <- data.frame()
  ks_test_results <- data.frame()
  
  for (c in 1:length(analytes_to_test)) {
    citokyne <- analytes_to_test[c]
    citokine_results <- input_table_melt[input_table_melt$variable==citokyne,]
    
    #Removing the lines with NA for the analyse:
    citokine_results <- citokine_results[!is.na(citokine_results$value),]
    
    V1_values = citokine_results[ citokine_results[[feature_to_evaluate]] == 0,"value"]
    V2_values = citokine_results[ citokine_results[[feature_to_evaluate]] == 1,"value"]
    
    Max_value <- max(citokine_results$value, na.rm = TRUE)
    
    if(length(V1_values)>=2 & length(V2_values) >=2 ) {
      V1_V2_wilcoxon_not_paired <- wilcox.test(V1_values,V2_values, paired = FALSE)
      temp_df_V1V2 <- data.frame(group1="Class_0", group2="Class_1", n1=length(V1_values), n2=length(V2_values), statistic=V1_V2_wilcoxon_not_paired$statistic, 
                                 n1_median=median(V1_values, na.rm = TRUE), 
                                 n2_median=median(V2_values, na.rm = TRUE), 
                                 n1_mean=mean(V1_values, na.rm = TRUE), 
                                 n2_mean=mean(V2_values, na.rm = TRUE),  
                                 p=V1_V2_wilcoxon_not_paired$p.value, cytokine=citokyne, Max.value=Max_value)
      
      wilcoxon_results <- rbind(wilcoxon_results, temp_df_V1V2)
      
      
      V1_V2_ks_test <- ks.test(V1_values,V2_values)
      ks_temp_V1V2 <- data.frame(group1="Class_0", group2="Class_1", n1=length(V1_values), n2=length(V2_values), statistic=V1_V2_ks_test$statistic, 
                                 n1_median=median(V1_values, na.rm = TRUE), 
                                 n2_median=median(V2_values, na.rm = TRUE), 
                                 n1_mean=mean(V1_values, na.rm = TRUE), 
                                 n2_mean=mean(V2_values, na.rm = TRUE),  
                                 p=V1_V2_ks_test$p.value, cytokine=citokyne, Max.value=Max_value)
      
      ks_test_results <- rbind(ks_test_results, ks_temp_V1V2)
      
    }

  }
  
  
  #Note, I will use the p-value instead of the adjusted -pvalue:
  wilcoxon_results <- wilcoxon_results[!is.na(wilcoxon_results$n1_median),]
  wilcoxon_results <- wilcoxon_results[!is.na(wilcoxon_results$n2_median),]
  wilcoxon_results$p.adj = p.adjust(wilcoxon_results$p, method = "BH")
  wilcoxon_results$p.adj.signif <- "ns"
  wilcoxon_results$p.adj.signif[wilcoxon_results$p.adj <=0.05] <- "*"
  wilcoxon_results$p.adj.signif[wilcoxon_results$p.adj <=0.01] <- "**"
  wilcoxon_results$p.adj.signif[wilcoxon_results$p.adj <=0.001] <- "***"
  wilcoxon_results$p.adj.signif[wilcoxon_results$p.adj <=0.0001] <- "****"
  
  wilcoxon_results <- wilcoxon_results[!is.na(wilcoxon_results$n1_median),]
  wilcoxon_results <- wilcoxon_results[!is.na(wilcoxon_results$n2_median),]
  wilcoxon_results$p.adj = p.adjust(wilcoxon_results$p, method = "BH")
  wilcoxon_results$p.signif <- "ns"
  wilcoxon_results$p.signif[wilcoxon_results$p <=0.05] <- "*"
  wilcoxon_results$p.signif[wilcoxon_results$p <=0.01] <- "**"
  wilcoxon_results$p.signif[wilcoxon_results$p <=0.001] <- "***"
  wilcoxon_results$p.signif[wilcoxon_results$p <=0.0001] <- "****"
  
  write.csv(wilcoxon_results, paste(outname, feature_to_evaluate, "wilcoxon.csv", sep = "."), row.names = FALSE, quote = FALSE)
  write.csv(ks_test_results, paste(outname, feature_to_evaluate, "KStest.csv", sep = "."), row.names = FALSE, quote = FALSE)
  
  
  #Formating the statistical results
  statistical_All4 <- wilcoxon_results[,c("cytokine", "group1", "group2", "p", "p.adj", "p.adj", "p.signif", "statistic", "Max.value")]
  colnames(statistical_All4) <- c(".y.", "group1", "group2", "p", "p.adj", "p.format", "p.signif", "method", "Max.value")
  
  statistical_All4$Norm <- 1.2
  statistical_All4$Max.value2 <- statistical_All4$Max.value*statistical_All4$Norm
  statistical_All4_2 <- statistical_All4[,c("group1", "group2", "p.signif", "Max.value2", ".y.")]
  colnames(statistical_All4_2) <- c("group1", "group2", "p.signif", "y.position", "variable")
  
  statistical_All4_2$group1 <- gsub("Class_0","0", statistical_All4_2$group1)
  statistical_All4_2$group2 <- gsub("Class_1","1", statistical_All4_2$group2)
  
  #Separating the data for the clinical results or immune markers
  colnames(input_table_melt)[2] <- "Class"
  input_table_melt[,"Class"] <- as.character(input_table_melt[,"Class"])
  
  #Selecting the order and using factors:
  
  levels_to_select <-  c("Age", "Height", "Weight", "Temperature", "SSystolicBP",
                         "SDiastolicBP", "Pulse", "Aspirate_grade", "Spleen_size", "Albumin", "Hemoglobin",
                         "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                         "Lymphocyte", "CX3CL1",  "CXCL12", "PTX3", "TGF.B1", "sCD40L", "sRAGE", "sCD25",
                         "sST2", "sTNF.RI", "sTNF.RII", "sTREM.1" )
  
  input_table_melt$variable <- factor(input_table_melt$variable, levels =levels_to_select)
  
  statistical_All4_2$variable <- factor(statistical_All4_2$variable, levels = levels_to_select)
  statistical_All4_2 <- statistical_All4_2[order(statistical_All4_2$variable), ]
  
  input_table_melt <- input_table_melt[order(input_table_melt$variable),]
  
  
  pdf(paste(outname, feature_to_evaluate, "violin_plots.pdf", sep = "_"), height = 18, width = 22)
  print(ggplot(input_table_melt, aes(x=Class, y=value)) + geom_violin(alpha=0.7, aes(fill=Class)) + 
          geom_boxplot(width=0.1, fill="white", color="black") +
          theme_minimal() + scale_y_log10() + ylab("Value") + facet_wrap(~variable, scales = "free_y") +
          scale_fill_manual(values = c("#DC143C", "#87CEEB"), breaks = c("1","0")) + 
          theme(strip.text.x = element_text(size = 30)) + add_pvalue(statistical_All4_2,  tip.length = 0, label.size=5))
  dev.off()
  
  #Plotting the distributions:
  pdf(paste(outname, feature_to_evaluate, "ks_plots.pdf", sep = "_"), height = 18, width = 22)
  print(ggplot(input_table_melt, aes(fill=Class, x=value)) + geom_density(alpha=0.5) +  
          theme_minimal() + ylab("Value") + facet_wrap(~variable, scales = "free") +
          scale_fill_manual(values = c("#DC143C", "#87CEEB"), breaks = c("1","0")) + 
          theme(strip.text.x = element_text(size = 30)))
  dev.off()
  
  return(statistical_All4_2)
  
}

#Generates the PCA and UMAP analysis and plots for the legendplex data:
pca_and_umap_funcion <- function(data_table, output_plots, umap_number_of_neighbors) {
  #PCA and UMAP
  #all_data_samples_melt2_temp
  #https://cran.r-project.org/web/packages/umap/vignettes/umap.html
  #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
  #https://www.datacamp.com/tutorial/pca-analysis-r
  #https://plotly.com/r/t-sne-and-umap-projections/
  
  #First for the PCA:
  
  #subseting the DF based on columns of interest:
  relevant_data_citokines<- c("CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", "sRAGE", "sST2",
                              "sTNF.RI", "sTNF.RII", "sTREM.1")
  
  input_table_umap1 <- data_table[, colnames(data_table) %in% c("UID", "Patient", "Class", "Country", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", "sRAGE", "sST2",
                                                                "sTNF.RI", "sTNF.RII", "sTREM.1")]
  #Removing all NAs
  input_table_umap2 <- input_table_umap1[complete.cases(input_table_umap1), ]
  
  #Removing columns if all values are equal for the numeric columns
  
  input_table_umap2 <- input_table_umap2[, sapply(input_table_umap2, function(col) {
    if (is.numeric(col)) {
      length(unique(col)) > 1
    } else {
      TRUE  # keep non-numeric columns as is
    }
  })]
  
  
  
  #lognorm before PCA - UMAP
  num_cols <- sapply(input_table_umap2, is.numeric)
  input_table_umap2[num_cols]<- lapply(input_table_umap2[num_cols], function(x) log(x + 0.001))
  
  #Scalling the dataset set for PCA and UMAP####
  input_table_umap3 <- input_table_umap2[,5:ncol(input_table_umap2)]
  input_table_umap3_scaled <- data.frame(scale(input_table_umap3, center=TRUE, scale=TRUE))
  input_table_umap3_scaled$group <- input_table_umap2$group
  input_table_umap3_scaled$UID <- input_table_umap2$UID
  rownames(input_table_umap3_scaled) <- input_table_umap3_scaled$UID
  
  input_table_umap3_scaled$Country <- input_table_umap2$Country
  input_table_umap3_scaled$Class <- input_table_umap2$Class
  
  
  input_table2 <- input_table_umap3_scaled[,!colnames(input_table_umap3_scaled) %in% c("UID", "Country", "Class"),]
  res.pca <- prcomp(input_table2, scale = FALSE)
  
  pdf(paste(output_plots,"PCAs_relevance.pdf", sep = "."), height = 5, width = 5)
  print(fviz_eig(res.pca))
  dev.off()
  
  pdf(paste(output_plots,"PC_directions.pdf", sep = "."), height = 7, width = 7)
  print(fviz_pca_var(res.pca,
                     col.var = "contrib", # Color by contributions to the PC
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                     repel = TRUE     # Avoid text overlapping
  ))
  dev.off()
  
  pdf(paste(output_plots,"elevance_each_marker.pdf", sep = "."), height = 5, width = 5)
  print(fviz_cos2(res.pca, choice = "var", axes = 1:2))
  dev.off()
  
  pdf(paste(output_plots,"PCA2_dim_biplot_legendplex_Inf.pdf", sep = "."), height = 8, width = 8)
  print(autoplot(res.pca, size=5, data = input_table_umap3_scaled, colour = 'Class', shape="Country", loadings = TRUE, loadings.label = TRUE, alpha=0.7,
                 loadings.label.size = 3, loadings.colour = 'blue') + 
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")) +
          theme_bw())
  dev.off()

  
  #Individuals coordinates:
  ind.coord <- data.frame(res.pca$x)
  ind.coord$Class <- input_table_umap3_scaled$Class
  ind.coord$UID <- input_table_umap3_scaled$UID
  
  #Summary of the variance:
  summary_PCA <- summary(res.pca)
  summary_PCA <- summary_PCA$importance
  pca1_relevance <- round((summary_PCA[2,1])*100,2)
  pca2_relevance <- round((summary_PCA[2,2])*100,2)
  
  pdf(paste(output_plots,"PCA2_dim.pdf", sep = "."), height = 10, width = 10)
  print(ggplot(ind.coord, aes(x=PC1, y=PC2, color=Class, label = UID)) + geom_point() + theme_bw() + geom_text_repel(size=3, show.legend = FALSE, max.overlaps = Inf) +
          geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
          geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + xlab(paste("PCA1(", pca1_relevance,"%)", sep = "")) +
          ylab(paste("PCA2(", pca2_relevance,"%)", sep = "")) + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV"))) 
  dev.off()
  
  
  if(length(data_table$Country) > 0 ) {
    ind.coord2 <- ind.coord
    ind.coord2$Coutry <- input_table_umap3_scaled$Country
    
    pdf(paste(output_plots,"PCA2_dim.pdf", sep = "."), height = 10, width = 10)
    print(ggplot(ind.coord2, aes(x=PC1, y=PC2, color=Class, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + xlab(paste("PCA1(", pca1_relevance,"%)", sep = "")) +
            ylab(paste("PCA2(", pca2_relevance,"%)", sep = "")) + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV"))) 
    dev.off()
    
  }
  
  
  #PCA inverted, with the arrows pointing to the mean of the samples:
  #Getting the mean of each citokine for each visit
  
  input_table_melt <- melt(input_table_umap3_scaled, id.vars = c("UID", "Country", "Class"))
  
  input_table_mel_median <- data.frame(input_table_melt %>% group_by(Class, variable) %>% summarise(Median_value = median(value)))
  input_table_mel_median_pivoted <- data.frame(pivot_wider(input_table_mel_median, names_from = variable, values_from = Median_value))
  rownames(input_table_mel_median_pivoted) <- input_table_mel_median_pivoted$Class
  input_table_mel_median_pivoted <- input_table_mel_median_pivoted[,colnames(input_table_mel_median_pivoted) !="Class"]
  input_table_mel_median_pivoted_t <- t(input_table_mel_median_pivoted)
  
  res.pca_inv <- prcomp(input_table_mel_median_pivoted_t, scale = FALSE)
  
  my.col.var <- c("#56B4E9", "#CC79A7", "#E69F00")
  
  pdf(paste(output_plots,"PC_directions_inverted_citokines.pdf", sep = "."), height = 7, width = 7)
  print(fviz_pca_biplot(res.pca_inv,
                        col.var = c("V1", "V2", "HV"), # Color by contributions to the PC
                        palette = my.col.var, 
                        arrowsize = 2, labelsize = 3,
  )+theme(legend.position = "none") + theme(axis.text = element_text(size = 14)) 
  )
  dev.off()
  
  pdf(paste(output_plots,"PC_directions_inverted_citokines_1_3.pdf", sep = "."), height = 7, width = 7)
  print(fviz_pca_biplot(res.pca_inv, axes = c(1, 3),
                        col.var = c("V1", "V2", "HV"), # Color by contributions to the PC
                        palette = my.col.var, 
                        arrowsize = 2, labelsize = 3,
  )+theme(legend.position = "none") + theme(axis.text = element_text(size = 14)) 
  )
  dev.off()
  
  #UMAP
  umap_test1 <- umap(input_table2, n_epochs=1000, n_neighbors=umap_number_of_neighbors)
  umap_test2 <- data.frame(umap_test1$layout)
  
  colnames(umap_test2) <- c("dim1", "dim2")
  umap_test2$Class <- input_table_umap2$Class
  umap_test2$UID <- rownames(umap_test2)
  
  
  pdf(paste(output_plots,"UMAP_2dim.pdf", sep = "."), height = 9, width = 9)
  print(ggplot(umap_test2, aes(x=dim1, y=dim2, color=Class, label = UID)) + geom_point() + theme_bw() + geom_text_repel(size=3, show.legend = FALSE, max.overlaps = Inf) +
          geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
          geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
  dev.off()
  
  
  if(length(data_table$Country) > 0 ) {
    umap_test2_2<- umap_test2
    umap_test2_2$Coutry <- umap_test2_2$UID
    umap_test2_2$Coutry <- gsub(".*_","",umap_test2_2$Coutry)
    
    pdf(paste(output_plots,"UMAP_2dim_contries_immune.pdf", sep = "."), height = 9, width = 9)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
    dev.off()
    
    pdf(paste(output_plots,"UMAP_2dim_coutries_split_immune.pdf", sep = "."), height = 9, width = 9)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) + facet_wrap(~Coutry) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
    dev.off()
    
    #Umap_short:
    pdf(paste(output_plots,"UMAP_2dim_contries_short_immune.pdf", sep = "."), height = 5, width = 5)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
    dev.off()
    
    
    hepatomegaly_Ids <- data_table$UID[which(data_table$Hepatomegaly == 1)]
    splenomegaly_Ids <- data_table$UID[which(data_table$Splenomegaly == 1)]
    hepatosplenomegaly_Ids <- data_table$UID[which(data_table$Hepatomegaly == 1 & data_table$Splenomegaly == 1) ]
    
    hepatomegaly_Ids <- setdiff(hepatomegaly_Ids, hepatosplenomegaly_Ids)
    splenomegaly_Ids <- setdiff(splenomegaly_Ids, hepatosplenomegaly_Ids)
    
    umap_test2_2$Clinical <- "No"
    umap_test2_2$Clinical[umap_test2_2$UID %in% splenomegaly_Ids] <- "Splenomegaly"
    umap_test2_2$Clinical[umap_test2_2$UID %in% hepatomegaly_Ids] <- "Hepatomegaly"
    umap_test2_2$Clinical[umap_test2_2$UID %in% hepatosplenomegaly_Ids] <- "Hepatosplenomegaly"
    
    pdf(paste(output_plots,"UMAP_2dim_contries_short_immune_hepatospleno.pdf", sep = "."), height = 5, width = 5)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                     breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")))
    dev.off()
    
    pdf(paste(output_plots,"UMAP_2dim_contries_short_immune_hepatospleno.pdf", sep = "."), height = 3, width = 10)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                     breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")) +
            facet_wrap(~Class))
    dev.off()
    
    pdf(paste(output_plots,"UMAP_2dim_contries_short_immune_hepatospleno2.pdf", sep = "."), height = 5, width = 4.5)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=3, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                     breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")) +
            facet_wrap(~Class, ncol = 1))
    dev.off()
    
  }
  
  
  #Pheatmap
  to_pca_umap3_scaled_t <- data.frame(t(input_table2))
  tocol_column <- input_table_umap2[,c("Class"), drop =FALSE]
  rownames(tocol_column) <- input_table_umap2$UID
  my_colour = list(
    Class = c(HV = "#CC79A7", V1 = "#E69F00", V2 ="#56B4E9")
  )
  
  input_table_umap2_to_heatmap <- input_table_umap2[, !colnames(input_table_umap2) %in% c("UID", "Patient", "Class", "Country")]
  rownames(input_table_umap2_to_heatmap) <- input_table_umap2$UID
  input_table_umap2_to_heatmap2 <- t(input_table_umap2_to_heatmap)
  
  pdf(paste(output_plots,"heatmap_scaled.pdf", sep = "."), height = 9, width = 18)
  print(pheatmap(input_table_umap2_to_heatmap2, clustering_distance_cols = "manhattan", clustering_method = "average",
                 scale="row", annotation_col=tocol_column, annotation_colors=my_colour
  ))
  dev.off()
  
  if(length(data_table$Country) > 0 ) {
    pdf(paste(output_plots,"heatmap_scaled_large.pdf", sep = "."), height = 9, width = 22)
    print(pheatmap(input_table_umap2_to_heatmap2, clustering_distance_cols = "manhattan", clustering_method = "average",
                   fontsize_col=3, scale="row", annotation_col=tocol_column, annotation_colors=my_colour))
    dev.off()
    
  }
  
}

#Generates the PCA analysis and plots for the clinical data:
pca_and_umap_funcion_clnical <- function(data_table, output_plots, umap_number_of_neighbors) {
  #PCA and UMAP
  #all_data_samples_melt2_temp
  #https://cran.r-project.org/web/packages/umap/vignettes/umap.html
  #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
  #https://www.datacamp.com/tutorial/pca-analysis-r
  #https://plotly.com/r/t-sne-and-umap-projections/
  
  #First for the PCA:
  #subseting the DF based on columns of interest:
  relevant_data_citokines<- c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin")
  
  input_table_umap1 <- data_table[c("UID", "Patient", "Class", "Country", "Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin")]
  #Removing all NAs
  input_table_umap2 <- input_table_umap1[complete.cases(input_table_umap1), ]
  
  #Scalling the dataset set for PCA and UMAP####
  input_table_umap3 <- input_table_umap2[,5:ncol(input_table_umap2)]
  input_table_umap3_scaled <- data.frame(scale(input_table_umap3, center=TRUE, scale=TRUE))
  input_table_umap3_scaled$group <- input_table_umap2$group
  input_table_umap3_scaled$UID <- input_table_umap2$UID
  rownames(input_table_umap3_scaled) <- input_table_umap3_scaled$UID
  
  input_table_umap3_scaled$Country <- input_table_umap2$Country
  input_table_umap3_scaled$Class <- input_table_umap2$Class
  
  
  input_table2 <- input_table_umap3_scaled[,!colnames(input_table_umap3_scaled) %in% c("UID", "Country", "Class"),]
  res.pca <- prcomp(input_table2, scale = FALSE)
  
  pdf(paste(output_plots,"PCAs_clinical_relevance.pdf", sep = "."), height = 5, width = 5)
  print(fviz_eig(res.pca))
  dev.off()
  
  pdf(paste(output_plots,"PC_directions_clinical.pdf", sep = "."), height = 7, width = 7)
  print(fviz_pca_var(res.pca,
                     col.var = "contrib", # Color by contributions to the PC
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                     repel = TRUE     # Avoid text overlapping
  ))
  dev.off()
  
  pdf(paste(output_plots,"elevance_each_marker_clinical.pdf", sep = "."), height = 5, width = 5)
  print(fviz_cos2(res.pca, choice = "var", axes = 1:2))
  dev.off()
  
  pdf(paste(output_plots,"PCA2_dim_biplot_clinical.pdf", sep = "."), height = 8, width = 8)
  print(autoplot(res.pca, size=5, data = input_table_umap3_scaled, colour = 'Class', shape="Country", loadings = TRUE, loadings.label = TRUE, alpha=0.7,
                 loadings.label.size = 3, loadings.colour = 'blue') + 
          scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")) +
          theme_bw())
  dev.off()
  #Individuals coordinates:
  ind.coord <- data.frame(res.pca$x)
  ind.coord$Class <- input_table_umap3_scaled$Class
  ind.coord$UID <- input_table_umap3_scaled$UID
  
  #Summary of the variance:
  summary_PCA <- summary(res.pca)
  summary_PCA <- summary_PCA$importance
  pca1_relevance <- round((summary_PCA[2,1])*100,2)
  pca2_relevance <- round((summary_PCA[2,2])*100,2)
  
  pdf(paste(output_plots,"PCA2_dim_clinical.pdf", sep = "."), height = 10, width = 10)
  print(ggplot(ind.coord, aes(x=PC1, y=PC2, color=Class, label = UID)) + geom_point() + theme_bw() + geom_text_repel(size=3, show.legend = FALSE, max.overlaps = Inf) +
          geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
          geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + xlab(paste("PCA1(", pca1_relevance,"%)", sep = "")) +
          ylab(paste("PCA2(", pca2_relevance,"%)", sep = "")) + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV"))) 
  dev.off()
  
  
  if(length(data_table$Country) > 0 ) {
    ind.coord2 <- ind.coord
    ind.coord2$Coutry <- input_table_umap3_scaled$Country
    
    pdf(paste(output_plots,"PCA2_dim_clinical.pdf", sep = "."), height = 10, width = 10)
    print(ggplot(ind.coord2, aes(x=PC1, y=PC2, color=Class, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + xlab(paste("PCA1(", pca1_relevance,"%)", sep = "")) +
            ylab(paste("PCA2(", pca2_relevance,"%)", sep = "")) + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV"))) 
    dev.off()
    
  }
  
  
  #PCA inverted, with the arrows pointing to the mean of the samples:
  #Getting the mean of each citokine for each visit
  
  input_table_melt <- melt(input_table_umap3_scaled, id.vars = c("UID", "Country", "Class"))
  
  input_table_mel_median <- data.frame(input_table_melt %>% group_by(Class, variable) %>% summarise(Median_value = median(value)))
  input_table_mel_median_pivoted <- data.frame(pivot_wider(input_table_mel_median, names_from = variable, values_from = Median_value))
  rownames(input_table_mel_median_pivoted) <- input_table_mel_median_pivoted$Class
  input_table_mel_median_pivoted <- input_table_mel_median_pivoted[,colnames(input_table_mel_median_pivoted) !="Class"]
  input_table_mel_median_pivoted_t <- t(input_table_mel_median_pivoted)
  
  res.pca_inv <- prcomp(input_table_mel_median_pivoted_t, scale = FALSE)
  
  my.col.var <- c("#56B4E9", "#CC79A7", "#E69F00")
  
  pdf(paste(output_plots,"PC_directions_inverted_clinical.pdf", sep = "."), height = 7, width = 7)
  print(fviz_pca_biplot(res.pca_inv,
                        col.var = c("V1", "V2", "HV"), # Color by contributions to the PC
                        palette = my.col.var, 
                        arrowsize = 2, labelsize = 3,
  )+theme(legend.position = "none") + theme(axis.text = element_text(size = 14)) 
  )
  dev.off()
  
  pdf(paste(output_plots,"PC_directions_inverted_clinical_1_3.pdf", sep = "."), height = 7, width = 7)
  print(fviz_pca_biplot(res.pca_inv, axes = c(1, 3),
                        col.var = c("V1", "V2", "HV"), # Color by contributions to the PC
                        palette = my.col.var, 
                        arrowsize = 2, labelsize = 3,
  )+theme(legend.position = "none") + theme(axis.text = element_text(size = 14)) 
  )
  dev.off()
  
  #UMAP
  umap_test1 <- umap(input_table2, n_epochs=1000, n_neighbors=umap_number_of_neighbors)
  umap_test2 <- data.frame(umap_test1$layout)
  
  colnames(umap_test2) <- c("dim1", "dim2")
  umap_test2$Class <- input_table_umap2$Class
  umap_test2$UID <- rownames(umap_test2)
  
  
  pdf(paste(output_plots,"UMAP_2dim_clinical.pdf", sep = "."), height = 9, width = 9)
  print(ggplot(umap_test2, aes(x=dim1, y=dim2, color=Class, label = UID)) + geom_point() + theme_bw() + geom_text_repel(size=3, show.legend = FALSE, max.overlaps = Inf) +
          geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
          geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
  dev.off()
  
  
  if(length(data_table$Country) > 0 ) {
    umap_test2_2<- umap_test2
    umap_test2_2$Coutry <- umap_test2_2$UID
    umap_test2_2$Coutry <- gsub(".*_","",umap_test2_2$Coutry)
    
    pdf(paste(output_plots,"UMAP_2dim_contries_clinical.pdf", sep = "."), height = 9, width = 9)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
    dev.off()
    
    pdf(paste(output_plots,"UMAP_2dim_coutries_split_clinical.pdf", sep = "."), height = 9, width = 9)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) + facet_wrap(~Coutry) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
    dev.off()
    
    #Umap_short:
    pdf(paste(output_plots,"UMAP_2dim_contries_short_clinical.pdf", sep = "."), height = 5, width = 5)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
    dev.off()
    
    #Umap colouring by hepatomegaly and splenomegaly
    hepatomegaly_Ids <- data_table$UID[which(data_table$Hepatomegaly == 1)]
    splenomegaly_Ids <- data_table$UID[which(data_table$Splenomegaly == 1)]
    hepatosplenomegaly_Ids <- data_table$UID[which(data_table$Hepatomegaly == 1 & data_table$Splenomegaly == 1) ]
    
    hepatomegaly_Ids <- setdiff(hepatomegaly_Ids, hepatosplenomegaly_Ids)
    splenomegaly_Ids <- setdiff(splenomegaly_Ids, hepatosplenomegaly_Ids)
    
    umap_test2_2$Clinical <- "No"
    umap_test2_2$Clinical[umap_test2_2$UID %in% splenomegaly_Ids] <- "Splenomegaly"
    umap_test2_2$Clinical[umap_test2_2$UID %in% hepatomegaly_Ids] <- "Hepatomegaly"
    umap_test2_2$Clinical[umap_test2_2$UID %in% hepatosplenomegaly_Ids] <- "Hepatosplenomegaly"
    
    pdf(paste(output_plots,"UMAP_2dim_contries_short_clinical_hepatospleno.pdf", sep = "."), height = 5, width = 5)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                     breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")))
    dev.off()
    
    pdf(paste(output_plots,"UMAP_2dim_contries_short_clinical_hepatospleno.pdf", sep = "."), height = 3, width = 10)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                     breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")) +
            facet_wrap(~Class))
    dev.off()
    
    pdf(paste(output_plots,"UMAP_2dim_contries_short_clinical_hepatospleno2.pdf", sep = "."), height = 5, width = 4.5)
    print(ggplot(umap_test2_2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=3, alpha=0.5) + theme_bw() +
            geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
            geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                     breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")) +
            facet_wrap(~Class, ncol = 1))
    dev.off()
    
  }

  #Pheatmap
  to_pca_umap3_scaled_t <- data.frame(t(input_table2))
  
  
  tocol_column <- input_table_umap2[,c("Class"), drop =FALSE]
  rownames(tocol_column) <- input_table_umap2$UID
  my_colour = list(
    Class = c(HV = "#CC79A7", V1 = "#E69F00", V2 ="#56B4E9")
  )
  
  input_table_umap2_to_heatmap <- input_table_umap2[, !colnames(input_table_umap2) %in% c("UID", "Patient", "Class", "Country")]
  rownames(input_table_umap2_to_heatmap) <- input_table_umap2$UID
  input_table_umap2_to_heatmap2 <- t(input_table_umap2_to_heatmap)
  
  pdf(paste(output_plots,"heatmap_scaled_clinical.pdf", sep = "."), height = 9, width = 18)
  print(pheatmap(input_table_umap2_to_heatmap2, clustering_distance_cols = "manhattan", clustering_method = "average",
                 scale="row", annotation_col=tocol_column, annotation_colors=my_colour
  ))
  dev.off()
  
  if(length(data_table$Country) > 0 ) {
    pdf(paste(output_plots,"heatmap_scaled_large_clinical.pdf", sep = "."), height = 9, width = 22)
    print(pheatmap(input_table_umap2_to_heatmap2, clustering_distance_cols = "manhattan", clustering_method = "average",
                   fontsize_col=3, scale="row", annotation_col=tocol_column, annotation_colors=my_colour))
    dev.off()
    
  }
}

# This is the function that does the Logistic regression, Odd Ratio analysis:
mixomics_logReg_scalling_2_new <- function(input_table, test_trait, vector_traits_evalute, country, ncomp_value, markers_number, boot_number, outname) {
  
  directory_name <- paste(outname, test_trait, country, "test.dir", sep = ".")
  dir.create(directory_name)  
  
  #Subsetting the table to work with the citokine data
  class <- input_table[, test_trait, drop=FALSE]      # Only the last column
  
  features <- input_table[, colnames(input_table) %in% paste(c("UID", vector_traits_evalute))]  
  
  #Lognorm the data:
  df_log_transformed <- features  # copy original
  
  # Identify numeric columns (this excludes non-numeric ones)
  num_cols <- sapply(features, is.numeric)
  
  # Apply log(x + 0.001) only to numeric columns
  df_log_transformed[num_cols] <- lapply(df_log_transformed[num_cols], function(x) log(x + 0.001))
  
  features <- df_log_transformed
  
  #Scalling the table:
  features[,2:ncol(features)] <- scale(features[,2:ncol(features)], center = TRUE, scale = TRUE)
  
  rownames(features) <- features$UID
  features <- features[,colnames(features) != "UID",]
  
  # Convert features to a matrix and class to a vector
  feature_matrix <- as.matrix(features)
  class_vector <- class[,test_trait]
  
  X <- feature_matrix # use the gene expression data as the X matrix
  Y <- class_vector # use the class data as the Y matrix
  Y <-  factor(Y)
  
  ## Initial sPLS-DA
  ## Initial sPLS-DA with all data:
  all_samples.plsda <- mixOmics::splsda(X, Y, ncomp = ncomp_value, keepX =   length(vector_traits_evalute), scale = FALSE)  # set ncomp based on the number of traits. As we have few, 2 is enough
  
  # Get loadings for Component 1
  comp1_loadings_all <- loadings(all_samples.plsda)$X[, 1]
  
  # Obtain their importance:
  top_vars1_all <- data.frame(comp1_loadings_all)
  top_vars1_all$marker <- rownames(top_vars1_all)
  colnames(top_vars1_all) <- c("relevanceComp1", "marker")
  write.csv(top_vars1_all, paste(directory_name,"/", "sPLSDA.representative_components_full_set.csv", sep ="" ), row.names = FALSE)
  
  png(paste(directory_name,"/",test_trait, "_",  "sPLSDA_Loadings_1.png", sep ="" ), width = 2500, height = 2500, res = 300)
  plotLoadings(all_samples.plsda, comp = 1, method = 'mean', contrib = 'max')
  dev.off()
  
  #Selecting the top3 for the logist regression withthe full sample:
  top_vars1_all_top3_splsda <- head(rownames(top_vars1_all[order(abs(top_vars1_all$relevanceComp1), decreasing = TRUE),]),markers_number)
  
  #Running the univariate logistic regression for each trait:
  #To make the balanced bootstrap:
  #Add back the outcome:
  feature_matrix_df <- data.frame(feature_matrix)
  feature_matrix_df[[test_trait]] <- Y
  feature_matrix_df$age <- input_table$Age
  
  #Itertating trough all predictors
  predictors <- unique(colnames(feature_matrix))
  temp_df_for_logisti_regression <- data.frame()
  
  
  for (j in predictors) {
    
    formula_all_sing <- as.formula(paste(test_trait, " ~", j ))
    
    univariate_log <- logistf(formula_all_sing, data = feature_matrix_df, family = "binomial" , method = "logistf",  lcontrol = logistf.control(maxit = 1000))
    coef_table_uni <- data.frame(
      CI_lower = univariate_log$ci.lower,
      CI_uper = univariate_log$ci.upper,
      Coef = univariate_log$coefficients,
      Pvalue = univariate_log$prob
    )
    
    coef_table_uni <- coef_table_uni[rownames(coef_table_uni) != "(Intercept)",]
    
    temp_df_for_logisti_regression <- rbind(temp_df_for_logisti_regression, coef_table_uni)
    
  }
  temp_df_for_logisti_regression
  temp_df_for_logisti_regression$Class <- rownames(temp_df_for_logisti_regression)
  
  temp_df_for_logisti_regression <- temp_df_for_logisti_regression[order(temp_df_for_logisti_regression$Pvalue, decreasing = TRUE),]
  temp_df_for_logisti_regression$Class <- factor(temp_df_for_logisti_regression$Class, levels=temp_df_for_logisti_regression$Class)
  temp_df_for_logisti_regression$Sig <- "No"
  temp_df_for_logisti_regression$Sig[temp_df_for_logisti_regression$Pvalue<=0.05] <- "Yes"
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_fulldata.png", sep ="" ), width = 1200, height = 1500, res = 300)
  print(ggplot(temp_df_for_logisti_regression, aes(x=Class, y=Coef, colour=Pvalue, shape=Sig)) + geom_point(size = 5) + 
          geom_errorbar(aes(ymin = CI_lower, ymax = CI_uper), width = 0.1) + theme_bw() + coord_flip() +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) + scale_colour_gradient2(low = "red", mid="green", high ="blue", midpoint = 0.5) + scale_shape_manual(values = c(16,17), breaks = c("No", "Yes")))
  dev.off()
  
  
  # Refit standard logistic regression
  #Runnning the Bootstrap Analysis manually, for boeth sPLSda and :
  scoring_data_frame <- data.frame()
  relevance_values_data_frame <- data.frame()
  logistic_regression_all <-  data.frame()
  iterations_that_failed <- data.frame()
  
  markers_dataframe_iteration_1_all <- data.frame(Features=colnames(feature_matrix_df))
  markers_dataframe_iteration_1_all <- markers_dataframe_iteration_1_all[which(markers_dataframe_iteration_1_all$Features != test_trait), , drop=FALSE]
  
  #Separate the sample in affected and non-affected:
  input_table2 <- feature_matrix_df[, c(test_trait, vector_traits_evalute)]
  
  #For binaru classiffyrs
  input_table2_0 <- input_table2[which(input_table2[,test_trait] == 0),]
  input_table2_1 <- input_table2[which(input_table2[,test_trait] == 1),]
  
  sample_size_0 <- nrow(input_table2_0)
  sample_size_1 <- nrow(input_table2_1)
  
  #For reproductibility
  set.seed(123)
  temp_df_for_logisti_regression_boot <- data.frame()
  
  #Doing 100 bootstrap replicates
  for (i in 1:boot_number) {
    
    sample_0 <- input_table2_0[sample(1:nrow(input_table2_0),sample_size_0, replace =TRUE), ]
    sample_1 <- input_table2_1[sample(1:nrow(input_table2_1),sample_size_1, replace =TRUE), ]
    
    bootstraped_data <- rbind(sample_0, sample_1)
    
    bootstraped_data_2 <- bootstraped_data[, vector_traits_evalute]
    bootstraped_data_outcome <- bootstraped_data[,test_trait]
    
    #Running splssa with n markers
    plsda_model <- mixOmics::splsda(bootstraped_data_2, bootstraped_data_outcome, ncomp = ncomp_value,  keepX = c(markers_number), scale = FALSE) 
    
    #Storing the relevant values:
    vars_comp1 <- selectVar(plsda_model, comp = 1)  

    # Get loadings for Component 1
    comp1_loadings <- loadings(plsda_model)$X[, 1]

    # Obtain their importance:
    top_vars1 <- data.frame(comp1_loadings)
    top_vars1$marker <- rownames(top_vars1)
    colnames(top_vars1) <- c("relevanceComp1", "marker")
    
    relevance_values_data_frame <- rbind(relevance_values_data_frame, top_vars1)
    
    #Selecting the Topx from each:
    top_vars1_top3_splsda <- head(rownames(top_vars1[order(abs(top_vars1$relevanceComp1), decreasing = TRUE),]),markers_number)
    
    #Comp1    
    markers_dataframe_iteration_1 <- data.frame(Features = colnames(feature_matrix))
    markers_dataframe_iteration_1[[paste("Boot_rep", i, sep = "_")]] <- 0
    markers_dataframe_iteration_1[which(markers_dataframe_iteration_1$Features %in% top_vars1_top3_splsda),2] <- 1
    markers_dataframe_iteration_1_all <- left_join(markers_dataframe_iteration_1_all, markers_dataframe_iteration_1, by="Features")
    
    # Refit standard logistic regression
    x_train2 <- data.frame(bootstraped_data_2)
    x_train2[[test_trait]] <- bootstraped_data_outcome
    
    
    for (k in predictors) {
      
      formula_all_sing_boot <- as.formula(paste(test_trait, " ~", k ))
      
      fit <- tryCatch({
        
        univariate_log_boot <- logistf(formula_all_sing_boot, data = x_train2, family = "binomial" , method = "logistf",  lcontrol = logistf.control(maxit = 1000))
      }, error = function(e) {
        message("Skipping this iteration due to error: ", e$message)
        temp_failed <- data.frame(Iteration_number = paste("iteration", i, sep = "_"))
        iterations_that_failed <<- rbind(iterations_that_failed, temp_failed)
      })
      if (!is.null(fit)) {
        coef_table_uni_boot <- data.frame(
          CI_lower = univariate_log_boot$ci.lower,
          CI_uper = univariate_log_boot$ci.upper,
          Coef = univariate_log_boot$coefficients,
          Pvalue = univariate_log_boot$prob
        )
        
        coef_table_uni_boot$Class <- rownames(coef_table_uni_boot)
        coef_table_uni_boot <- coef_table_uni_boot[rownames(coef_table_uni_boot) != "(Intercept)",]
        coef_table_uni_boot$sample <- paste("Rep", i, sep = "_")
        
        temp_df_for_logisti_regression_boot <- rbind(temp_df_for_logisti_regression_boot, coef_table_uni_boot)
      }
      else {
        next  # skips to next iteration in the for loop
      }
      
    }
    
    
    
    
  }
  #Checking the relevance of each parameter in the iterations:
  top_values_comp1 <- data.frame(apply(markers_dataframe_iteration_1_all[2:ncol(markers_dataframe_iteration_1_all)], 1, sum))
  top_values_comp1$Features <- markers_dataframe_iteration_1_all$Features
  colnames(top_values_comp1) <- c("Comp1", "Features")
  
  #Plotting the number of times that each marker was in the top 5:
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Used_markersTOP3.png", sep ="" ), width = 1500, height = 1500, res = 300)
  print(ggplot(top_values_comp1, aes(x=Features, y=Comp1)) + geom_col() + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of occurences in the top3"))
  dev.off()
  
  write.csv(top_values_comp1, paste(directory_name,"/",test_trait, "_", country, ".", "Comp", "1e2", ".", "representative_components.csv", sep ="" ), row.names = FALSE)
  write.csv(markers_dataframe_iteration_1_all, paste(directory_name,"/",test_trait, "_", country, ".", "Comp", "1", ".", "representative_components.csv", sep ="" ), row.names = FALSE)

  #Relevance of each value to each replicate:
  write.csv(relevance_values_data_frame, paste(directory_name,"/",test_trait, "_", country, ".", "Comp", "1", ".", "trait_importance_each_iteration.csv", sep ="" ), row.names = FALSE)
  
  
  #Now for the results from the logistic regression
  write.csv(temp_df_for_logisti_regression_boot, paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression", ".csv", sep ="" ), row.names = TRUE)
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression.png", sep ="" ), width = 2500, height = 2200, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(x=-log10(Pvalue), y=Coef, color=Class)) + geom_point(alpha=0.5, size=2) +
          geom_vline(xintercept = 1.30103, color = "red", linetype = "dashed", size = 0.5) + theme_bw() +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + ylab("LR Coef") +
          xlab("-log10(p-value)"))
  dev.off()
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_wrap.png", sep ="" ), width = 2500, height = 2200, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(x=-log10(Pvalue), y=Coef, color=Class)) + geom_point(alpha=0.5, size=2) +
          geom_vline(xintercept = 1.30103, color = "red", linetype = "dashed", size = 0.5) + theme_bw() +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + ylab("LR Coef") +
          xlab("-log10(p-value)") + facet_wrap(~Class))
  dev.off()
  
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_wrap2.png", sep ="" ), width = 2500, height = 2200, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(x = -log10(Pvalue), y = Coef)) +
          geom_pointdensity(size = 2, alpha=0.5) +
          geom_vline(xintercept = 1.30103, color = "red", linetype = "dashed", size = 0.5) +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
          scale_color_viridis_c() +
          facet_wrap(~Class) +
          theme_bw() +
          ylab("LR Coef") +
          xlab("-log10(p-value)"))
  dev.off()
  
  
  options(scipen=999)
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_OR.png", sep ="" ), width = 2200, height = 1800, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(y=Class, x=Coef, fill=Class, colour=Class)) + geom_boxplot(alpha=0.5, size=1) +
          theme_bw() + geom_vline(xintercept = 0,  color = "red", linetype = "dashed", size = 0.5) + xlab("Coefficient"))
  dev.off()
  
  #For the Confidence Intervals:
  logistic_regression_all2_ci <- temp_df_for_logisti_regression_boot[, c("sample", "Class", "CI_lower", "CI_uper")]
  logistic_regression_all2_ci_melt <- melt(logistic_regression_all2_ci, id.vars = c("sample", "Class"))

  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_ORCI.png", sep ="" ), width = 2500, height = 2000, res = 300)
  print(ggplot(logistic_regression_all2_ci_melt, aes(y=variable, x=value, color=Class, group=sample)) + geom_point(size=1, alpha=0.7) +
          geom_line(size=0.1) +  facet_wrap(~Class, scales = "free") + theme_bw() + ylab("Confidence_Interval") +
          geom_vline(xintercept = 0,  color = "red", linetype = "dashed", size = 0.5) + xlab("Coefficient"))
  dev.off()
  
  if (nrow(iterations_that_failed) > 0) {
    write.csv(iterations_that_failed,
              paste(directory_name, "/", test_trait, "_", country, ".Iterations_that_failed.csv", sep = ""),
              row.names = FALSE)
  }
  
  
  return(logistic_regression_all2_ci_melt)
  
  
}
#Logistic regression. Same as above, but using age as a co variate:
mixomics_logReg_scalling_2_new_age_covariate <- function(input_table, test_trait, vector_traits_evalute, country, ncomp_value, markers_number, boot_number, outname) {
  
  directory_name <- paste(outname, test_trait, country, "test.dir", sep = ".")
  dir.create(directory_name)  
  
  #Subsetting the table to work with the citokine data
  class <- input_table[, test_trait, drop=FALSE]      # Only the last column
  
  features <- input_table[, colnames(input_table) %in% paste(c("UID", vector_traits_evalute))]  
  
  #Lognorm the data, without transforming age:
  features$Age <- as.character(features$Age)
  
  df_log_transformed <- features  # copy original
  
  # Identify numeric columns (this excludes non-numeric ones)
  num_cols <- sapply(features, is.numeric)
  
  # Apply log(x + 0.001) only to numeric columns
  df_log_transformed[num_cols] <- lapply(df_log_transformed[num_cols], function(x) log(x + 0.001))
  
  features <- df_log_transformed
  features$Age <- as.numeric(as.character(features$Age))
  
  
  #Scalling the table:
  features[,2:ncol(features)] <- scale(features[,2:ncol(features)], center = TRUE, scale = TRUE)
  
  rownames(features) <- features$UID
  features <- features[,colnames(features) != "UID",]
  
  # Convert features to a matrix and class to a vector
  feature_matrix <- as.matrix(features)
  #feature_matrix <- features
  class_vector <- class[,test_trait]
  
  # factor <- factor(class_vector)
  X <- feature_matrix # use the gene expression data as the X matrix
  Y <- class_vector # use the class data as the Y matrix
  Y <-  factor(Y)
  
  ## Initial sPLS-DA
  ## Initial sPLS-DA with all data:
  all_samples.plsda <- mixOmics::splsda(X, Y, ncomp = ncomp_value, keepX =   length(vector_traits_evalute), scale = FALSE)  
  
  # Get loadings for Component 1
  comp1_loadings_all <- loadings(all_samples.plsda)$X[, 1]
  
  # Obtain their importance:
  top_vars1_all <- data.frame(comp1_loadings_all)
  top_vars1_all$marker <- rownames(top_vars1_all)
  colnames(top_vars1_all) <- c("relevanceComp1", "marker")
  write.csv(top_vars1_all, paste(directory_name,"/", "sPLSDA.representative_components_full_set.csv", sep ="" ), row.names = FALSE)
  
  png(paste(directory_name,"/",test_trait, "_",  "sPLSDA_Loadings_1.png", sep ="" ), width = 2500, height = 2500, res = 300)
  plotLoadings(all_samples.plsda, comp = 1, method = 'mean', contrib = 'max')
  dev.off()
  
  #Selecting the top3 for the logist regression withthe full sample:
  top_vars1_all_top3_splsda <- head(rownames(top_vars1_all[order(abs(top_vars1_all$relevanceComp1), decreasing = TRUE),]),markers_number)
  
  #Running the univariate logistic regression for each trait:
  #To make the balanced bootstrap:
  #Add back the outcome:
  feature_matrix_df <- data.frame(feature_matrix)
  feature_matrix_df[[test_trait]] <- Y
  feature_matrix_df$age <- input_table$Age
  
  #Itertating trough all predictors
  predictors <- unique(colnames(feature_matrix))
  temp_df_for_logisti_regression <- data.frame()
  
  predictors <- predictors[predictors != "Age"]
  
  for (j in predictors) {
    
    #Including age as covariate
    formula_all_sing <- as.formula(paste(paste(test_trait, " ~", j ), " ", "Age", sep = " +" ))
    
    univariate_log <- logistf(formula_all_sing, data = feature_matrix_df, family = "binomial" , method = "logistf",  lcontrol = logistf.control(maxit = 1000))
    coef_table_uni <- data.frame(
      CI_lower = univariate_log$ci.lower,
      CI_uper = univariate_log$ci.upper,
      Coef = univariate_log$coefficients,
      Pvalue = univariate_log$prob
    )
    
    coef_table_uni <- coef_table_uni[!rownames(coef_table_uni) %in% c("(Intercept)", "Age"),]
    
    temp_df_for_logisti_regression <- rbind(temp_df_for_logisti_regression, coef_table_uni)
    
  }
  temp_df_for_logisti_regression
  temp_df_for_logisti_regression$Class <- rownames(temp_df_for_logisti_regression)
  
  temp_df_for_logisti_regression <- temp_df_for_logisti_regression[order(temp_df_for_logisti_regression$Pvalue, decreasing = TRUE),]
  temp_df_for_logisti_regression$Class <- factor(temp_df_for_logisti_regression$Class, levels=temp_df_for_logisti_regression$Class)
  temp_df_for_logisti_regression$Sig <- "No"
  temp_df_for_logisti_regression$Sig[temp_df_for_logisti_regression$Pvalue<=0.05] <- "Yes"
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_fulldata.png", sep ="" ), width = 1200, height = 1500, res = 300)
  print(ggplot(temp_df_for_logisti_regression, aes(x=Class, y=Coef, colour=Pvalue, shape=Sig)) + geom_point(size = 5) + 
          geom_errorbar(aes(ymin = CI_lower, ymax = CI_uper), width = 0.1) + theme_bw() + coord_flip() +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5) + scale_colour_gradient2(low = "red", mid="green", high ="blue", midpoint = 0.5) + scale_shape_manual(values = c(16,17), breaks = c("No", "Yes")))
  dev.off()
  
  
  # Refit standard logistic regression
  #Runnning the Bootstrap Analysis manually, for boeth sPLSda and :
  scoring_data_frame <- data.frame()
  relevance_values_data_frame <- data.frame()
  logistic_regression_all <-  data.frame()
  iterations_that_failed <- data.frame()
  
  markers_dataframe_iteration_1_all <- data.frame(Features=colnames(feature_matrix_df))
  markers_dataframe_iteration_1_all <- markers_dataframe_iteration_1_all[which(markers_dataframe_iteration_1_all$Features != test_trait), , drop=FALSE]
  
  #Separate the sample in affected and non-affected:
  input_table2 <- feature_matrix_df[, c(test_trait, vector_traits_evalute)]
  
  #For binaru classiffyrs
  input_table2_0 <- input_table2[which(input_table2[,test_trait] == 0),]
  input_table2_1 <- input_table2[which(input_table2[,test_trait] == 1),]
  
  sample_size_0 <- nrow(input_table2_0)
  sample_size_1 <- nrow(input_table2_1)
  
  #For reproductibility
  set.seed(123)
  temp_df_for_logisti_regression_boot <- data.frame()
  
  #Doing bootstrap replicates
  for (i in 1:boot_number) {
    
    sample_0 <- input_table2_0[sample(1:nrow(input_table2_0),sample_size_0, replace =TRUE), ]
    sample_1 <- input_table2_1[sample(1:nrow(input_table2_1),sample_size_1, replace =TRUE), ]
    
    bootstraped_data <- rbind(sample_0, sample_1)
    
    bootstraped_data_2 <- bootstraped_data[, vector_traits_evalute]
    bootstraped_data_outcome <- bootstraped_data[,test_trait]
    
    #Running splssa with n markers
    plsda_model <- mixOmics::splsda(bootstraped_data_2, bootstraped_data_outcome, ncomp = ncomp_value,  keepX = c(markers_number), scale = FALSE) 
    
    #Storing the relevant values:
    vars_comp1 <- selectVar(plsda_model, comp = 1)  
    # Get loadings for Component 1
    comp1_loadings <- loadings(plsda_model)$X[, 1]
    # Obtain their importance:
    top_vars1 <- data.frame(comp1_loadings)
    top_vars1$marker <- rownames(top_vars1)
    colnames(top_vars1) <- c("relevanceComp1", "marker")
    
    relevance_values_data_frame <- rbind(relevance_values_data_frame, top_vars1)
    
    #Selecting the Topx from each:
    top_vars1_top3_splsda <- head(rownames(top_vars1[order(abs(top_vars1$relevanceComp1), decreasing = TRUE),]),markers_number)
    
    #Comp1    
    markers_dataframe_iteration_1 <- data.frame(Features = colnames(feature_matrix))
    markers_dataframe_iteration_1[[paste("Boot_rep", i, sep = "_")]] <- 0
    markers_dataframe_iteration_1[which(markers_dataframe_iteration_1$Features %in% top_vars1_top3_splsda),2] <- 1
    markers_dataframe_iteration_1_all <- left_join(markers_dataframe_iteration_1_all, markers_dataframe_iteration_1, by="Features")
    
    # Refit standard logistic regression
    x_train2 <- data.frame(bootstraped_data_2)
    x_train2[[test_trait]] <- bootstraped_data_outcome
    
    
    for (k in predictors) {
      
      formula_all_sing_boot <- as.formula(paste(test_trait, " ~", k ))
      
      fit <- tryCatch({
        
        univariate_log_boot <- logistf(formula_all_sing_boot, data = x_train2, family = "binomial" , method = "logistf",  lcontrol = logistf.control(maxit = 1000))
      }, error = function(e) {
        message("Skipping this iteration due to error: ", e$message)
        temp_failed <- data.frame(Iteration_number = paste("iteration", i, sep = "_"))
        iterations_that_failed <<- rbind(iterations_that_failed, temp_failed)
      })
      if (!is.null(fit)) {
        coef_table_uni_boot <- data.frame(
          CI_lower = univariate_log_boot$ci.lower,
          CI_uper = univariate_log_boot$ci.upper,
          Coef = univariate_log_boot$coefficients,
          Pvalue = univariate_log_boot$prob
        )
        
        coef_table_uni_boot$Class <- rownames(coef_table_uni_boot)
        coef_table_uni_boot <- coef_table_uni_boot[rownames(coef_table_uni_boot) != "(Intercept)",]
        coef_table_uni_boot$sample <- paste("Rep", i, sep = "_")
        
        temp_df_for_logisti_regression_boot <- rbind(temp_df_for_logisti_regression_boot, coef_table_uni_boot)
      }
      else {
        next  # skips to next iteration in the for loop
      }
      
    }
    
    
    
    
  }
  #Checking the relevance of each parameter in the iterations:
  top_values_comp1 <- data.frame(apply(markers_dataframe_iteration_1_all[2:ncol(markers_dataframe_iteration_1_all)], 1, sum))
  top_values_comp1$Features <- markers_dataframe_iteration_1_all$Features
  colnames(top_values_comp1) <- c("Comp1", "Features")
  
  #Plotting the number of times that each marker was in the top 5:
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Used_markersTOP3.png", sep ="" ), width = 1500, height = 1500, res = 300)
  print(ggplot(top_values_comp1, aes(x=Features, y=Comp1)) + geom_col() + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of occurences in the top3"))
  dev.off()
  
  write.csv(top_values_comp1, paste(directory_name,"/",test_trait, "_", country, ".", "Comp", "1e2", ".", "representative_components.csv", sep ="" ), row.names = FALSE)
  write.csv(markers_dataframe_iteration_1_all, paste(directory_name,"/",test_trait, "_", country, ".", "Comp", "1", ".", "representative_components.csv", sep ="" ), row.names = FALSE)
  #Relevance of each value to each replicate:
  write.csv(relevance_values_data_frame, paste(directory_name,"/",test_trait, "_", country, ".", "Comp", "1", ".", "trait_importance_each_iteration.csv", sep ="" ), row.names = FALSE)
  
  
  #Now for the results from the logistic regression
  write.csv(temp_df_for_logisti_regression_boot, paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression", ".csv", sep ="" ), row.names = TRUE)
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression.png", sep ="" ), width = 2500, height = 2200, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(x=-log10(Pvalue), y=Coef, color=Class)) + geom_point(alpha=0.5, size=2) +
          geom_vline(xintercept = 1.30103, color = "red", linetype = "dashed", size = 0.5) + theme_bw() +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + ylab("LR Coef") +
          xlab("-log10(p-value)"))
  dev.off()
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_wrap.png", sep ="" ), width = 2500, height = 2200, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(x=-log10(Pvalue), y=Coef, color=Class)) + geom_point(alpha=0.5, size=2) +
          geom_vline(xintercept = 1.30103, color = "red", linetype = "dashed", size = 0.5) + theme_bw() +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + ylab("LR Coef") +
          xlab("-log10(p-value)") + facet_wrap(~Class))
  dev.off()
  
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_wrap2.png", sep ="" ), width = 2500, height = 2200, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(x = -log10(Pvalue), y = Coef)) +
          geom_pointdensity(size = 2, alpha=0.5) +
          geom_vline(xintercept = 1.30103, color = "red", linetype = "dashed", size = 0.5) +
          geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
          scale_color_viridis_c() +
          facet_wrap(~Class) +
          theme_bw() +
          ylab("LR Coef") +
          xlab("-log10(p-value)"))
  dev.off()
  
  
  options(scipen=999)
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_OR.png", sep ="" ), width = 2200, height = 1800, res = 300)
  print(ggplot(temp_df_for_logisti_regression_boot, aes(y=Class, x=Coef, fill=Class, colour=Class)) + geom_boxplot(alpha=0.5, size=1) +
          theme_bw() + geom_vline(xintercept = 0,  color = "red", linetype = "dashed", size = 0.5) + xlab("Coefficient"))
  dev.off()
  
  #For the Confidence Intervals:
  logistic_regression_all2_ci <- temp_df_for_logisti_regression_boot[, c("sample", "Class", "CI_lower", "CI_uper")]
  logistic_regression_all2_ci_melt <- melt(logistic_regression_all2_ci, id.vars = c("sample", "Class"))
  
  png(paste(directory_name,"/",test_trait, "_", country, ".", "Logistic_regression_ORCI.png", sep ="" ), width = 2500, height = 2000, res = 300)
  print(ggplot(logistic_regression_all2_ci_melt, aes(y=variable, x=value, color=Class, group=sample)) + geom_point(size=1, alpha=0.7) +
          geom_line(size=0.1) +  facet_wrap(~Class, scales = "free") + theme_bw() + ylab("Confidence_Interval") +
          geom_vline(xintercept = 0,  color = "red", linetype = "dashed", size = 0.5) + xlab("Coefficient"))
  dev.off()
  
  if (nrow(iterations_that_failed) > 0) {
    write.csv(iterations_that_failed,
              paste(directory_name, "/", test_trait, "_", country, ".Iterations_that_failed.csv", sep = ""),
              row.names = FALSE)
  }
  
  
  return(logistic_regression_all2_ci_melt)
  
  
}

# This function does the PLSDA analysis
# The PLSDA already scales the data internally
mixomics_sPLSDA_all <- function(input_table, column_to_evaluate, vector_traits_evalute, country, ncomp_value,  markers_number, outname) {
  set.seed(123) # for reproductibilitty
  
  directory_name <- paste(outname, column_to_evaluate, country, "mixo_sPLSDA.dir", sep = ".")
  dir.create(directory_name)  
  
  #Sub-setting the table to work with the cytokine data
  
  class <- input_table[, column_to_evaluate, drop=FALSE]
  
  features <- input_table[, colnames(input_table) %in% c("UID", vector_traits_evalute)]  
  
  rownames(features) <- features$UID
  features <- features[,colnames(features) != "UID",]
  
  # Identify numeric columns (this excludes non-numeric ones)
  num_cols <- sapply(features, is.numeric)
  features[num_cols] <- lapply(features[num_cols], function(x) log(x + 0.001))
  
  size_PCA <- length(features)
  # Convert features to a matrix and class to a vector
  feature_matrix <- features
  class_vector <- class[,column_to_evaluate]
  
  X <- feature_matrix # use the gene expression data as the X matrix
  Y <- class_vector # use the class data as the Y matrix
  Y <- as.factor(as.character(Y))
  
  ## Initial sPLS-DA
  ethiopia.plsda <- mixOmics::splsda(X, Y, ncomp = ncomp_value) 
  
  png(paste(directory_name,"/",column_to_evaluate, "_", country, ".", "sPLSDA_Loadings_1.png", sep ="" ), width = 2500, height = 2500, res = 300)
  plotLoadings(ethiopia.plsda, comp = 1, method = 'mean', contrib = 'max')
  dev.off()
  
  # Get loadings for Component 1
  comp1_loadings <- loadings(ethiopia.plsda)$X[, 1]
  
  # Rank them by absolute importance
  top_vars <- data.frame(sort(abs(comp1_loadings), decreasing = TRUE))
  top_vars$marker <- rownames(top_vars)
  colnames(top_vars) <- c("relevance", "marker")
  
  # View top 10 contributing variables, and remove the ones with low values less than 0.1
  top_vars2 <- top_vars[order(top_vars$relevance, decreasing = TRUE),]
  top_vars2 <- top_vars2[1:markers_number,]
  top_vars2 <- top_vars2[!is.na(top_vars2$marker),]
  top_vars$Class <- "Not_selected"
  top_vars$Class[top_vars$marker %in% top_vars2$marker] <- "Selected"
  
  write.csv(top_vars, paste(directory_name,"/",column_to_evaluate, "_", country, ".", "sPLSDA_Relevance.csv", sep ="" ))
  
  png(paste(directory_name,"/",column_to_evaluate, "_", country, ".", "PLSDA_relevance_scores.png", sep ="" ), width = 2000, height = 1500, res = 300)
  print(ggplot(top_vars, aes(x=marker, y=relevance, fill=Class)) + geom_col(alpha=0.8) + ylab("Relevance") +
          theme_bw() +  geom_hline(yintercept=0.1, linetype="dashed", color = "red", size=0.5, alpha= 0.7) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = c("#808080","#DC143C"), breaks = c("Not_selected","Selected")))
  dev.off()  
  
  vip_scores <- data.frame(vip(ethiopia.plsda))
  vip_scores$marker <- rownames(vip_scores)
  vip_scores$Class <- "Not_selected"
  vip_scores$Class[vip_scores$marker %in% top_vars2$marker] <- "Selected"
  
  png(paste(directory_name,"/",column_to_evaluate, "_", country, ".", "PLSDA_VIP_scores.png", sep ="" ), width = 2000, height = 1500, res = 300)
  print(ggplot(vip_scores, aes(x=marker, y=comp1, fill=Class)) + geom_col(alpha=0.8) + ylab("Relevance") +
          theme_bw() +  geom_hline(yintercept=1.3, linetype="dashed", color = "red", size=0.5, alpha= 0.7) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = c("#808080","#DC143C"), breaks = c("Not_selected","Selected")))
  dev.off()  
  
  #Testing a model with the data:
  #Selecting and generating the model based on the number of selected markers
  marker_names <- rownames(top_vars2)
  X_selected <- X[, marker_names ]
  marker_names_number <- length(marker_names)
  
  final.splsda2 <- mixOmics::splsda(X_selected, Y)
  
  png(paste(directory_name,"/",column_to_evaluate, "_", country, ".", "Comp", "1", ".", "AUC_1PC.png", sep ="" ), width = 2000, height = 2000, res = 300)
  auc.splsda = auroc(final.splsda2, roc.comp = 1, print = FALSE) # AUROC for the first component
  dev.off()
  
  png(paste(directory_name,"/",column_to_evaluate, "_", country, ".", "sPLSDA_Loadings_Final_1.png", sep ="" ), width = 2500, height = 2500, res = 300)
  plotLoadings(final.splsda2, comp = 1, method = 'mean', contrib = 'max')
  dev.off()
  
  
  #Plotting the violin plot of the selected groups:
  X_selected_long <- data.frame(X_selected)
  X_selected_long$Id <- rownames(X_selected_long)
  X_selected_long <- melt(X_selected_long, id.vars = c("Id"))
  X_selected_long$Class <- 0
  
  altered_samples <- input_table$UID[input_table[,column_to_evaluate] ==1]
  
  X_selected_long$Class[X_selected_long$Id %in% altered_samples ] <- 1
  X_selected_long$Class <- as.factor(as.character(X_selected_long$Class))
  
  png(paste(directory_name,"/",column_to_evaluate, "_", country, ".", ".", "Selected_markers", ".png", sep ="" ), width = 2000, height = 2000, res = 300)
  print(ggplot(X_selected_long, aes(x=Class, y=value, group=Class, fill=Class)) +  geom_violin(alpha=0.5, aes(fill=Class)) + 
          geom_boxplot(size=0.2, fill="white", width=0.2) + facet_wrap(~variable, scales = "free") + theme_bw() +
          scale_fill_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("1","0")) +
          stat_compare_means(method = "wilcox.test", size=2, color="#B22222") + ylab("") + xlab(""))
  dev.off()
  
  #Runnning the LOO Analysis manually:
  scoring_data_frame <- data.frame()
  
  i=1
  for (i in 1:nrow(X_selected)) {
    # Create the training set by removing the i-th sample (row)
    x_train <- X_selected[-i, , drop = FALSE]  # Remove the i-th sample
    y_train <- Y[-i]                  # Remove the i-th sample's outcome
    
    # Create the test set by selecting the i-th sample
    x_test <- X_selected[i, , drop = FALSE]    # The i-th sample is the test set
    y_test <- Y[i]                    # The true label for the i-th sample
    
    # Fit the PLS-DA model using the training set
    # train the model on the training data
    plsda_model <- mixOmics::splsda(x_train, y_train, ncomp = ncomp_value)  # Use appropriate ncomp
    
    # Make predictions for the i-th sample
    pred <- predict(plsda_model, newdata = x_test)$class  # Extract the predicted class
    
    # Store the prediction
    temp_df_classifyer <- data.frame(Patient = rownames(x_test),C1_MaxDist = pred$max.dist[1], C1_CentDist = pred$centroids.dist[1],  C1_MahaDist = pred$mahalanobis.dist[1],
                                     Real=y_test)
    scoring_data_frame <- rbind(scoring_data_frame, temp_df_classifyer)
    
    
  }
  
  
  scoring_data_frame$Real <- as.character(scoring_data_frame$Real)
  
  scoring_data_frame2 <- scoring_data_frame
  scoring_data_frame2[scoring_data_frame2=="1"] <- "Affected"
  scoring_data_frame2[scoring_data_frame2=="0"] <- "Not_Affected"
  
  # Manually compute the confusion matrix
  
  conf_matrix_MAH <- caret::confusionMatrix(factor(scoring_data_frame2$C1_MahaDist), factor(scoring_data_frame2$Real))
  conf_matrix_CEN <- caret::confusionMatrix(factor(scoring_data_frame2$C1_CentDist), factor(scoring_data_frame2$Real))
  conf_matrix_MAX <- caret::confusionMatrix(factor(scoring_data_frame2$C1_MaxDist), factor(scoring_data_frame2$Real))
  
  write.csv(scoring_data_frame2, paste(directory_name,"/",column_to_evaluate, "_", country, ".", "Comp", "1", ".", "Predictions_table_all.csv", sep ="" ), row.names = TRUE)
  
  cm_table <- as.table(conf_matrix_MAH$table)
  write.csv(cm_table, paste(directory_name,"/",column_to_evaluate, "_", country, ".", "Comp", "1", ".", "Conf_table_rowPred_Coul_Ref.csv", sep ="" ), row.names = TRUE)
  
  #as.table(conf_matrix_MAX$table)
  
  # Extract summary statistics of the confusion matrix
  cm_summary_MH <- as.data.frame(conf_matrix_MAH$byClass)
  colnames(cm_summary_MH) <- c("MAHDist_Results")
  cm_summary_CEN <- as.data.frame(conf_matrix_CEN$byClass)
  colnames(cm_summary_CEN) <- c("CENDist_Results")
  cm_summary_MAX <- as.data.frame(conf_matrix_MAX$byClass)
  colnames(cm_summary_MAX) <- c("MAXDist_Results")
  
  three_metrics <- cbind(cm_summary_MH, cm_summary_CEN, cm_summary_MAX)
  three_metrics_rounded <- as.data.frame(lapply(three_metrics, function(x) round(x, 2)))
  rownames(three_metrics_rounded) <- rownames(three_metrics)
  
  three_metrics_rounded$Features <- markers_number
  
  # Export to CSV
  write.csv(three_metrics_rounded, paste(directory_name,"/",column_to_evaluate, "_", country, ".", "Comp", "1", ".", "metrics.csv", sep ="" ), row.names = TRUE)
  
  return(three_metrics_rounded)
}
# This function automatizes the run of PLSDA, selecting 2,3,4,5,6,10 or 15 markers: 
run_splsda_several_values <- function(input_table, trait,  vector_traits_evalute, country, outname) {
  
  set.seed(123) # For reproductibility
  
  Temp_PLSDA_2 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute, country, 1, 2,  paste(outname, "run2", sep = ""))
  Temp_PLSDA_3 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute, country, 1, 3,  paste(outname, "run3", sep = ""))
  Temp_PLSDA_4 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute,  country, 1, 4,  paste(outname, "run4", sep = ""))
  Temp_PLSDA_5 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute,  country, 1, 5,  paste(outname, "run5", sep = ""))
  Temp_PLSDA_6 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute,  country, 1, 6,  paste(outname, "run6", sep = ""))
  Temp_PLSDA_10 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute,  country, 1, 10,  paste(outname, "run10", sep = ""))
  Temp_PLSDA_15 <- mixomics_sPLSDA_all(input_table, trait, vector_traits_evalute,  country, 1, 15,  paste(outname, "run15", sep = ""))
  Temp_PLSDA_2$Class <- "PLSDA"
  Temp_PLSDA_3$Class <- "PLSDA"
  Temp_PLSDA_4$Class <- "PLSDA"
  Temp_PLSDA_5$Class <- "PLSDA"
  Temp_PLSDA_6$Class <- "PLSDA"
  Temp_PLSDA_10$Class <- "PLSDA"
  Temp_PLSDA_15$Class <- "PLSDA"
 
  Temp_Hep_scores <- rbind(Temp_PLSDA_2, Temp_PLSDA_3, Temp_PLSDA_4, Temp_PLSDA_5, Temp_PLSDA_6, Temp_PLSDA_10, Temp_PLSDA_15)
  Temp_Hep_scores2 <- Temp_Hep_scores[,c("CENDist_Results", "Features", "Class")]
  Temp_Hep_scores2$Data <- rownames(Temp_Hep_scores2)
  
  Temp_Hep_scores2$Data <- sub("[0-9]+$", "", Temp_Hep_scores2$Data)
  Temp_Hep_scores2wide <- data.frame(pivot_wider(Temp_Hep_scores2, names_from = Data, values_from = CENDist_Results))
  Temp_Hep_scores2wide$Features <- as.character(Temp_Hep_scores2wide$Features)
  
  colnames(Temp_Hep_scores2wide)[colnames(Temp_Hep_scores2wide) =="F"] <- "F1"
  
  Temp_Hep_scores2wide$Features <- factor(Temp_Hep_scores2wide$Features , levels = c("2", "3","4","5","6","10","15"))
  
  A <- ggplot(Temp_Hep_scores2wide, aes(x=Features, y=Balanced.Accuracy, colour = Class, group=Class)) + geom_point() + ylim(0,1) + geom_line() + theme_bw() +
    ggtitle("Balanced Accuracy")
  
  B <- ggplot(Temp_Hep_scores2wide, aes(x=Features, y=Sensitivity, colour = Class, group=Class)) + geom_point() + ylim(0,1) + geom_line() + theme_bw() +
    ggtitle("Sensitivity")
  
  C <- ggplot(Temp_Hep_scores2wide, aes(x=Features, y=Specificity, colour = Class, group=Class)) + geom_point() + ylim(0,1) + geom_line() + theme_bw() +
    ggtitle("Specificity")
  
  D <- ggplot(Temp_Hep_scores2wide, aes(x=Features, y=Pos.Pred.Value, colour = Class, group=Class)) + geom_point() + ylim(0,1) + geom_line() + theme_bw() +
    ggtitle("Pos.Pred.Value")
  
  E <- ggplot(Temp_Hep_scores2wide, aes(x=Features, y=Neg.Pred.Value, colour = Class, group=Class)) + geom_point() + ylim(0,1) + geom_line() + theme_bw() +
    ggtitle("Neg.Pred.Value")
  
  G <- ggplot(Temp_Hep_scores2wide, aes(x=Features, y=F1, colour = Class, group=Class)) + geom_point() + ylim(0,1) + geom_line() + theme_bw() +
    ggtitle("F1")
  
  png(paste(outname, "png", sep = "."),  width = 3500,height = 2000, res = 300)
  print(plot_grid(A,B,C,D,E,G,  nrow=2))
  dev.off()
  
  
  #Plotting each separatelly, all metrics:
  Temp_Hep_scores2wide_Leaking <- Temp_Hep_scores2wide[which(Temp_Hep_scores2wide$Class=="PLSDA"),]

  png(paste(outname, "Leaking_predictions.png", sep ="" ), width = 1500, height = 1200, res = 300)
  print(ggplot(Temp_Hep_scores2wide_Leaking) + 
          geom_point(aes(x=Features, y=Sensitivity, colour = "Sensitivity", group=Class),  shape = 8, size=3)  + geom_line(aes(x=Features, y=Sensitivity, colour = "Sensitivity", group=Class)) +
          geom_point(aes(x=Features, y=Specificity, colour = "Specificity", group=Class),  shape = 17, size=3)  + geom_line(aes(x=Features, y=Specificity, colour = "Specificity", group=Class)) +
          geom_point(aes(x=Features, y=Pos.Pred.Value, colour = "PPV", group=Class), shape = 18, size=3)  + geom_line(aes(x=Features, y=Pos.Pred.Value, colour = "PPV", group=Class)) +
          geom_point(aes(x=Features, y=Neg.Pred.Value, colour = "NPV", group=Class), shape = 3, size=3)  + geom_line(aes(x=Features, y=Neg.Pred.Value, colour = "NPV", group=Class)) +
          geom_point(aes(x=Features, y=Balanced.Accuracy, colour = "Balanced.Accuracy", group=Class), shape = 16, size=3)  + geom_line(aes(x=Features, y=Balanced.Accuracy, colour = "Balanced.Accuracy", group=Class)) +
          geom_point(aes(x=Features, y=F1, colour = "F1", group=Class), shape = 4, size=3)  + geom_line(aes(x=Features, y=F1, colour = "F1", group=Class)) +
          scale_colour_manual(values = c("purple", "blue", "green", "brown", "red", "orange"), breaks = c("Sensitivity", "Specificity", "PPV", "NPV", "Balanced.Accuracy", "F1")) +
          ylim(0,1) + ylab("Score") + xlab("Number of features") +
          theme_bw())
  dev.off()
  
}

#########################################
## Evaluation the patient clinical data #
#########################################

#Import and process, clean and format the clinical data:
#Importing the clinical metadata:####
Ethiopia_data <- read.table("Ethiopia_metadata.table.csv", sep = ",", header = TRUE)
Kenya_data <- read.table("Kenya_metadata.table.csv", sep = ",", header = TRUE)
Sudan_data <- read.table("Sudan_metadata.table.csv", sep = ",", header = TRUE)
Uganda_data <- read.table("Uganda_metadata.table.csv", sep = ",", header = TRUE)

#correcting the Kenya weight in patient VL023:
Kenya_data$Weight[which(Kenya_data$Weight==0.4 & Kenya_data$Patient=="VL023")] <- 40

#Correcting IDs to match the citometry data:
Ethiopia_data$Patient <- gsub("HV1","HV0",Ethiopia_data$Patient)
Ethiopia_data$Patient <- gsub("P1","VL0",Ethiopia_data$Patient)
Kenya_data$Patient <- gsub("HV","HV0",Kenya_data$Patient)
Uganda_data$Patient <- gsub("A0","VL0", Uganda_data$Patient)
Uganda_data$Patient <- gsub("HC","HV", Uganda_data$Patient)
Sudan_data$Patient <- gsub("HC","HV", Sudan_data$Patient)

#formating the data:
Ethiopia_data[Ethiopia_data=="deceased"] <- NA
Ethiopia_data[Ethiopia_data=="NA"] <- NA
Kenya_data[Kenya_data=="NA"] <- NA
Uganda_data[Uganda_data=="NA"] <- NA
Sudan_data[Sudan_data=="NA"] <- NA
Kenya_data <- Kenya_data[,colnames(Kenya_data) != "ID"]

#Correcting different measurements:
#This was needed to standardize the measurements in some markers that were using different unities between sites
#Kenya Corrections
#1µmol/L=0.011312mg/dL
Kenya_data$Creatinine <- Kenya_data$Creatinine * 0.011312 #1g/L=0.1g/dL
Kenya_data$Albumin <- Kenya_data$Albumin * 0.1 #1 µmol/L of bilirubin=0.058467 mg/dL
Kenya_data$Total_Bilirubin <- Kenya_data$Total_Bilirubin * 0.058467
#Uganda Corrections
Uganda_data$Creatinine <- Uganda_data$Creatinine * 0.011312 # 1 µmol/L=0.011312mg/dL
Uganda_data$Total_Bilirubin <- Uganda_data$Total_Bilirubin * 0.058467 #1 µmol/L of bilirubin=0.058467 mg/dL

#Sudan Correction
Sudan_data$Height <- Sudan_data$Height/100 # cm to m

#Names corrections:
cols_to_rename <- c("Sex", "Hepatomegaly_York", "Splenomegaly_York", "SittingSystolicBloodPressure", "SittingDiastoliccBloodPressure")
new_names <- c("Males", "Hepatomegaly", "Splenomegaly", "SSystolicBP", "SDiastolicBP")
#Renaming Ethiopia
idx <- match(cols_to_rename, names(Ethiopia_data))
names(Ethiopia_data)[idx] <- new_names
#Renaming Kenya
idx <- match(cols_to_rename, names(Kenya_data))
names(Kenya_data)[idx] <- new_names
#Renaming Uganda
idx <- match(cols_to_rename, names(Uganda_data))
names(Uganda_data)[idx] <- new_names
#Renaming Sudan
idx <- match(cols_to_rename, names(Sudan_data))
names(Sudan_data)[idx] <- new_names

#Selecting the clinical metadata to be used in the downstream analysis
relevant_data <- c("Patient", "Class", "Age", "Height", "Weight", "SSystolicBP",
                   "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin",  "WBCells", "Neutrophil", "Aspirate_grade",
                   "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin", "Country", 
                   "Males","Hepatomegaly","Splenomegaly","Auxiliar_Lymphnodes", "Spleen_size")

#Structuring the patient metadata:
significant_plot_ethiopia2 <- Ethiopia_data
significant_plot_kenya2 <- Kenya_data
significant_plot_sudan2 <- Sudan_data
significant_plot_Uganda2 <- Uganda_data

significant_plot_ethiopia2$Country <- "Ethiopia"
significant_plot_kenya2$Country <- "Kenya"
significant_plot_sudan2$Country <- "Sudan"
significant_plot_Uganda2$Country <- "Uganda"

significant_plot_ethiopia3 <- significant_plot_ethiopia2[,colnames(significant_plot_ethiopia2) %in% relevant_data]
significant_plot_kenya3 <- significant_plot_kenya2[,colnames(significant_plot_kenya2) %in% relevant_data]
significant_plot_sudan3 <- significant_plot_sudan2[,colnames(significant_plot_sudan2) %in% relevant_data]
significant_plot_Uganda3 <- significant_plot_Uganda2[,colnames(significant_plot_Uganda2) %in% relevant_data]

#Combining all countries in a single data frame ####
Clinical_data_all <- rbind(significant_plot_ethiopia3, significant_plot_kenya3, significant_plot_sudan3, significant_plot_Uganda3)
Clinical_data_all$UID <- paste(Clinical_data_all$Patient, Clinical_data_all$Class, Clinical_data_all$Country, sep = "_")

#Correcting some labels that were misnamed:
Clinical_data_all_melt <- melt(Clinical_data_all, id.vars = c("Patient", "UID", "Class","Country"))
Clinical_data_all_melt$value[Clinical_data_all_melt$value =="NA"] <- NA
Clinical_data_all_melt$value[Clinical_data_all_melt$value =="NR"] <- NA
Clinical_data_all_melt$value[Clinical_data_all_melt$value =="NaN"] <- NA
Clinical_data_all_melt$value[Clinical_data_all_melt$value =="not recorded"] <- NA
Clinical_data_all_melt$value <- as.numeric(as.character(Clinical_data_all_melt$value))
###########################################################
# End of the Clinical data data import and processing #####
###########################################################


###################################
# Reading the LEGENDplex data: ####
###################################
# Reading and formatting the Legendplex tables: ####
Inflam_ethiopia1 <- read.table("Ethiophia_Inf_plat1.csv", sep = ",", header = TRUE)
Inflam_ethiopia2 <- read.table("Ethiophia_Inf_plat2.csv", sep = ",", header = TRUE)
Inflam_ethiopia3 <- read.table("Ethiophia_Inf_plat3.csv", sep = ",", header = TRUE)

Inflam_Kenya1 <- read.table("kenya_Inf_plat1.csv", sep = ",", header = TRUE)
Inflam_Kenya2 <- read.table("kenya_Inf_plat2.csv", sep = ",", header = TRUE)
Inflam_Kenya3 <- read.table("kenya_Inf_plat3.csv", sep = ",", header = TRUE)
Inflam_Kenya4 <- read.table("kenya_Inf_plat4.csv", sep = ",", header = TRUE)

Inflam_Uganda1 <- read.table("Uganda_Inf_plat1.csv", sep = ",", header = TRUE)
Inflam_Uganda2 <- read.table("Uganda_Inf_plat2.csv", sep = ",", header = TRUE)

Inflam_Sudan1 <- read.table("Sudan_Inf_plat1.csv", sep = ",", header = TRUE)
Inflam_Sudan2 <- read.table("Sudan_Inf_plat2.csv", sep = ",", header = TRUE)

#removing the wells that have the VL_pool 2 and HV_Pool 2: (02-Well-H12 for VL and 02-Well-F12 for HV)
#These are not usefull:
Inflam_Uganda2 <- Inflam_Uganda2[Inflam_Uganda2$sample != "02-Well-H12",]
Inflam_Uganda2 <- Inflam_Uganda2[Inflam_Uganda2$sample != "02-Well-F12",]

#Removing the erroneous duplicate value based on Karen's assessment:
#Kenya plate1
Inflam_Kenya1[Inflam_Kenya1$well =="02-Well-E5.fcs",  4:ncol(Inflam_Kenya1)] <- Inflam_Kenya1[Inflam_Kenya1$well =="02-Well-F5.fcs",  4:ncol(Inflam_Kenya1)]
Inflam_Kenya1[Inflam_Kenya1$well =="02-Well-G5.fcs",  4:ncol(Inflam_Kenya1)] <- Inflam_Kenya1[Inflam_Kenya1$well =="02-Well-H5.fcs", 4:ncol(Inflam_Kenya1)]
Inflam_Kenya1[Inflam_Kenya1$well =="02-Well-F8.fcs",  4:ncol(Inflam_Kenya1)] <- Inflam_Kenya1[Inflam_Kenya1$well =="02-Well-E8.fcs", 4:ncol(Inflam_Kenya1)]
##Kenya plate2
Inflam_Kenya2[Inflam_Kenya2$well =="02-Well-E10.fcs",  4:ncol(Inflam_Kenya2)] <- Inflam_Kenya2[Inflam_Kenya2$well =="02-Well-F10.fcs",  4:ncol(Inflam_Kenya2)]
Inflam_Kenya2[Inflam_Kenya2$well =="02-Well-A11.fcs",  4:ncol(Inflam_Kenya2)] <- Inflam_Kenya2[Inflam_Kenya2$well =="02-Well-B11.fcs", 4:ncol(Inflam_Kenya2)]

#Assigning the plate name for each plate, has to match the LOD information:
Inflam_ethiopia1$plate <- "plate1"
Inflam_ethiopia2$plate <- "plate2"
Inflam_ethiopia3$plate <- "plate3"

Inflam_Kenya1$plate <- "plate1" 
Inflam_Kenya2$plate <- "plate2" 
Inflam_Kenya3$plate <- "plate3" 
Inflam_Kenya4$plate <- "plate4" 

Inflam_Uganda1$plate <- "plate1" 
Inflam_Uganda2$plate <- "plate2" 

Inflam_Sudan1$plate <- "plate1" 
Inflam_Sudan2$plate <- "plate2" 

#Assigning the country:
Inflam_ethiopia1$class <- "Ethiopia"
Inflam_ethiopia2$class <- "Ethiopia"
Inflam_ethiopia3$class <- "Ethiopia"

Inflam_Kenya1$class <- "Kenya" 
Inflam_Kenya2$class <- "Kenya" 
Inflam_Kenya3$class <- "Kenya" 
Inflam_Kenya4$class <- "Kenya" 

Inflam_Uganda1$class <- "Uganda" 
Inflam_Uganda2$class <- "Uganda" 

Inflam_Sudan1$class <- "Sudan" 
Inflam_Sudan2$class <- "Sudan" 

# Correcting by the limit of detection amd combining plates from the same country ####
#Reading the limits of detection. The Plate name has to patch the plate ID
Ethiopia_LOD <- read.table("Ethiopia_LOQ.csv", sep = ",", header = TRUE)
Kenya_LOD <- read.table("Kenya_LOQ.csv", sep = ",", header = TRUE)
Uganda_LOD <- read.table("uganda_LOQ.csv", sep = ",", header = TRUE)
Sudan_LOD <- read.table("Sudan_LOQ.csv", sep = ",", header = TRUE)

#Correctin the limit of detection and names, and formating to the "long" format
Inflam_ethiopia1_corr <- correct_measurements_limit(Inflam_ethiopia1, Ethiopia_LOD, 2)
Inflam_ethiopia2_corr <- correct_measurements_limit(Inflam_ethiopia2, Ethiopia_LOD, 2)
Inflam_ethiopia3_corr <- correct_measurements_limit(Inflam_ethiopia3, Ethiopia_LOD, 2)

Inflam_Kenya1_corr <- correct_measurements_limit(Inflam_Kenya1, Kenya_LOD, 2)
Inflam_Kenya2_corr <- correct_measurements_limit(Inflam_Kenya2, Kenya_LOD, 2)
Inflam_Kenya3_corr <- correct_measurements_limit(Inflam_Kenya3, Kenya_LOD, 2)
Inflam_Kenya4_corr <- correct_measurements_limit(Inflam_Kenya4, Kenya_LOD, 2)


#Note: Uganda samples Patient 1 to 6 were diluted. All the others were non-diluted
Inflam_Uganda1_corr <- correct_measurements_limit_UGANDA(Inflam_Uganda1, Uganda_LOD, 1)
Inflam_Uganda2_corr <- correct_measurements_limit_UGANDA(Inflam_Uganda2, Uganda_LOD, 1)

Inflam_Sudan1_corr <- correct_measurements_limit(Inflam_Sudan1, Sudan_LOD, 2)
Inflam_Sudan2_corr <- correct_measurements_limit(Inflam_Sudan2, Sudan_LOD, 2)

#Combining the values of each individual country
Ethiopia_all <- rbind(Inflam_ethiopia1_corr, Inflam_ethiopia2_corr, Inflam_ethiopia3_corr)
Kenya_all <- rbind(Inflam_Kenya1_corr, Inflam_Kenya2_corr, Inflam_Kenya3_corr, Inflam_Kenya4_corr)
Uganda_all <- rbind(Inflam_Uganda1_corr, Inflam_Uganda2_corr)
Sudan_all <- rbind(Inflam_Sudan1_corr, Inflam_Sudan2_corr)

#Correcting the citokine names to avoid issues wirh R formating
Ethiopia_all$variable <- gsub("TGF.β1","TGF.B1", Ethiopia_all$variable)
Kenya_all$variable <- gsub("TGF.β1","TGF.B1", Kenya_all$variable)
Uganda_all$variable <- gsub("TGF.β1","TGF.B1", Uganda_all$variable)
Sudan_all$variable <- gsub("TGF.β1","TGF.B1", Sudan_all$variable)

#removing the markers that are not trustworthy based on report from Legendplex, and experts evaluation:
Ethiopia_all <- Ethiopia_all[!Ethiopia_all$variable %in% c("PAI.1","sCD130"),]
Kenya_all <- Kenya_all[!Kenya_all$variable %in% c("PAI.1","sCD130"),]
Uganda_all <- Uganda_all[!Uganda_all$variable %in% c("PAI.1","sCD130"),]
Sudan_all <- Sudan_all[!Sudan_all$variable %in% c("PAI.1","sCD130"),]

#Plotting the number of samples below the limit of detection:
plot_limit_of_detection(Ethiopia_all, "Ethiopia")
plot_limit_of_detection(Kenya_all, "Kenya")
plot_limit_of_detection(Uganda_all, "Uganda")
plot_limit_of_detection(Sudan_all, "Sudan")

#Generating the standardized results, for individual plates and the heatmap of the normalizers
Ethiopia_all_norm <- standardizing_plates(Ethiopia_all, "Ethiopia")
Kenya_all_norm <- standardizing_plates(Kenya_all, "Kenya")
Uganda_all_norm <- standardizing_plates(Uganda_all, "Uganda")
Sudan_all_norm <- standardizing_plates(Sudan_all, "Sudan")

#Exporting the normalized tables:
write.table(Ethiopia_all_norm, "Ethiopia_infl_all_norm.csv", sep = ",", row.names = FALSE, quote = FALSE)
write.table(Kenya_all_norm, "Kenya_infl_all_norm.csv", sep = ",", row.names = FALSE, quote = FALSE)
write.table(Uganda_all_norm, "Uganda_infl_all_norm.csv", sep = ",", row.names = FALSE, quote = FALSE)
write.table(Sudan_all_norm, "Sudan_infl_all_norm.csv", sep = ",", row.names = FALSE, quote = FALSE)

#Comparing the results for the duplicates
comparing_replicates(Ethiopia_all_norm, "Ethiopia")
comparing_replicates(Kenya_all_norm, "Kenya")
comparing_replicates(Uganda_all_norm, "Uganda")
comparing_replicates(Sudan_all_norm, "Sudan")

#Double checking the sample assignments:
ethiopia_IDS <- unique(Ethiopia_all_norm[,c("Sample_name", "group", "well", "plate", "class")])
kenya_IDS <- unique(Kenya_all_norm[,c("Sample_name", "group", "well", "plate", "class")])
uganda_IDS <- unique(Uganda_all_norm[,c("Sample_name", "group", "well", "plate", "class")])
sudan_IDS <- unique(Sudan_all_norm[,c("Sample_name", "group", "well", "plate", "class")])

allcountries_IDS <- rbind(ethiopia_IDS, kenya_IDS, uganda_IDS, sudan_IDS)
write.csv(allcountries_IDS, "immstat_Male_data.dir/All_samples_Ids.csv", quote = FALSE, row.names = FALSE)

###Cleaning the table and organizing for the analysis and plots
#Ethiophia
#Removing the controls:
Ethiopia_all_norm_to_plot <- Ethiopia_all_norm[!Ethiopia_all_norm$Sample_name %in% c("VL_Pool", "HV_Pool"),]
Ethiopia_all_norm_to_plot <- Ethiopia_all_norm_to_plot[!is.na(Ethiopia_all_norm_to_plot$Sample_name),]
Ethiopia_all_norm_to_plot_summarized <- data.frame(Ethiopia_all_norm_to_plot %>% group_by(Sample_name, group, variable) 
                                                   %>% summarise(Mean_rep=mean(norm_value), SD_rep=sd(norm_value)))
#Kenya:
#Removing the controls and summarizing the data:
Kenya_all_norm_to_plot <- Kenya_all_norm[!Kenya_all_norm$Sample_name %in% c("VL_Pool", "HV_Pool"),]
Kenya_all_norm_to_plot <- Kenya_all_norm_to_plot[!is.na(Kenya_all_norm_to_plot$Sample_name),]
Kenya_all_norm_to_plot_summarized <- data.frame(Kenya_all_norm_to_plot %>% group_by(Sample_name, group, variable) 
                                                %>% summarise(Mean_rep=mean(norm_value), SD_rep=sd(norm_value)))
#Uganda:
#Removing the controls and summarizing the data:
Uganda_all_norm_to_plot <- Uganda_all_norm[!Uganda_all_norm$Sample_name %in% c("VL_Pool", "HV_Pool"),]
Uganda_all_norm_to_plot <- Uganda_all_norm_to_plot[!is.na(Uganda_all_norm_to_plot$Sample_name),]
Uganda_all_norm_to_plot_summarized <- data.frame(Uganda_all_norm_to_plot %>% group_by(Sample_name, group, variable) 
                                                 %>% summarise(Mean_rep=mean(norm_value), SD_rep=sd(norm_value)))

#Sudan
Sudan_all_norm_to_plot <- Sudan_all_norm[!Sudan_all_norm$Sample_name %in% c("VL_Pool", "HV_Pool"),]
Sudan_all_norm_to_plot <- Sudan_all_norm_to_plot[!is.na(Sudan_all_norm_to_plot$Sample_name),]
Sudan_all_norm_to_plot_summarized <- data.frame(Sudan_all_norm_to_plot %>% group_by(Sample_name, group, variable) 
                                                 %>% summarise(Mean_rep=mean(norm_value), SD_rep=sd(norm_value)))

#Plotting the normalizer values for all countries in a single plot:
table_toplot_normalizers <- rbind(Ethiopia_all_norm, Kenya_all_norm_to_plot, Uganda_all_norm_to_plot, Sudan_all_norm_to_plot)
table_toplot_normalizers_unique <- unique(table_toplot_normalizers[,c("plate", "normalizer", "class", "variable")])

png(paste("Test", "Countries_boxplot_normalizer_LOD.png", sep = "_"), width = 3500,height = 3500, res = 300)
print(ggplot(table_toplot_normalizers_unique, aes(x=plate, y=normalizer, fill=class, group=class)) + geom_col(alpha=0.5, position = "dodge2") +  
        theme_bw() + ylab("Normalizer") + facet_wrap(~variable) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

png(paste("Test", "Countries_boxplot_normalizer_LOD_2.png", sep = "_"), width = 3500,height = 3500, res = 300)
print(ggplot(table_toplot_normalizers_unique, aes(x=class, y=normalizer, fill=plate, group=plate)) + geom_col(alpha=0.5, position = "dodge2") +  
        theme_bw() + ylab("Normalizer") + facet_wrap(~variable) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()


#Combining all data:
Ethiopia_all_norm_to_plot_summarized2 <- Ethiopia_all_norm_to_plot_summarized
Kenya_all_norm_to_plot_summarized2 <- Kenya_all_norm_to_plot_summarized
Uganda_all_norm_to_plot_summarized2 <- Uganda_all_norm_to_plot_summarized
Sudan_all_norm_to_plot_summarized2 <- Sudan_all_norm_to_plot_summarized

Ethiopia_all_norm_to_plot_summarized2$Country <- "Ethiopia"
Kenya_all_norm_to_plot_summarized2$Country <- "Kenya"
Uganda_all_norm_to_plot_summarized2$Country <- "Uganda"
Sudan_all_norm_to_plot_summarized2$Country <- "Sudan"

Ethiopia_all_norm_to_plot_summarized2$UID <- paste(Ethiopia_all_norm_to_plot_summarized2$Sample_name, Ethiopia_all_norm_to_plot_summarized2$group ,Ethiopia_all_norm_to_plot_summarized2$Country, sep = "_")
Kenya_all_norm_to_plot_summarized2$UID <- paste(Kenya_all_norm_to_plot_summarized2$Sample_name, Kenya_all_norm_to_plot_summarized2$group, Kenya_all_norm_to_plot_summarized2$Country, sep = "_")
Uganda_all_norm_to_plot_summarized2$UID <- paste(Uganda_all_norm_to_plot_summarized2$Sample_name, Uganda_all_norm_to_plot_summarized2$group, Uganda_all_norm_to_plot_summarized2$Country, sep = "_")
Sudan_all_norm_to_plot_summarized2$UID <- paste(Sudan_all_norm_to_plot_summarized2$Sample_name, Sudan_all_norm_to_plot_summarized2$group, Sudan_all_norm_to_plot_summarized2$Country, sep = "_")

All_countries_Legendplex_data <- rbind(Ethiopia_all_norm_to_plot_summarized2, Kenya_all_norm_to_plot_summarized2, Uganda_all_norm_to_plot_summarized2, Sudan_all_norm_to_plot_summarized2)


#Combining the full data from the metadata and Inflammation Panel: ####
#Expanding the tables to join:
Clinical_data_all_wide <- data.frame(pivot_wider(Clinical_data_all_melt, names_from = variable, values_from = value))
All_countries_Legendplex_data2 <- All_countries_Legendplex_data[,!colnames(All_countries_Legendplex_data) %in% c("SD_rep", "Country", "Sample_name", "group")]
All_countries_Legendplex_data_wide <- data.frame(pivot_wider(All_countries_Legendplex_data2, names_from = variable, values_from = Mean_rep))

#Obtaining a list of columns for each set:
Legendplex_data <- as.character(unique(All_countries_Legendplex_data2$variable))
Legendplex_data_ids <- c("Patient", "UID", "Class",  "Country", Legendplex_data)

Clinical_data <- as.character(unique(Clinical_data_all_melt$variable))
Clinical_data_ids <- as.vector(c("Patient", "UID", "Class",  "Country", Clinical_data))

#Combining the legendplex data and clinical metadata:
All_countries_clinical_Legendplex_data <- merge(Clinical_data_all_wide, All_countries_Legendplex_data_wide, by="UID", all = TRUE)
#Estimating the missing data in the full set:
all_data_counts <- data.frame(All_countries_clinical_Legendplex_data %>% group_by(Class, Country) %>% summarise(Count=n()))
all_data_counts2 <- data.frame(pivot_wider(all_data_counts, names_from = Class, values_from = Count))

#Removing any line with no patient ID:
All_countries_clinical_Legendplex_data <- All_countries_clinical_Legendplex_data[!is.na(All_countries_clinical_Legendplex_data$Patient),]


#Generating the table with the median and quantiles of each part of the data:
All_countries_clinical_Legendplex_data_melt <- melt(All_countries_clinical_Legendplex_data, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))
All_countries_clinical_Legendplex_data_melt <- 
  All_countries_clinical_Legendplex_data_melt[!All_countries_clinical_Legendplex_data_melt$variable %in% c("Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes"),]
All_countries_clinical_Legendplex_data_melt2 <- All_countries_clinical_Legendplex_data_melt[complete.cases(All_countries_clinical_Legendplex_data_melt),]

#Including in the table the units for the plot:
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Age","Age(years)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Height","Height(meter)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Weight","Weight(Kg)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("SSystolicBP","SSystolicBP(mmHg)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("SDiastolicBP","SDiastolicBP(mmHg)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Temperature","Temperature(°C)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Pulse","Pulse(Beats/min)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Hemoglobin","Hemoglobin(g/dl)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("WBCells","WBCells(10E9/L)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Neutrophil","Neutrophil(10E9/L)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Lymphocyte","Lymphocyte(10E9/L)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Platelets","Platelets(10E9/L)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("ALT_Results","ALT_Results(U/L)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Creatinine","Creatinine(mg/dl)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("Albumin","Albumin(g/dl)", All_countries_clinical_Legendplex_data_melt2$variable)

All_countries_clinical_Legendplex_data_melt2$variable <- gsub("CX3CL1","CX3CL1(pg/mL)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("CXCL12","CXCL12(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("PTX3","PTX3(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("TGF.B1","TGF.B1(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sCD25","sCD25(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sCD40L","sCD40L(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sRAGE","sRAGE(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sST2","sST2(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sTNF.RI","sTNF.RI(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sTNF.RII","sTNF.RII(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)
All_countries_clinical_Legendplex_data_melt2$variable <- gsub("sTREM.1","sTREM.1(pg/ml)", All_countries_clinical_Legendplex_data_melt2$variable)

All_countries_clinical_Legendplex_data_melt2_males <- All_countries_clinical_Legendplex_data_melt2[which(All_countries_clinical_Legendplex_data_melt2$Males==1),]
All_countries_clinical_Legendplex_data_melt2_females <- All_countries_clinical_Legendplex_data_melt2[which(All_countries_clinical_Legendplex_data_melt2$Males==0),]

#Males:
All_countries_clinical_Legendplex_data_melt2_males_summary <- data.frame(All_countries_clinical_Legendplex_data_melt2_males %>% 
                group_by(Class,  Country, variable) %>% 
                  summarise(Min= round(min(value),2), Q1 = round(quantile(value, 0.25, na.rm = TRUE),2), Median = round(median(value),2), 
                            Q3 = round(quantile(value, 0.75, na.rm = TRUE),2), Max= round(max(value),2), Patients = n()) )

write.csv(All_countries_clinical_Legendplex_data_melt2_males_summary, "immstat_Male_data.dir/All_data_males.overall.table.csv", quote = FALSE, row.names = FALSE)


#Females:
All_countries_clinical_Legendplex_data_melt2_females_summary <- data.frame(All_countries_clinical_Legendplex_data_melt2_females %>% 
                                                                           group_by(Class, Country, variable) %>% 
                                                                           summarise(Min= round(min(value),2), Q1 = round(quantile(value, 0.25, na.rm = TRUE),2), Median = round(median(value),2), 
                                                                                     Q3 = round(quantile(value, 0.75, na.rm = TRUE),2), Max= round(max(value),2),  Patients = n()))

write.csv(All_countries_clinical_Legendplex_data_melt2_females_summary, "immstat_female_data.dir/All_data_females.overall.table.csv",  quote = FALSE, row.names = FALSE)



#Subseting the data on Males and females:
#Males tables ####
All_countries_data_males <- All_countries_clinical_Legendplex_data[which(All_countries_clinical_Legendplex_data$Males==1),]

All_countries_data_males_melt <- melt(All_countries_data_males, id.vars = c("UID", "Patient", "Class",  "Country", "Males", "Age", "Height", "Weight",
                                                                            "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes", "Aspirate_grade", "Spleen_size"))


#Plotting the differences in markers V1 an V2:
#Subtracting V2 from V1
All_countries_clinical_Legendplex_data_melt <- melt(All_countries_clinical_Legendplex_data, id.vars = c("UID", "Patient", "Class",  "Country", "Males", "Age", "Height", "Weight",
                                                                            "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes", "Aspirate_grade", "Spleen_size"))

All_countries_clinical_Legendplex_data_to_diferences_V2_V1 <- All_countries_clinical_Legendplex_data_melt[which(All_countries_clinical_Legendplex_data_melt$Class != "HV"),
                                                                                   colnames(All_countries_clinical_Legendplex_data_melt) %in% c("Patient", "Class", "Country",
                                                                                                                                                      "variable", "value")]

All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide <- data.frame(pivot_wider(All_countries_clinical_Legendplex_data_to_diferences_V2_V1, names_from = "Class", values_from = "value"))
All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide <- All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide[complete.cases(All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide),]

All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide$V2minusV1 <- All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide$V2-All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide$V1 
All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide$NormDiff <- All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide$V2minusV1/All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide$V1


png(paste("immstat_Male_data.dir/Alldata", "Countries_V1V2_diff.boxplot.png", sep = "_"), width = 6000,height = 3500, res = 300)
ggplot(All_countries_clinical_Legendplex_data_to_diferences_V2_V1_wide, aes(x=Country, y=NormDiff, fill =Country)) + geom_point(alpha=0.7) + facet_wrap(~variable, scales = "free") + 
  coord_flip() + geom_boxplot() + theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=0.5, alpha= 0.7)
dev.off()

#Male tables ####
All_countries_data_males_ethiopia <- All_countries_data_males[which(All_countries_data_males$Country=="Ethiopia"),]
All_countries_data_males_kenya <- All_countries_data_males[which(All_countries_data_males$Country=="Kenya"),]
All_countries_data_males_uganda <- All_countries_data_males[which(All_countries_data_males$Country=="Uganda"),]
All_countries_data_males_sudan <- All_countries_data_males[which(All_countries_data_males$Country=="Sudan"),]

#Females tables ####
All_countries_data_females <- All_countries_clinical_Legendplex_data[All_countries_clinical_Legendplex_data$Males==0,]
All_countries_data_females_ethiopia <- All_countries_data_females[which(All_countries_data_females$Country=="Ethiopia"),]
All_countries_data_females_kenya <- All_countries_data_females[which(All_countries_data_females$Country=="Kenya"),]
All_countries_data_females_uganda <- All_countries_data_females[which(All_countries_data_females$Country=="Uganda"),]
All_countries_data_females_sudan <- All_countries_data_females[which(All_countries_data_females$Country=="Sudan"),]

# Working with male data: ####
#Note!!! For the Uganda clinical discrete data, Patient V2 22 has no data (only sex).
#For this reason it is 15 rather than 16 for the whole plot.

#Estimating the missing data
estimate_missing_data(All_countries_data_males, "immstat_Male_data.dir/Missing_all_males")
estimate_missing_data(All_countries_clinical_Legendplex_data, "immstat_Male_data.dir/Missing_all_data")

#Poting the Clinical data distribution
plot_continuous_and_disctreta_data(All_countries_data_males, "immstat_Male_data.dir/All_males_clinical")
#For All data:
plot_continuous_and_disctreta_data(All_countries_clinical_Legendplex_data, "immstat_Male_data.dir/All_samples_clinical")
#For the females
plot_continuous_and_disctreta_data(All_countries_data_females, "immstat_female_data.dir/All_females_clinical")

#Creating directories to store the statistical data analysis:
dir.create("immstat_Male_data.dir/Statistical_analysis_all")
dir.create("immstat_Male_data.dir/Statistical_analysis_all/All")
dir.create("immstat_Male_data.dir/Statistical_analysis_all/Ethiopia")
dir.create("immstat_Male_data.dir/Statistical_analysis_all/Kenya")
dir.create("immstat_Male_data.dir/Statistical_analysis_all/Sudan")
dir.create("immstat_Male_data.dir/Statistical_analysis_all/Uganda")

#Statistical comparison 
stats_all <- statistical_comparison_all(All_countries_data_males, "immstat_Male_data.dir/Statistical_analysis_all/All/All_males_statistical")
stats_ethiopia <- statistical_comparison_all(All_countries_data_males_ethiopia, "immstat_Male_data.dir/Statistical_analysis_all/Ethiopia/Ethiopia_males_statistical")
stats_kenya <- statistical_comparison_all(All_countries_data_males_kenya, "immstat_Male_data.dir/Statistical_analysis_all/Kenya/Kenya_males_statistical")
stats_uganda <- statistical_comparison_all(All_countries_data_males_uganda, "immstat_Male_data.dir/Statistical_analysis_all/Sudan/Uganda_males_statistical")
stats_sudan <- statistical_comparison_all(All_countries_data_males_sudan, "immstat_Male_data.dir/Statistical_analysis_all/Uganda/Sudan_males_statistical")

#Creating directories to store the PCA data analysis:
dir.create("immstat_Male_data.dir/PCA_analysis_all")
#dir.create("immstat_Male_data.dir/PCA_analysis_all/All")
dir.create("immstat_Male_data.dir/PCA_analysis_all/Ethiopia")
dir.create("immstat_Male_data.dir/PCA_analysis_all/Kenya")
dir.create("immstat_Male_data.dir/PCA_analysis_all/Sudan")
dir.create("immstat_Male_data.dir/PCA_analysis_all/Uganda")

#PCA of the clinical data:
pca_and_umap_funcion(All_countries_data_males_ethiopia, "immstat_Male_data.dir/PCA_analysis_all/Ethiopia/Ethiopia_males_PCA", 15)
pca_and_umap_funcion(All_countries_data_males_kenya, "immstat_Male_data.dir/PCA_analysis_all/Kenya/Kenya_males_PCA", 15)
pca_and_umap_funcion(All_countries_data_males_uganda, "immstat_Male_data.dir/PCA_analysis_all/Uganda/Uganda_males_PCA", 8)
pca_and_umap_funcion(All_countries_data_males_sudan, "immstat_Male_data.dir/PCA_analysis_all/Sudan/Sudan_males_PCA", 8)

#Correlation plot

#Creating directories to store the PCA data analysis:
# dir.create("immstat_Male_data.dir/PCA_analysis_all")
dir.create("immstat_Male_data.dir/Cor_analysis_all")
#dir.create("immstat_Male_data.dir/Cor_analysis_all/All")
dir.create("immstat_Male_data.dir/Cor_analysis_all/Ethiopia")
dir.create("immstat_Male_data.dir/Cor_analysis_all/Kenya")
dir.create("immstat_Male_data.dir/Cor_analysis_all/Sudan")
dir.create("immstat_Male_data.dir/Cor_analysis_all/Uganda")

sample_order_to_plot <- c("Age", "Height", "Weight", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes", 
                          "SSystolicBP", "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", 
                          "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin", "Aspirate_grade", "Spleen_size", 
                          "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", "sRAGE", "sST2", "sTNF.RI", "sTNF.RII", "sTREM.1" )

#Correlation using only V1 or V2, for each country:
#Ethiopia
All_countries_data_males_ethiopia_V1 <- All_countries_data_males_ethiopia[All_countries_data_males_ethiopia$Class=="V1",]
All_countries_data_males_ethiopia_V2 <- All_countries_data_males_ethiopia[All_countries_data_males_ethiopia$Class=="V2",]
ethiopiaV1_coor <- correlation_plots_with_metadata(All_countries_data_males_ethiopia_V1, sample_order_to_plot, "Ethiopia", "immstat_Male_data.dir/Cor_analysis_all/Ethiopia/Ethiopia_V1_males_correlation")
ethiopiaV2_coor <- correlation_plots_with_metadata(All_countries_data_males_ethiopia_V2, sample_order_to_plot, "Ethiopia", "immstat_Male_data.dir/Cor_analysis_all/Ethiopia/Ethiopia_V2_males_correlation")
#Kenya
# View(All_countries_data_males_kenya)
All_countries_data_males_kenia_V1 <- All_countries_data_males_kenya[All_countries_data_males_kenya$Class=="V1",]
All_countries_data_males_kenia_V2 <- All_countries_data_males_kenya[All_countries_data_males_kenya$Class=="V2",]
kenyaV1_coor <- correlation_plots_with_metadata(All_countries_data_males_kenia_V1, sample_order_to_plot, "Kenya", "immstat_Male_data.dir/Cor_analysis_all/Kenya/Kenya_V1_males_correlation")
kenyaV2_coor <- correlation_plots_with_metadata(All_countries_data_males_kenia_V2, sample_order_to_plot, "Kenya", "immstat_Male_data.dir/Cor_analysis_all/Kenya/Kenya_V2_males_correlation")

#Uganda:
All_countries_data_males_uganda_V1 <- All_countries_data_males_uganda[All_countries_data_males_uganda$Class=="V1",]
All_countries_data_males_uganda_V2 <- All_countries_data_males_uganda[All_countries_data_males_uganda$Class=="V2",]
ugandaV1_coor <- correlation_plots_with_metadata(All_countries_data_males_uganda_V1, sample_order_to_plot, "Uganda", "immstat_Male_data.dir/Cor_analysis_all/Uganda/Uganda_V1_males_correlation")
ugandaV2_coor <- correlation_plots_with_metadata(All_countries_data_males_uganda_V2, sample_order_to_plot, "Uganda", "immstat_Male_data.dir/Cor_analysis_all/Uganda/Uganda_V2_males_correlation")

#Sudan
All_countries_data_males_sudan_V1 <- All_countries_data_males_sudan[All_countries_data_males_sudan$Class=="V1",]
All_countries_data_males_sudan_V2 <- All_countries_data_males_sudan[All_countries_data_males_sudan$Class=="V2",]
sudanV1_coor <- correlation_plots_with_metadata(All_countries_data_males_sudan_V1, sample_order_to_plot, "Sudan", "immstat_Male_data.dir/Cor_analysis_all/Sudan/Sudan_V1_males_correlation")
sudanV2_coor <-correlation_plots_with_metadata(All_countries_data_males_sudan_V2, sample_order_to_plot, "Sudan", "immstat_Male_data.dir/Cor_analysis_all/Sudan/Sudan_V2_males_correlation")

#Comparing the values of clinical traits with reference values for each country
compare_clinical_results_to_referece(All_countries_data_males, "immstat_Male_data.dir/Males_Compare_standard")
compare_clinical_results_to_referece(All_countries_clinical_Legendplex_data, "immstat_Male_data.dir/All_data")
compare_clinical_results_to_referece(All_countries_data_females, "immstat_female_data.dir/Fenale_Compare_standard")


#Female data for Kenya, Ugganda and Ethiopia:

#Creating directories to store the statistical data analysis:
dir.create("immstat_female_data.dir/Statistical_analysis_all")
# dir.create("immstat_female_data.dir/Statistical_analysis_all/All")
# dir.create("immstat_female_data.dir/Statistical_analysis_all/Ethiopia")
dir.create("immstat_female_data.dir/Statistical_analysis_all/Kenya")
dir.create("immstat_female_data.dir/Statistical_analysis_all/Sudan")
dir.create("immstat_female_data.dir/Statistical_analysis_all/Uganda")

#Creating directories to store the PCA data analysis:
dir.create("immstat_female_data.dir/PCA_analysis_all")
#dir.create("immstat_Male_data.dir/PCA_analysis_all/All")
# dir.create("immstat_female_data.dir/PCA_analysis_all/Ethiopia")
dir.create("immstat_female_data.dir/PCA_analysis_all/Kenya")
dir.create("immstat_female_data.dir/PCA_analysis_all/Sudan")
dir.create("immstat_female_data.dir/PCA_analysis_all/Uganda")


stats_kenya_fem <- statistical_comparison_all_noHV_nostats_inPlot(All_countries_data_females_kenya, "immstat_female_data.dir/Statistical_analysis_all/Kenya/Kenya_females_statistical")
stats_uganda_fem <- statistical_comparison_all(All_countries_data_females_uganda, "immstat_female_data.dir/PCA_analysis_all/Uganda/Uganda_females_statistical")
stats_sudan_fem <- statistical_comparison_all(All_countries_data_females_sudan, "immstat_female_data.dir/Statistical_analysis_all/Sudan/Sudan_females_statistical")
pca_and_umap_funcion(All_countries_data_females_kenya, "immstat_female_data.dir/PCA_analysis_all/Kenya/Kenya_females_PCA", 7)
pca_and_umap_funcion(All_countries_data_females_uganda, "immstat_female_data.dir/PCA_analysis_all/Uganda/Uganda_females_PCA", 7)
pca_and_umap_funcion(All_countries_data_females_sudan, "immstat_female_data.dir/PCA_analysis_all/Sudan/Sudan_females_PCA", 7)

#Comparing males and females

#Creating directories to store the PCA data analysis:
dir.create("immstat_compare_males_females.dir")
dir.create("immstat_compare_males_females.dir/Kenya")
dir.create("immstat_compare_males_females.dir/Sudan")
dir.create("immstat_compare_males_females.dir/Uganda")



All_countries_clinical_Legendplex_data_uganda <- All_countries_clinical_Legendplex_data[which(All_countries_clinical_Legendplex_data$Country == "Uganda"),]
All_countries_clinical_Legendplex_data_uganda_V1 <- All_countries_clinical_Legendplex_data_uganda[which(All_countries_clinical_Legendplex_data_uganda$Class=="V1"),]
All_countries_clinical_Legendplex_data_uganda_V2 <- All_countries_clinical_Legendplex_data_uganda[which(All_countries_clinical_Legendplex_data_uganda$Class=="V2"),]
All_countries_clinical_Legendplex_data_sudan <- All_countries_clinical_Legendplex_data[which(All_countries_clinical_Legendplex_data$Country == "Sudan"),]
All_countries_clinical_Legendplex_data_sudan_V1 <- All_countries_clinical_Legendplex_data_sudan[which(All_countries_clinical_Legendplex_data_sudan$Class=="V1"),]
All_countries_clinical_Legendplex_data_sudan_V2 <- All_countries_clinical_Legendplex_data_sudan[which(All_countries_clinical_Legendplex_data_sudan$Class=="V2"),]
All_countries_clinical_Legendplex_data_kenya <- All_countries_clinical_Legendplex_data[which(All_countries_clinical_Legendplex_data$Country == "Kenya"),]
All_countries_clinical_Legendplex_data_kenya_V1 <- All_countries_clinical_Legendplex_data_kenya[which(All_countries_clinical_Legendplex_data_kenya$Class=="V1"),]
All_countries_clinical_Legendplex_data_kenya_V2 <- All_countries_clinical_Legendplex_data_kenya[which(All_countries_clinical_Legendplex_data_kenya$Class=="V2"),]

statsitical_compare_plot_group(All_countries_clinical_Legendplex_data_uganda_V1, "Males", "immstat_compare_males_females.dir/Uganda/Uganda_male_female_V1")
statsitical_compare_plot_group(All_countries_clinical_Legendplex_data_uganda_V2, "Males", "immstat_compare_males_females.dir/Uganda/Uganda_male_female_V2")
statsitical_compare_plot_group(All_countries_clinical_Legendplex_data_sudan_V1, "Males", "immstat_compare_males_females.dir/Sudan/sudan_male_female_V1")
statsitical_compare_plot_group(All_countries_clinical_Legendplex_data_sudan_V2, "Males", "immstat_compare_males_females.dir/Sudan/sudan_male_female_V2")
statsitical_compare_plot_group(All_countries_clinical_Legendplex_data_kenya_V1, "Males", "immstat_compare_males_females.dir/Kenya/kenya_male_female_V1")
statsitical_compare_plot_group(All_countries_clinical_Legendplex_data_kenya_V2, "Males", "immstat_compare_males_females.dir/Kenya/kenya_male_female_V2")

#Checking the patients that Relapsed:
ethiopia_relapse <- c("VL006","VL016")
All_countries_data_males_ethiopia_V1_relapse <- All_countries_data_males_ethiopia_V1
All_countries_data_males_ethiopia_V2_relapse <- All_countries_data_males_ethiopia_V2

All_countries_data_males_ethiopia_V1_relapse$Class[All_countries_data_males_ethiopia_V1_relapse$UID %in% c("VL006_V1_Ethiopia","VL016_V1_Ethiopia")] <- "Relapse"
All_countries_data_males_ethiopia_V2_relapse$Class[All_countries_data_males_ethiopia_V2_relapse$UID %in% c("VL006_V2_Ethiopia","VL016_V2_Ethiopia")] <- "Relapse"

All_countries_data_males_ethiopia_V1_relapse_melt <- melt(All_countries_data_males_ethiopia_V1_relapse, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))
All_countries_data_males_ethiopia_V2_relapse_melt <- melt(All_countries_data_males_ethiopia_V2_relapse, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))

All_countries_data_males_ethiopia_relapse <- All_countries_data_males_ethiopia
All_countries_data_males_ethiopia_relapse$Class[All_countries_data_males_ethiopia_relapse$UID %in% c("VL006_V1_Ethiopia","VL016_V1_Ethiopia")] <- "Relapse"
All_countries_data_males_ethiopia_relapse$Class[All_countries_data_males_ethiopia_relapse$UID %in% c("VL006_V2_Ethiopia","VL016_V2_Ethiopia")] <- "Relapse"

All_countries_data_males_ethiopia_relapse_melt <- melt(All_countries_data_males_ethiopia, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))

All_countries_data_males_ethiopia_relapse_melt2 <- All_countries_data_males_ethiopia_relapse_melt[which(All_countries_data_males_ethiopia_relapse_melt$Class != "HV"),]
All_countries_data_males_ethiopia_relapse_melt2$Relapse <- "No"
All_countries_data_males_ethiopia_relapse_melt2$Relapse[All_countries_data_males_ethiopia_relapse_melt2$UID %in% c("VL006_V1_Ethiopia","VL016_V1_Ethiopia",
                                                                                                                   "VL006_V2_Ethiopia","VL016_V2_Ethiopia")] <- "Yes"

#To reorder samples of interest to the front:
# Assume you want groups "A" and "B" to be in front
front_groups <- c("VL006", "VL016")
All_countries_data_males_ethiopia_relapse_melt2$Patient <- factor(All_countries_data_males_ethiopia_relapse_melt2$Patient)  # Make sure it's a factor for consistency
# Reorder: non-front groups first, front groups last
All_countries_data_males_ethiopia_relapse_melt2_ordered <- All_countries_data_males_ethiopia_relapse_melt2[!(All_countries_data_males_ethiopia_relapse_melt2$Patient %in% front_groups), ]
All_countries_data_males_ethiopia_relapse_melt2_ordered2 <- rbind(All_countries_data_males_ethiopia_relapse_melt2_ordered,
                                                                  All_countries_data_males_ethiopia_relapse_melt2[All_countries_data_males_ethiopia_relapse_melt2$Patient %in% front_groups, ])

All_countries_data_males_ethiopia_relapse_melt2_ordered2$Relapse[which(All_countries_data_males_ethiopia_relapse_melt2_ordered2$Patient =="VL006")] <- "VL006"
All_countries_data_males_ethiopia_relapse_melt2_ordered2$Relapse[which(All_countries_data_males_ethiopia_relapse_melt2_ordered2$Patient =="VL016")] <- "VL016"

All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot <- 
  All_countries_data_males_ethiopia_relapse_melt2_ordered2[which(!All_countries_data_males_ethiopia_relapse_melt2_ordered2$variable %in% 
                                                                    c("Age", "Height", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes", "Aspirate_grade",
                                                                      "SSystolicBP","SDiastolicBP","Temperature", "Pulse", "Weight")),]
#To correct: Include the dashed lines values:
dashed_lines <- data.frame(variable = c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin"),
  yintercept = c(13, 4, 1.5, 1.0, 150, 10, 0.7, 3.5, 0.1 )  # Different y-values for each facet
)

dashed_lines_max <- data.frame(variable = c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin"),
  yintercept = c(17, 11, 8, 3.5, 400, 40, 1.3, 5.5, 1.2 )  # Different y-values for each facet
)
dashed_lines$variable <- factor(dashed_lines$variable, c( "Albumin", "Hemoglobin",
                                                          "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                          "Lymphocyte"))
dashed_lines_max$variable <- factor(dashed_lines_max$variable, c( "Albumin", "Hemoglobin",
                                                                  "Platelets", "Creatinine", "Total_Bilirubin", "ALT_Results", "WBCells", "Neutrophil", 
                                                                  "Lymphocyte"))

png(paste("immstat_Male_data.dir/Ethiopia_relapse_continuous.png", sep = "_"), width = 2300,height = 2300, res = 300)
ggplot(All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot, aes(x=Class, y=value, color=Relapse, group=Patient, shape=Relapse)) + geom_point(size=3) +
        geom_line(size=0.1) +  facet_wrap(~variable, scales = "free_y", ncol = 3) + scale_shape_manual(values = c(15,16,1), breaks = c("VL006", "VL016", "No")) + theme_bw() +
        geom_hline(data = dashed_lines, aes(yintercept = yintercept),  linetype = "dashed", color = "#DC143C", size = 0.7) +
        geom_hline(data = dashed_lines_max, aes(yintercept = yintercept),  linetype = "dashed", color = "#228B22", size = 0.7) +
  scale_colour_manual(values = c("red", "orange", "grey"), breaks = c("VL006", "VL016", "No")) + scale_y_log10() + ylab("")
dev.off()

All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot2 <- 
  All_countries_data_males_ethiopia_relapse_melt2_ordered2[which(All_countries_data_males_ethiopia_relapse_melt2_ordered2$variable %in% 
        c("Aspirate_grade")), , drop = F]

All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot2 <- 
  All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot2[which(All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot2$Class=="V1"),]

png(paste("immstat_Male_data.dir/Ethiopia_relapse_discrete.png", sep = "_"), width = 800,height = 2300, res = 300)
ggplot(All_countries_data_males_ethiopia_relapse_melt2_ordered2_to_plot2, aes(x=value, y=Patient, fill=Relapse, group=Relapse, shape=Relapse)) + geom_col() +
  facet_wrap(~variable, scales = "free", ncol = 3) + scale_fill_manual(values = c("red", "orange", "gray"), breaks = c("VL006", "VL016", "No")) + theme_bw()
dev.off()

# Now, the relapse from patient from Sudan:
All_countries_data_males_sudan_V1_relapse <- All_countries_data_males_sudan_V1
All_countries_data_males_sudan_V2_relapse <- All_countries_data_males_sudan_V2
All_countries_data_males_sudan_V1_relapse$Class[All_countries_data_males_sudan_V1_relapse$UID %in% c("VL017_V1_Sudan")] <- "Relapse"
All_countries_data_males_sudan_V2_relapse$Class[All_countries_data_males_sudan_V2_relapse$UID %in% c("VL017_V2_Sudan")] <- "Relapse"
All_countries_data_males_sudan_V1_relapse_melt <- melt(All_countries_data_males_sudan_V1_relapse, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))
All_countries_data_males_sudan_V2_relapse_melt <- melt(All_countries_data_males_sudan_V2_relapse, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))
All_countries_data_males_sudan_relapse <- All_countries_data_males_sudan
All_countries_data_males_sudan_relapse$Class[All_countries_data_males_sudan_relapse$UID %in% c("VL017_V1_Sudan")] <- "Relapse"
All_countries_data_males_sudan_relapse$Class[All_countries_data_males_sudan_relapse$UID %in% c("VL017_V2_Sudan")] <- "Relapse"
All_countries_data_males_sudan_relapse_melt <- melt(All_countries_data_males_sudan, id.vars = c("UID", "Patient", "Class",  "Country", "Males"))
All_countries_data_males_sudan_relapse_melt2 <- All_countries_data_males_sudan_relapse_melt[which(All_countries_data_males_sudan_relapse_melt$Class != "HV"),]
All_countries_data_males_sudan_relapse_melt2$Relapse <- "No"
All_countries_data_males_sudan_relapse_melt2$Relapse[All_countries_data_males_sudan_relapse_melt2$UID %in% c("VL017_V1_Sudan","VL017_V2_Sudan")] <- "Yes"


#To reorder samples of interest to the front:
# Assume you want groups "A" and "B" to be in front
front_groups <- c("VL017")
All_countries_data_males_sudan_relapse_melt2$Patient <- factor(All_countries_data_males_sudan_relapse_melt2$Patient)  # Make sure it's a factor for consistency
# Reorder: non-front groups first, front groups last
All_countries_data_males_sudan_relapse_melt2_ordered <- All_countries_data_males_sudan_relapse_melt2[!(All_countries_data_males_sudan_relapse_melt2$Patient %in% front_groups), ]
All_countries_data_males_sudan_relapse_melt2_ordered2 <- rbind(All_countries_data_males_sudan_relapse_melt2_ordered,
                                                                  All_countries_data_males_sudan_relapse_melt2[All_countries_data_males_sudan_relapse_melt2$Patient %in% front_groups, ])
All_countries_data_males_sudan_relapse_melt2_ordered2$Relapse[which(All_countries_data_males_sudan_relapse_melt2_ordered2$Patient =="VL017")] <- "VL017"
All_countries_data_males_sudan_relapse_melt2_ordered2_to_plot <- 
  All_countries_data_males_sudan_relapse_melt2_ordered2[which(!All_countries_data_males_sudan_relapse_melt2_ordered2$variable %in% 
                                                                   c("Age", "Height", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes", "Aspirate_grade",
                                                                     "SSystolicBP","SDiastolicBP","Temperature", "Pulse", "Weight")),]

png(paste("immstat_Male_data.dir/sudan_relapse_continuous.png", sep = "_"), width = 2300,height = 2300, res = 300)
ggplot(All_countries_data_males_sudan_relapse_melt2_ordered2_to_plot, aes(x=Class, y=value, color=Relapse, group=Patient, shape=Relapse)) + geom_point(size=3) +
  geom_line(size=0.1) +  facet_wrap(~variable, scales = "free_y", ncol = 3) + scale_shape_manual(values = c(15,1), breaks = c("VL017", "No")) + theme_bw() +
  geom_hline(data = dashed_lines, aes(yintercept = yintercept),  linetype = "dashed", color = "#DC143C", size = 0.7) +
  geom_hline(data = dashed_lines_max, aes(yintercept = yintercept),  linetype = "dashed", color = "#228B22", size = 0.7) +
  scale_colour_manual(values = c("blue", "grey"), breaks = c("VL017", "No")) + scale_y_log10() + ylab("")
dev.off()

##########################
#PCA and UMAP all samples#
#########################
#Generating the UMAP for the individually scaled data:

All_countries_data_males_ethiopia_scaled <- All_countries_data_males_ethiopia
All_countries_data_males_kenya_scaled <- All_countries_data_males_kenya
All_countries_data_males_uganda_scaled <- All_countries_data_males_uganda
All_countries_data_males_sudan_scaled <- All_countries_data_males_sudan

#Scaling the data for the PCA and UMAP with all samples:
All_countries_data_males_ethiopia_scaled[,12:ncol(All_countries_data_males_ethiopia_scaled)] <-lapply(All_countries_data_males_ethiopia_scaled[,12:ncol(All_countries_data_males_ethiopia_scaled)], function(x) log(x + 0.001))
All_countries_data_males_kenya_scaled[,12:ncol(All_countries_data_males_kenya_scaled)] <-lapply(All_countries_data_males_kenya_scaled[,12:ncol(All_countries_data_males_kenya_scaled)], function(x) log(x + 0.001))
All_countries_data_males_uganda_scaled[,12:ncol(All_countries_data_males_uganda_scaled)] <-lapply(All_countries_data_males_uganda_scaled[,12:ncol(All_countries_data_males_uganda_scaled)], function(x) log(x + 0.001))
All_countries_data_males_sudan_scaled[,12:ncol(All_countries_data_males_sudan_scaled)] <-lapply(All_countries_data_males_sudan_scaled[,12:ncol(All_countries_data_males_sudan_scaled)], function(x) log(x + 0.001))
#Scale:
All_countries_data_males_ethiopia_scaled[,12:ncol(All_countries_data_males_ethiopia_scaled)] <- scale(All_countries_data_males_ethiopia_scaled[,12:ncol(All_countries_data_males_ethiopia_scaled)], center = TRUE, scale = TRUE)
All_countries_data_males_kenya_scaled[,12:ncol(All_countries_data_males_kenya_scaled)] <- scale(All_countries_data_males_kenya_scaled[,12:ncol(All_countries_data_males_kenya_scaled)], center = TRUE, scale = TRUE)
All_countries_data_males_uganda_scaled[,12:ncol(All_countries_data_males_uganda_scaled)] <- scale(All_countries_data_males_uganda_scaled[,12:ncol(All_countries_data_males_uganda_scaled)], center = TRUE, scale = TRUE)
All_countries_data_males_sudan_scaled[,12:ncol(All_countries_data_males_sudan_scaled)] <- scale(All_countries_data_males_sudan_scaled[,12:ncol(All_countries_data_males_sudan_scaled)], center = TRUE, scale = TRUE)
v1v2all_countries_scaled <- rbind(All_countries_data_males_ethiopia_scaled, All_countries_data_males_kenya_scaled, All_countries_data_males_uganda_scaled, All_countries_data_males_sudan_scaled)

### Now for the UMAP with all markers:
v1v2all_countries_scaled_toUMAP_all <- v1v2all_countries_scaled[c("UID", "Patient", "Class", "Country", 
                                                                  "SSystolicBP", "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin", "WBCells", "Neutrophil",
                                                                  "Lymphocyte", "Platelets", "ALT_Results",  "Creatinine", "Albumin", "Total_Bilirubin",
                                                                  "CX3CL1", "CXCL12", "PTX3", "sCD25", "sCD40L", "sRAGE", "sST2",
                                                                  "sTNF.RI", "sTNF.RII", "sTREM.1")]

#Removing NAs
v1v2all_countries_scaled_toUMAP_all <- v1v2all_countries_scaled_toUMAP_all[complete.cases(v1v2all_countries_scaled_toUMAP_all),]
rownames(v1v2all_countries_scaled_toUMAP_all) <- v1v2all_countries_scaled_toUMAP_all$UID
v1v2all_countries_scaled_toUMAP_all2 <-  v1v2all_countries_scaled_toUMAP_all[c("SSystolicBP", "SDiastolicBP", "Temperature", "Pulse", "Hemoglobin", "WBCells", "Neutrophil",
                                                                               "Lymphocyte", "Platelets", "ALT_Results",  "Creatinine", "Albumin", "Total_Bilirubin",
                                                                               "CX3CL1", "CXCL12", "PTX3", "sCD25", "sCD40L", "sRAGE", "sST2",
                                                                               "sTNF.RI", "sTNF.RII", "sTREM.1")]

v1v2all_countries_scaled_toPCA_all2 <- v1v2all_countries_scaled_toUMAP_all2
rownames(v1v2all_countries_scaled_toPCA_all2) <- v1v2all_countries_scaled_toUMAP_all$UID
#Acessing the number of samples used in the PCA and UMAP:
nrow(v1v2all_countries_scaled_toPCA_all2) #206
v1v2all_countries_scaled_toUMAP_all %>% group_by(Class, Country) %>% summarise(count=n())
# HV    Ethiopia    26
# HV    Kenya       28
# HV    Sudan       10
# V1    Ethiopia    31
# V1    Kenya       31
# V1    Sudan       13
# V1    Uganda       9
# V2    Ethiopia    26
# V2    Kenya       21
# V2    Sudan       11



#PCA all data:
res.pca_all <- prcomp(v1v2all_countries_scaled_toPCA_all2, scale = FALSE)
pdf("immstat_Male_data.dir/PCAs_relevance_all.pdf", height = 5, width = 5)
print(fviz_eig(res.pca_all))
dev.off()
pdf("immstat_Male_data.dir/relevance_each_marker_all.pdf", height = 5, width = 5)
print(fviz_cos2(res.pca_all, choice = "var", axes = 1:2))
dev.off()
#Supporting table:
v1v2all_countries_scaled_toUMAP_all_temp <- v1v2all_countries_scaled_toUMAP_all
rownames(v1v2all_countries_scaled_toUMAP_all_temp) <- v1v2all_countries_scaled_toUMAP_all_temp$UID
pdf("immstat_Male_data.dir/PCA2_dim_biplot_all_Inf.pdf", height = 8, width = 8)
print(autoplot(res.pca_all, size=5, data = v1v2all_countries_scaled_toUMAP_all_temp, colour = 'Class', shape="Country", loadings = TRUE, loadings.label = TRUE, alpha=0.7,
               loadings.label.size = 3, loadings.colour = 'blue') + 
        scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")) +
        theme_bw())
dev.off()

#UMAP all data
set.seed(123)
v1v2all_countries_scaled_toUMAP_all2est1 <- umap::umap(v1v2all_countries_scaled_toUMAP_all2, n_epochs=1000, n_neighbors=15,  scale = FALSE)
v1v2all_countries_scaled_toUMAP_all2est2 <- data.frame(v1v2all_countries_scaled_toUMAP_all2est1$layout)
colnames(v1v2all_countries_scaled_toUMAP_all2est2) <- c("dim1", "dim2")
v1v2all_countries_scaled_toUMAP_all2est2$Class <- v1v2all_countries_scaled_toUMAP_all$Class
v1v2all_countries_scaled_toUMAP_all2est2$UID <- v1v2all_countries_scaled_toUMAP_all$UID
v1v2all_countries_scaled_toUMAP_all2est2$Coutry <- v1v2all_countries_scaled_toUMAP_all2est2$UID
v1v2all_countries_scaled_toUMAP_all2est2$Coutry <- gsub(".*_","",v1v2all_countries_scaled_toUMAP_all2est2$Coutry)

#Umap_short:
pdf("immstat_Male_data.dir/UMAP_2dim_contries_short_all.pdf", height = 5, width = 5)
print(ggplot(v1v2all_countries_scaled_toUMAP_all2est2, aes(x=dim1, y=dim2, color=Class, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
        geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
        geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), breaks = c("V1","V2","HV")))
dev.off()

hepatomegaly_Ids <- v1v2all_countries_scaled$UID[which(v1v2all_countries_scaled$Hepatomegaly == 1)]
splenomegaly_Ids <- v1v2all_countries_scaled$UID[which(v1v2all_countries_scaled$Splenomegaly == 1)]
hepatosplenomegaly_Ids <- v1v2all_countries_scaled$UID[which(v1v2all_countries_scaled$Hepatomegaly == 1 & v1v2all_countries_scaled$Splenomegaly == 1) ]
hepatomegaly_Ids <- setdiff(hepatomegaly_Ids, hepatosplenomegaly_Ids)
splenomegaly_Ids <- setdiff(splenomegaly_Ids, hepatosplenomegaly_Ids)
v1v2all_countries_scaled_toUMAP_all2est2$Clinical <- "No"
v1v2all_countries_scaled_toUMAP_all2est2$Clinical[v1v2all_countries_scaled_toUMAP_all2est2$UID %in% splenomegaly_Ids] <- "Splenomegaly"
v1v2all_countries_scaled_toUMAP_all2est2$Clinical[v1v2all_countries_scaled_toUMAP_all2est2$UID %in% hepatomegaly_Ids] <- "Hepatomegaly"
v1v2all_countries_scaled_toUMAP_all2est2$Clinical[v1v2all_countries_scaled_toUMAP_all2est2$UID %in% hepatosplenomegaly_Ids] <- "Hepatosplenomegaly"

pdf("immstat_Male_data.dir/UMAP_2dim_contries_short_immune_hepatospleno_all.pdf", height = 5, width = 5)
print(ggplot(v1v2all_countries_scaled_toUMAP_all2est2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
        geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
        geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                 breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")))
dev.off()

pdf("immstat_Male_data.dir/UMAP_2dim_contries_short_immune_hepatospleno_all_split.pdf", height = 3, width = 10)
print(ggplot(v1v2all_countries_scaled_toUMAP_all2est2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=4, alpha=0.5) + theme_bw() +
        geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
        geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                 breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")) +
        facet_wrap(~Class))
dev.off()

pdf("immstat_Male_data.dir/UMAP_2dim_contries_short_immune_hepatospleno2_all.pdf", height = 5, width = 4.5)
print(ggplot(v1v2all_countries_scaled_toUMAP_all2est2, aes(x=dim1, y=dim2, color=Clinical, label = UID, shape=Coutry)) + geom_point(size=3, alpha=0.5) + theme_bw() +
        geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7) + theme(legend.title=element_blank()) +
        geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.5, alpha= 0.7)  + scale_color_manual(values = c("#B8860B", "#56B4E9", "#FF6347", "darkgray"), 
                                                                                                                 breaks = c("Splenomegaly","Hepatomegaly","Hepatosplenomegaly", "No")) +
        facet_wrap(~Class, ncol = 1))
dev.off()


#Saving the session:
save.image(file = "Seession_prior_PLSDA.RData")


#sPLS-DA analysis:
#https://mixomicsteam.github.io/mixOmics-Vignette/id_05.html#id_05:splsda
#Removing missing data from Hepatomegaly

#Creating a directory for each comparison analysis
dir.create("immstat_Male_data.dir/Ethiopia_HepatoV1_PLSDA")
dir.create("immstat_Male_data.dir/Kenya_SplenoV2_PLSDA")  #Splenomegaly in Kenya on V2
dir.create("immstat_Male_data.dir/Ethiopia_SplenoV2_PLSDA") # Splenomegaly in Ethiopia on V2
dir.create("immstat_Male_data.dir/Ethiopia_SplenoV1toV2_PLSDA") #Using the samples from V1 to predc the outcome in V2 in splenomegaly
dir.create("immstat_Male_data.dir/Kenya_SplenoV1toV2_PLSDA") #Kenya using V1 to predict recovery from splenomegaly in Kenya
dir.create("immstat_Male_data.dir/Sudan_hepV1_PLSDA") #Sudan hepatonomegaly 
dir.create("immstat_Male_data.dir/Sudan_SplenoV1_PLSDA") #Sudan splenomegaly 
dir.create("immstat_Male_data.dir/Uganda_splenV1_PLSDA") #uganda splenomegaly

#Removing unwanted columns, and moving the categorical columns to the beginning:
selec_and_order_traits <- c( "UID", "Patient", "Class", "Country", "Males", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes",
                             "Age", "Height", "Weight",  "SSystolicBP", "SDiastolicBP", "Temperature",  "Pulse", "Hemoglobin",
                             "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                             "Aspirate_grade", "Spleen_size", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", 
                             "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")

#Traits that will be evaluated:
vector_traits_evalute <- c("CX3CL1", "CXCL12", "PTX3", "sCD25", "sCD40L", "sRAGE", "sST2", "sTNF.RI", "sTNF.RII", "sTREM.1")
ids_to_keep <- c( "UID", "Patient", "Class", "Country", "Males", "Hepatomegaly", "Splenomegaly")

vector_traits_evalute2 <-  c( "Hemoglobin",
                             "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                             "CX3CL1", "CXCL12", "PTX3", "sCD25", "sCD40L",
                             "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")

vector_traits_evalute2_age  <- c(vector_traits_evalute2, "Age")
vector_traits_evalute2_AG <- c(vector_traits_evalute2, "Aspirate_grade")
vector_traits_evalute2_age_AG <- c(vector_traits_evalute2_age, "Aspirate_grade")

#############################
#Ethiopia splenomegaly in V1#
#############################
#Clinical traits:
clinical_traits_evaluate <- c("Hemoglobin", "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin")
All_countries_data_males_ethiopia_reorder <- All_countries_data_males_ethiopia[, selec_and_order_traits]
All_countries_data_males_ethiopia_reorder_clean <- All_countries_data_males_ethiopia_reorder[which(All_countries_data_males_ethiopia_reorder$Class != "HV"),]
All_countries_data_males_ethiopia_reorder_clean <- All_countries_data_males_ethiopia_reorder_clean[,!colnames(All_countries_data_males_ethiopia_reorder_clean) %in%
                                                                                                     c("Auxiliar_Lymphnodes")]
#Selecting V1
All_countries_data_males_ethiopia_reorder_clean_V1 <-
  All_countries_data_males_ethiopia_reorder_clean[All_countries_data_males_ethiopia_reorder_clean$Class=="V1",]
#Citokines
All_countries_data_males_ethiopia_reorder_clean_V1_citokines <- All_countries_data_males_ethiopia_reorder_clean_V1[,
                                        c(ids_to_keep,vector_traits_evalute)]
All_countries_data_males_ethiopia_reorder_clean_V1_citokines <- All_countries_data_males_ethiopia_reorder_clean_V1_citokines[
  complete.cases(All_countries_data_males_ethiopia_reorder_clean_V1_citokines),]
#Clinical markers:
All_countries_data_males_ethiopia_reorder_clean_V1_clinical <- 
  All_countries_data_males_ethiopia_reorder_clean_V1[, c(ids_to_keep,clinical_traits_evaluate)]
All_countries_data_males_ethiopia_reorder_clean_V1_clinical <- All_countries_data_males_ethiopia_reorder_clean_V1_clinical[
  complete.cases(All_countries_data_males_ethiopia_reorder_clean_V1_clinical),]
#All
All_countries_data_males_ethiopia_reorder_clean_V1_noNA_all <- All_countries_data_males_ethiopia_reorder_clean_V1[
  complete.cases(All_countries_data_males_ethiopia_reorder_clean_V1),]
#Plotting the difference between groups:
statsitical_compare_plot_group(All_countries_data_males_ethiopia_reorder_clean_V1_noNA_all, "Hepatomegaly", "immstat_Male_data.dir/Ethiopia_HepatoV1_PLSDA/Ethiopia_Hepato_V1")
#Running the predcitions of the OR
#With Aspirate grade:
mixomics_logReg_scalling_2_new(All_countries_data_males_ethiopia_reorder_clean_V1_noNA_all, "Hepatomegaly", vector_traits_evalute2_AG, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/Ethiopia_HepatoV1_PLSDA/Eth_HepaV1_boot_alldata_AG")
#This runs the predictions with 6-8 and 10, 15 traits
run_splsda_several_values(All_countries_data_males_ethiopia_reorder_clean_V1_noNA_all, "Hepatomegaly", vector_traits_evalute2_AG, "Ethiopia","immstat_Male_data.dir/Ethiopia_HepatoV1_PLSDA/Eth_HepaV1_boot_alldata_AG")
#Logistic regression Age covariate:
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_ethiopia_reorder_clean_V1_noNA_all, "Hepatomegaly", vector_traits_evalute2_age, "Ethiopia", 1, 3, 100, "immstat_Male_data.dir/Ethiopia_HepatoV1_PLSDA/Eth_HepaV1_boot_alldata_AGE")


##############################
#Splenomegaly in Kenya on V2:#
##############################
selec_and_order_traits <- c( "UID", "Patient", "Class", "Country", "Males", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes",
                             "Age", "Height", "Weight",  "SSystolicBP", "SDiastolicBP", "Temperature",  "Pulse", "Hemoglobin",
                             "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                             "Aspirate_grade", "Spleen_size", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", 
                             "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")

All_countries_data_males_kenya_reorder <- All_countries_data_males_kenya[, selec_and_order_traits]
All_countries_data_males_kenya_Spl <- All_countries_data_males_kenya_reorder[
  !is.na(All_countries_data_males_kenya_reorder$Splenomegaly),]
#Removing HV and selecting only v2
All_countries_data_males_kenya_Spl <- 
  All_countries_data_males_kenya_Spl[which(All_countries_data_males_kenya_Spl$Class=="V2"),]
#Removing Weight and Spleen size, as there are a lot of missing data for V2 in it
All_countries_data_males_kenya_Spl <- 
  All_countries_data_males_kenya_Spl[,!colnames(All_countries_data_males_kenya_Spl) %in% c("Weight", "Spleen_size")]
#Removing all samples with missing data
All_countries_data_males_kenya_Spl_clean <- All_countries_data_males_kenya_Spl[complete.cases(All_countries_data_males_kenya_Spl), ]

#Plotting the difference between groups:
statsitical_compare_plot_group(All_countries_data_males_kenya_Spl_clean, "Splenomegaly", "immstat_Male_data.dir/Kenya_SplenoV2_PLSDA/KEN_Spl_V2")
#statsitical_compare_plot_group(All_countries_data_males_ethiopia_reorder_clean_V1_citokines, "Hepatomegaly", "immstat_Male_data.dir/Ethiopia_Hepato_V1_citokine")
#With Aspirate grade:
mixomics_logReg_scalling_2_new(All_countries_data_males_kenya_Spl_clean, "Splenomegaly", vector_traits_evalute2_AG, "Kenya", 1, 3, 1000, "immstat_Male_data.dir/Kenya_SplenoV2_PLSDA/KEN_SplV2_boot_alldata_AG")
run_splsda_several_values(All_countries_data_males_kenya_Spl_clean, "Splenomegaly", vector_traits_evalute2_AG, "Kenya","immstat_Male_data.dir/Kenya_SplenoV2_PLSDA/KEN_SplV2_boot_alldata_AG")
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_kenya_Spl_clean, "Splenomegaly", vector_traits_evalute2_age_AG, "Kenya", 1, 3, 1000, "immstat_Male_data.dir/Kenya_SplenoV2_PLSDA/KEN_SplV2_boot_alldata_Age_AG")


#################################
#Splenomegaly in Ethiopia on V2:#
#################################
selec_and_order_traits <- c( "UID", "Patient", "Class", "Country", "Males", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes",
                             "Age", "Height", "Weight",  "SSystolicBP", "SDiastolicBP", "Temperature",  "Pulse", "Hemoglobin",
                             "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                             "Aspirate_grade", "Spleen_size", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", 
                             "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")

All_countries_data_males_Ethiopia_reorder <- All_countries_data_males_ethiopia[, selec_and_order_traits]
All_countries_data_males_Ethiopia_Spl <- All_countries_data_males_Ethiopia_reorder[
  !is.na(All_countries_data_males_Ethiopia_reorder$Splenomegaly),]
#Removing HV and selecting only v2
All_countries_data_males_Ethiopia_Spl <- 
  All_countries_data_males_Ethiopia_Spl[which(All_countries_data_males_Ethiopia_Spl$Class=="V2"),]
#Removing all samples with missing data
All_countries_data_males_Ethiopia_Spl_clean <- All_countries_data_males_Ethiopia_Spl[complete.cases(All_countries_data_males_Ethiopia_Spl), ]
#Citokines
All_countries_data_males_Ethiopia_Spl_citokines <- All_countries_data_males_Ethiopia_Spl[,c(ids_to_keep,vector_traits_evalute)]
All_countries_data_males_Ethiopia_Spl_citokines <- All_countries_data_males_Ethiopia_Spl_citokines[
  complete.cases(All_countries_data_males_Ethiopia_Spl_citokines),]

#Statistical Comparison
statsitical_compare_plot_group(All_countries_data_males_Ethiopia_Spl_clean, "Splenomegaly", "immstat_Male_data.dir/Ethiopia_SplenoV2_PLSDA/Ethiopia_Spleno_V2_all")

#With Aspirate grade:
mixomics_logReg_scalling_2_new(All_countries_data_males_Ethiopia_Spl_clean, "Splenomegaly", vector_traits_evalute2_AG, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/Ethiopia_SplenoV2_PLSDA/ETH_SplV2_boot_AG")
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_Ethiopia_Spl_clean, "Splenomegaly", vector_traits_evalute2_age_AG, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/Ethiopia_SplenoV2_PLSDA/ETH_SplV2_boot_Age_AG")
run_splsda_several_values(All_countries_data_males_Ethiopia_Spl_clean, "Splenomegaly", vector_traits_evalute2_AG, "Ethiopia","immstat_Male_data.dir/Ethiopia_SplenoV2_PLSDA/ETH_SplV2_boot_AG")

#Logistic regression Age covariate:
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_Ethiopia_Spl_clean, "Splenomegaly", vector_traits_evalute2_age, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/Ethiopia_SplenoV2_PLSDA/ETH_SplV2_boot_Age")

################################################################
#Ethiopia using V1 to predict recovery from splenomegaly on V2:#
################################################################
#Selecting only v2, to get the splenomegaly data:
All_countries_data_males_Ethiopia_V2_Spl <- 
  All_countries_data_males_Ethiopia_reorder[which(All_countries_data_males_Ethiopia_reorder$Class=="V2"),]

#Getting the results from Splenomegaly from V2, renaming as V1 to integrate with V1 data
All_countries_data_males_Ethiopia_V2_Spl_subset <- All_countries_data_males_Ethiopia_V2_Spl[,c("UID","Splenomegaly" )]
colnames(All_countries_data_males_Ethiopia_V2_Spl_subset) <- c("UID","Splenomegaly_V2")
All_countries_data_males_Ethiopia_V2_Spl_subset$UID <- gsub("_V2_","_V1_", All_countries_data_males_Ethiopia_V2_Spl_subset$UID)
#Replacing V1 hepatomegaly for V2
All_countries_data_males_Ethiopia_V1_Spl <- 
  All_countries_data_males_Ethiopia_reorder[which(All_countries_data_males_Ethiopia_reorder$Class=="V1"),]
All_countries_data_males_ethiopia_Hep_V1_splenoV2 <- merge(All_countries_data_males_Ethiopia_V1_Spl, All_countries_data_males_Ethiopia_V2_Spl_subset, by="UID")
All_countries_data_males_ethiopia_Hep_V1_splenoV2_nona <- All_countries_data_males_ethiopia_Hep_V1_splenoV2[complete.cases(All_countries_data_males_ethiopia_Hep_V1_splenoV2),]

#Citokines
All_countries_data_males_ethiopia_Hep_V1_splenoV2citokines <- All_countries_data_males_ethiopia_Hep_V1_splenoV2[,c(ids_to_keep,vector_traits_evalute, "Splenomegaly_V2")]
All_countries_data_males_ethiopia_Hep_V1_splenoV2citokines <- All_countries_data_males_ethiopia_Hep_V1_splenoV2citokines[
  complete.cases(All_countries_data_males_ethiopia_Hep_V1_splenoV2citokines),]

#With Aspirate grade:
mixomics_logReg_scalling_2_new(All_countries_data_males_ethiopia_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_AG, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/Ethiopia_SplenoV1toV2_PLSDA/ETH_SplV1V2_boot_AG")
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_ethiopia_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_age_AG, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/Ethiopia_SplenoV1toV2_PLSDA/ETH_SplV1V2_boot_Age_AG")
run_splsda_several_values(All_countries_data_males_ethiopia_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_AG, "Ethiopia","immstat_Male_data.dir/Ethiopia_SplenoV1toV2_PLSDA/ETH_SplV1V2_boot_AG")

#Logistic regression Age covariate:
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_ethiopia_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_age, "Ethiopia", 1, 3, 1000, "immstat_Male_data.dir/ETH_SplV1V2_boot_Age")

#############################################################
#Kenya using V1 to predict recovery from splenomegaly on V2:#
#############################################################
All_countries_data_males_Kenya_V2_Spl <- 
  All_countries_data_males_kenya[which(All_countries_data_males_kenya$Class=="V2"),]

#Getting the results from Splenomegaly from V2, renaming as V1 to integrate with V1 data
All_countries_data_males_Kenya_V2_Spl_subset <- All_countries_data_males_Kenya_V2_Spl[,c("UID","Splenomegaly" )]
colnames(All_countries_data_males_Kenya_V2_Spl_subset) <- c("UID","Splenomegaly_V2")
All_countries_data_males_Kenya_V2_Spl_subset$UID <- gsub("_V2_","_V1_", All_countries_data_males_Kenya_V2_Spl_subset$UID)
#Replacing V1 hepatomegaly for V2
All_countries_data_males_Kenya_V1_Spl <- 
  All_countries_data_males_kenya[which(All_countries_data_males_kenya$Class=="V1"),]
All_countries_data_males_Kenya_Hep_V1_splenoV2 <- merge(All_countries_data_males_Kenya_V1_Spl, All_countries_data_males_Kenya_V2_Spl_subset, by="UID")
# nrow(All_countries_data_males_Kenya_Hep_V1_splenoV2)
# nrow(All_countries_data_males_Kenya_Hep_V1)

All_countries_data_males_Kenya_Hep_V1_splenoV2_nona <- All_countries_data_males_Kenya_Hep_V1_splenoV2[complete.cases(All_countries_data_males_Kenya_Hep_V1_splenoV2),]

#Citokines
All_countries_data_males_Kenya_Hep_V1_splenoV2citokines <- All_countries_data_males_Kenya_Hep_V1_splenoV2[,c(ids_to_keep,vector_traits_evalute, "Splenomegaly_V2")]
All_countries_data_males_Kenya_Hep_V1_splenoV2citokines <- All_countries_data_males_Kenya_Hep_V1_splenoV2citokines[
  complete.cases(All_countries_data_males_Kenya_Hep_V1_splenoV2citokines),]


#Statistical Comparison
statsitical_compare_plot_group(All_countries_data_males_Kenya_Hep_V1_splenoV2_nona, "Splenomegaly_V2", "immstat_Male_data.dir/Kenya_SplenoV1toV2_PLSDA/KEN_SplV1V2_boot")

#With Aspirate grade:
mixomics_logReg_scalling_2_new(All_countries_data_males_Kenya_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_AG, "Kenya", 1, 3, 1000, "immstat_Male_data.dir/Kenya_SplenoV1toV2_PLSDA/KEN_SplV1V2_boot_AG")
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_Kenya_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_age_AG, "Kenya", 1, 3, 1000, "immstat_Male_data.dir/Kenya_SplenoV1toV2_PLSDA/KEN_SplV1V2_boot_Age_AG")
run_splsda_several_values(All_countries_data_males_Kenya_Hep_V1_splenoV2_nona, "Splenomegaly_V2", vector_traits_evalute2_AG, "Kenya","immstat_Male_data.dir/Kenya_SplenoV1toV2_PLSDA/KEN_SplV1V2_boot_AG")



###########################
#Sudan hepatomegaly on V1:#
###########################
selec_and_order_traits <- c( "UID", "Patient", "Class", "Country", "Males", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes",
                             "Age", "Height", "Weight",  "SSystolicBP", "SDiastolicBP", "Temperature",  "Pulse", "Hemoglobin",
                             "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                             "Aspirate_grade", "Spleen_size", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", 
                             "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")

All_countries_data_males_sudan_reorder <- All_countries_data_males_sudan[, selec_and_order_traits]
All_countries_data_males_sudan_reorder_clean <- All_countries_data_males_sudan_reorder[which(All_countries_data_males_sudan_reorder$Class != "HV"),]
All_countries_data_males_sudan_reorder_clean <- All_countries_data_males_sudan_reorder_clean[,!colnames(All_countries_data_males_sudan_reorder_clean) %in%
                                                                                                     c("Auxiliar_Lymphnodes", "Aspirate_grade", "Spleen_size")]
#Selecting V1
All_countries_data_males_sudan_reorder_clean_V1 <-
  All_countries_data_males_sudan_reorder_clean[All_countries_data_males_sudan_reorder_clean$Class=="V1",]
#Removing missing data
All_countries_data_males_sudan_reorder_clean_V1_noNA <- All_countries_data_males_sudan_reorder_clean_V1[
  complete.cases(All_countries_data_males_sudan_reorder_clean_V1),]

#Citokines
All_countries_data_males_sudan_reorder_clean_V1citokines <- All_countries_data_males_sudan_reorder_clean_V1[,c(ids_to_keep,vector_traits_evalute)]
All_countries_data_males_sudan_reorder_clean_V1citokines <- All_countries_data_males_sudan_reorder_clean_V1citokines[
  complete.cases(All_countries_data_males_sudan_reorder_clean_V1citokines),]

#Statistical Comparison
statsitical_compare_plot_group(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Hepatomegaly", "immstat_Male_data.dir/Sudan_hepV1_PLSDA/SUD_HepV1_boot")
mixomics_logReg_scalling_2_new(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Hepatomegaly", vector_traits_evalute2, "Sudan", 1, 3, 1000, "immstat_Male_data.dir/Sudan_hepV1_PLSDA/SUD_HepV1_boot")
#This runs the predictions with 6-8 and 10, 15 traits
run_splsda_several_values(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Hepatomegaly", vector_traits_evalute2, "Sudan","immstat_Male_data.dir/Sudan_hepV1_PLSDA/SUD_HepV1_boot")
#Logistic regression Age covariate:
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Hepatomegaly", vector_traits_evalute2_age, "Sudan", 1, 3, 1000, "immstat_Male_data.dir/Sudan_hepV1_PLSDA/SUD_HepV1_boot_Age")
#There are no aspirate grade in Sudan


###########################
#Sudan Splenomegaly on V1:#
###########################
mixomics_logReg_scalling_2_new(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Splenomegaly", vector_traits_evalute2, "Sudan", 1, 3, 1000, "immstat_Male_data.dir/Sudan_SplenoV1_PLSDA/SUD_SplV1_boot")
#This runs the predictions with 6-8 and 10, 15 traits
run_splsda_several_values(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Splenomegaly", vector_traits_evalute2, "Sudan","immstat_Male_data.dir/Sudan_SplenoV1_PLSDA/SUD_SplV1_boot")
#Statistical Comparison
statsitical_compare_plot_group(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Splenomegaly", "immstat_Male_data.dir/SUD_SplV1_boot")
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_sudan_reorder_clean_V1_noNA, "Splenomegaly", vector_traits_evalute2_age, "Sudan", 1, 3, 1000, "immstat_Male_data.dir/Sudan_SplenoV1_PLSDA/SUD_SplV1_boot_Age")
#There are no aspirate grade in Sudan

############################
#Uganda Splenomegaly on V1:#
############################
selec_and_order_traits <- c( "UID", "Patient", "Class", "Country", "Males", "Hepatomegaly", "Splenomegaly", "Auxiliar_Lymphnodes",
                             "Age", "Height", "Weight",  "SSystolicBP", "SDiastolicBP", "Temperature",  "Pulse", "Hemoglobin",
                             "WBCells", "Neutrophil", "Lymphocyte", "Platelets", "ALT_Results", "Creatinine", "Albumin", "Total_Bilirubin",
                             "Aspirate_grade", "Spleen_size", "CX3CL1", "CXCL12", "PTX3", "TGF.B1", "sCD25", "sCD40L", 
                             "sRAGE",  "sST2", "sTNF.RI", "sTNF.RII",  "sTREM.1")

All_countries_data_males_uganda_reorder <- All_countries_data_males_uganda[, selec_and_order_traits]
All_countries_data_males_uganda_reorder_clean <- All_countries_data_males_uganda_reorder[which(All_countries_data_males_uganda_reorder$Class != "HV"),]
All_countries_data_males_uganda_reorder_clean <- All_countries_data_males_uganda_reorder_clean[,!colnames(All_countries_data_males_uganda_reorder_clean) %in%
                                                                                               c("Auxiliar_Lymphnodes", "Spleen_size")]
#Selecting V1
All_countries_data_males_uganda_reorder_clean_V1 <-
  All_countries_data_males_uganda_reorder_clean[All_countries_data_males_uganda_reorder_clean$Class=="V1",]
#Removing missing data
All_countries_data_males_uganda_reorder_clean_V1_noNA <- All_countries_data_males_uganda_reorder_clean_V1[
  complete.cases(All_countries_data_males_uganda_reorder_clean_V1),]

#Citokines
All_countries_data_males_uganda_reorder_cleancitokines <- All_countries_data_males_uganda_reorder_clean[,c(ids_to_keep,vector_traits_evalute)]
All_countries_data_males_uganda_reorder_cleancitokines <- All_countries_data_males_uganda_reorder_cleancitokines[
  complete.cases(All_countries_data_males_uganda_reorder_cleancitokines),]

statsitical_compare_plot_group(All_countries_data_males_uganda_reorder_clean_V1_noNA, "Splenomegaly", "immstat_Male_data.dir/UGA_SplV1_boot")
#With aspirate grade:
mixomics_logReg_scalling_2_new(All_countries_data_males_uganda_reorder_clean_V1_noNA, "Splenomegaly", vector_traits_evalute2_AG, "Uganda", 1, 3, 1000, "immstat_Male_data.dir/Uganda_splenV1_PLSDA/UGA_SplV1_boot_AG")
mixomics_logReg_scalling_2_new_age_covariate(All_countries_data_males_uganda_reorder_clean_V1_noNA, "Splenomegaly", vector_traits_evalute2_age_AG, "Uganda", 1, 3, 1000, "immstat_Male_data.dir/Uganda_splenV1_PLSDA/UGA_SplV1_boot_Age_AG")
run_splsda_several_values(All_countries_data_males_uganda_reorder_clean_V1_noNA, "Splenomegaly", vector_traits_evalute2_AG, "Uganda","immstat_Male_data.dir/Uganda_splenV1_PLSDA/UGA_SplV1_boot_AG")


#######################################################################################
save.image(file = "Current_results_all.RData")
#load("Current_results_all.RData")










