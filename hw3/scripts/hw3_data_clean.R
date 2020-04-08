# hw3_data_clean.R
# Noah Stafford

# libraries
library(tidyverse)
library(readxl)

# Setwd where RawData files are stored
setwd("C:/Users/woody/Google Drive/github_repos/826_Advanced_Data_Analysis/hw3/RawData")

### Save data into accessible csv files from their current state in excel files

### final_rbm_data.xls
# RBM immunoassay data
final_rbm_data <- read_excel("final_rbm_data.xls", sheet = "Antigen Immunoassay", skip = 16, col_names = FALSE)

# Remove rows with non-data notes at the end
final_rbm_data <- final_rbm_data[1:553,]

# Read in the entire file and extract the column headers as a vector and assign the column names
raw_rbm <- read_excel("final_rbm_data.xls", sheet = "Antigen Immunoassay")
rbm_columns <- as.character(raw_rbm[11,])
rbm_columns[1] <- "samples"
colnames(final_rbm_data) <- rbm_columns
# remove empty columns
final_rbm_data <- final_rbm_data[, !(is.na(colnames(final_rbm_data)))]
colnames(final_rbm_data)[1] <- "sample_id"

# Data cleaning

# Print problem columns
for (c in 1:ncol(final_rbm_data)){
  print(colnames(final_rbm_data[,c]))
  class_vec <- sapply(final_rbm_data[,c],class)
  if(class_vec == "character"){
  print(table(final_rbm_data[,c]))
  }
}

# Change all '<LOW>' values to NA
final_rbm_data[final_rbm_data == '<LOW>'] <- NA

# Check if all rows have numeric data in each character column, and recast if they do
cast_to_numeric <- c()
for (c in 1:ncol(final_rbm_data)){
  print(colnames(final_rbm_data[,c]))
  class_vec <- sapply(final_rbm_data[,c],class)
  
  if(class_vec == "character"){
    #print(sapply(final_rbm_data[, c], str_remove, "[.]"))
    bin_vec <- !grepl("\\D", sapply(final_rbm_data[, c], str_remove, "[.]"))
    #print(bin_vec)
    n_numeric <- sum(bin_vec)
    if(n_numeric == nrow(final_rbm_data)){
      cast_to_numeric <- c(cast_to_numeric,colnames(final_rbm_data[,c]))
    }
  }
}

# Re-cast columns in cast_to_numeric
final_rbm_data<- cbind(final_rbm_data %>% select(cast_to_numeric) %>% mutate_if(is.character, as.numeric),
      final_rbm_data %>% select(colnames(final_rbm_data)[!colnames(final_rbm_data) %in% cast_to_numeric]))

# Fix Fibrinogen
final_rbm_data$Fibrinogen <- sapply(final_rbm_data$Fibrinogen,str_remove, "[>]")
final_rbm_data$Fibrinogen <- as.numeric(final_rbm_data$Fibrinogen)

# Fix Clusterin 
final_rbm_data$Clusterin <- sapply(final_rbm_data$Clusterin,str_remove, "[>]")
final_rbm_data$Clusterin <- as.numeric(final_rbm_data$Clusterin)

# Fix `Growth Hormone` 
final_rbm_data$`Growth Hormone` <- sapply(final_rbm_data$`Growth Hormone`,str_remove, "[>]")
final_rbm_data$`Growth Hormone` <- as.numeric(final_rbm_data$`Growth Hormone`)


### RBM_Tube_Number_Key.xls
# contains the mouse/sample ID correspondence
# for final_rbm_data.xls
rbm_tube <- read_excel("RBM_Tube_Number_Key.xls", sheet = "Sheet1", col_names = TRUE)
colnames(rbm_tube) <- c("mouse_id", "sample_id")

# Remove whitespace from mouse_id
rbm_tube$mouse_id <- gsub(" ", "", rbm_tube$mouse_id, fixed = TRUE)

# Extract number from mouse_id and cast as a number
rbm_tube$mouse_id <- as.numeric(str_split_fixed(rbm_tube$mouse_id, "#", 2)[,2])

### CPL_Rosetta_Lipids_FINAL.xls
# The trait measurements are in the last four columns
cpl_rosetta <- read_excel("CPL_Rosetta_Lipids_FINAL.xls", sheet = "BTBRobxB60B", col_names = TRUE)

# rermove rows with non-data notes at the end
cpl_rosetta <- cpl_rosetta[1:554, c(1,7:10)]

# create mouse_id
split_matrix <- str_split_fixed(cpl_rosetta$Mouse, "_", 4)
mouse_id_vec <- split_matrix[,2]
cpl_rosetta$mouse_id <- str_replace(mouse_id_vec, "Mouse", "Mouse#")
cpl_rosetta$mouse_id <- as.numeric(str_split_fixed(cpl_rosetta$mouse_id, "#", 2)[,2])

# Clean data
to_na_vec <- c('not received', 'QNS', 'less than 2.00')
cpl_rosetta[cpl_rosetta %in% to_na_vec] <- NA

cpl_rosetta$`NEFA mEq/L` <- as.numeric(cpl_rosetta$`NEFA mEq/L`)
cpl_rosetta$`LDL mg/dL` <- as.numeric(cpl_rosetta$`LDL mg/dL`)
cpl_rosetta$`HDL mg/dL` <- as.numeric(cpl_rosetta$`HDL mg/dL`)
cpl_rosetta$`Total CHOL. mg/dL` <- as.numeric(cpl_rosetta$`Total CHOL. mg/dL`)


### Complete F2 Liver TG Set.xls
# Relevant are Sex, Body Weight, Total.Liver.Weight,
# Glucose, Insulin, Plasma TG, ug TG/mg Protein
complete_f2 <- read_excel("Complete F2 Liver TG Set.xls", sheet = "TG", skip = 1)
colnames(complete_f2)[1] <- "mouse_id"
complete_f2  <- complete_f2 %>% select(mouse_id, Sex, `Body Weight`, Total.Liver.Weight.g., `10 wk Glucose mg/gl`, `10wk Insulin ng/ml`,
                                       `10wk Insulin ng/ml`, `10wk Plasma TG`, `ug TG/mg Protein`)
complete_f2$mouse_id <- gsub(" ", "", complete_f2$mouse_id, fixed = TRUE)
complete_f2$mouse_id <- as.numeric(str_split_fixed(complete_f2$mouse_id, "#", 2)[,2])

# Data cleaning
table(complete_f2$Total.Liver.Weight.g.)
complete_f2$Total.Liver.Weight.g. <- as.numeric(complete_f2$Total.Liver.Weight.g.)

### pheno_lipomics_bleeds.xls
# we're interested in sex,
# weeks 4,6,8,10: weight, length, glucose, insulin, triglyceride
# last few columns with week 10 insulin, c-peptide, and glucose

# Generate column names from two-row nested configuration of labeling in .xls file
pheno_lipomics_colnames <- read_excel("pheno_lipomics_bleeds.xls", col_names = FALSE)
pheno_header_row <- as.character(pheno_lipomics_colnames[1,])

lag_element <- ""
for (i in 1:length(pheno_header_row)){
  element <- pheno_header_row[i]
  
  if(is.na(element)){
    pheno_header_row[i] <- lag_element
  } else {
    lag_element <- element
  }
}

pheno_column_row <- as.character(pheno_lipomics_colnames[2,])
pheno_column_names <- paste(pheno_header_row, pheno_column_row)
pheno_column_names <- trimws(pheno_column_names)

pheno_lipomics <- read_excel("pheno_lipomics_bleeds.xls", col_names = FALSE, skip = 2)
colnames(pheno_lipomics) <- pheno_column_names

pheno_lipomics <- pheno_lipomics %>% select(`Mouse ID`, SEX, 
                                            `4 wk Orbital Eye Bleed WEIGHT (g)`,
                                            `4 wk Orbital Eye Bleed BODY LENGTH (cm)`,
                                            `4 wk Orbital Eye Bleed GLUCOSE (mg/dl)`, 
                                            `4 wk Orbital Eye Bleed INSULIN (ng/ml)`,
                                            `4 wk Orbital Eye Bleed TRIGLYCERIDE (mg/dl)`, 
                                            `6 wk Orbital Eye Bleed WEIGHT (g)`,
                                            `6 wk Orbital Eye Bleed BODY LENGTH (cm)`,
                                            `6 wk Orbital Eye Bleed GLUCOSE (mg/dl)`, 
                                            `6 wk Orbital Eye Bleed INSULIN (ng/ml)`,
                                            `6 wk Orbital Eye Bleed TRIGLYCERIDE (mg/dl)`, 
                                            `8 wk Orbital Eye Bleed WEIGHT (g)`,
                                            `8 wk Orbital Eye Bleed BODY LENGTH (cm)`,
                                            `8 wk Orbital Eye Bleed GLUCOSE (mg/dl)`, 
                                            `8 wk Orbital Eye Bleed INSULIN (ng/ml)`,
                                            `8 wk Orbital Eye Bleed TRIGLYCERIDE (mg/dl)`, 
                                            `10 wk Orbital Eye Bleed WEIGHT (g)`,
                                            `10 wk Orbital Eye Bleed BODY LENGTH (cm)`,
                                            `10 wk Orbital Eye Bleed GLUCOSE (mg/dl)`, 
                                            `10 wk Orbital Eye Bleed INSULIN (ng/ml)`,
                                            `10 wk Orbital Eye Bleed TRIGLYCERIDE (mg/dl)`,
                                            `TNB- Cholesterol 10 wk c peptide`,
                                            `TNB- Cholesterol 10 wk insulin`,
                                            `TNB- Cholesterol 10 wk glucose triplicate repeat`)
pheno_lipomics$mouse_id <- paste0("Mouse#", pheno_lipomics$`Mouse ID`)
pheno_lipomics$mouse_id <- as.numeric(str_split_fixed(pheno_lipomics$mouse_id, "#", 2)[,2])
pheno_lipomics$mouse_id <- ifelse(pheno_lipomics$mouse_id == 9.5, 3611, pheno_lipomics$mouse_id)

# Data cleaning
# Print problem columns
for (c in 1:ncol(pheno_lipomics)){
  print(colnames(pheno_lipomics[,c]))
  class_vec <- sapply(pheno_lipomics[,c],class)
  if(class_vec == "character"){
    print(table(pheno_lipomics[,c]))
  }
}

pheno_lipomics$`4 wk Orbital Eye Bleed GLUCOSE (mg/dl)` <- as.numeric(pheno_lipomics$`4 wk Orbital Eye Bleed GLUCOSE (mg/dl)`)
pheno_lipomics$`6 wk Orbital Eye Bleed WEIGHT (g)` <- as.numeric(pheno_lipomics$`6 wk Orbital Eye Bleed WEIGHT (g)`)
pheno_lipomics$`8 wk Orbital Eye Bleed WEIGHT (g)` <- as.numeric(pheno_lipomics$`8 wk Orbital Eye Bleed WEIGHT (g)`)
pheno_lipomics$`10 wk Orbital Eye Bleed BODY LENGTH (cm)` <- as.numeric(pheno_lipomics$`10 wk Orbital Eye Bleed BODY LENGTH (cm)`)

### necropsy_tracking_report.xls  
# we're interested in the organ weights
# plus number of harvested pancreatic islets ("Approx. # harvested")
necropsy_tracking <- read_excel("necropsy_tracking_report.xls", sheet = "Necropsy_Tracking_Report", skip = 1)
necropsy_tracking <- necropsy_tracking %>% select(`Mouse ID`, `Hypothalamus weight(mg)`, `Brain weight(mg)`, `Liver weight(mg)`, `Rt. Kidney weight(mg)`,
                             `Rt. Kidney weight(mg)`, `Lt. Adipose weight(g)`, `Total Liver Weight(g)`, `Remaining Liver Weight(g)`,
                             `Soleus Weight(mg)`, `Gastroc. Weight(mg)`, `Spleen weight(mg)`, `Heart weight(mg)`, `Sacrifice Weight`, 
                             `Weight(mg)`, `Approx. # harvested`)
necropsy_tracking$mouse_id <- gsub(" ", "", necropsy_tracking$`Mouse ID`, fixed = TRUE)
necropsy_tracking$mouse_id <- as.numeric(str_split_fixed(necropsy_tracking$mouse_id, "#", 2)[,2])

# Clean Data
summary(necropsy_tracking)
necropsy_tracking$`Hypothalamus weight(mg)` <- as.numeric(necropsy_tracking$`Hypothalamus weight(mg)`)
necropsy_tracking$`Rt. Kidney weight(mg)` <- as.numeric(necropsy_tracking$`Rt. Kidney weight(mg)`)
necropsy_tracking$`Lt. Adipose weight(g)` <- as.numeric(necropsy_tracking$`Lt. Adipose weight(g)`)
necropsy_tracking$`Total Liver Weight(g)` <- as.numeric(necropsy_tracking$`Total Liver Weight(g)`)
necropsy_tracking$`Remaining Liver Weight(g)` <- as.numeric(necropsy_tracking$`Remaining Liver Weight(g)`)
necropsy_tracking$`Heart weight(mg)` <- as.numeric(necropsy_tracking$`Heart weight(mg)`)


### Merge datasets
immunoassay <- merge(final_rbm_data, rbm_tube, by = "sample_id", all.x=TRUE)

# Merge immunoassay data with cpl_rosetta
immunoassay <- merge(immunoassay, cpl_rosetta, by = "mouse_id", all.x=TRUE)

# Merge immunoassay data with complete_f2
immunoassay <- merge(immunoassay, complete_f2, by = "mouse_id", all.x=TRUE)

# Merge immunoassay data with pheno_lipomics
immunoassay <- merge(immunoassay, pheno_lipomics, by = "mouse_id", all.x=TRUE)

# Merge immunoassay data with necropsy_tracking
immunoassay <- merge(immunoassay, necropsy_tracking, by = "mouse_id", all.x = TRUE)

# Rewrite some column names for clarity
colnames(immunoassay)[77] <- "sex_f2" 
colnames(immunoassay)[84] <- "Mouse ID_pl"
colnames(immunoassay)[85] <- "sex_pl" 
colnames(immunoassay)[109] <- "Mouse ID_nt" 

# Write cleaned data to file
save(immunoassay, file = "../cleanData/immunoassay.RData")
