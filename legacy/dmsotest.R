###########################Section 1###############################
## Cleaning commands and graphs & Setting working directory
cat("\014")
rm(list=ls())
setwd("D:/Masters_project/data/filter")
getwd()

#load packages
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyr)
suppressPackageStartupMessages(library(plotly))
library(ggplot2)
library(ggdendro)
library(tibble)
suppressPackageStartupMessages(library(permute))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(DescTools))
library(plot3D)
library(berryFunctions)
library(plotrix)
library(gridExtra)
library(MASS)
library(shape)


# function to remove the temperature
remove.temp <- function(x){
  # remove the temperature column
  y <- x[,c(1,3:98)]
  return(y)
}

# function to make the time into hours
clean.time <- function(x){
  x$Time <- sapply(strsplit(x$Time, ":"),
                   function(z) {
                     z <- as.numeric(z)
                     z[1]+(z[2]+z[3]/60)/60
                   })
  return(x) 
}

## Defining function to calculate auc
AreaUnderSpline <- function(Time, OD600, minT, maxT){
  spline.area <- stats::integrate(stats::splinefun(Time, OD600, method = "natural"), 
                                  lower = minT, upper = maxT)$value
  return(spline.area)
}


###########################Section 2####################################
### read data of all plates and calculate auc for each well
# Setting strain, chemical library, and replica No.
strain <- c("85","398","436")
library <- c("1", "2")
rep <- c("1", "2", "3")
experiment <- c("_1", "_2")

# creating data frame for auc values
all_data <- data.frame()
columns <- c("Strain", "Library", "Replica", "Well", "AUC", "Exp")
auc_values <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(auc_values) <- columns

# read data and calculate auc data for all wells
for(i in strain){
  for(j in library){
    for(k in rep){
      for(ex in experiment){
        filepath <- paste(i, "_p_", j, "_r_", k, ex, ".txt", sep = "")
        print(filepath) 
        
        # read in the data
        plate_data <-  clean.time(remove.temp(read.table(filepath,  sep = '\t', header = TRUE, check.names = FALSE)))
        # rearrange the data a bit first
        # make the data from wide to long format
        clean_data <- pivot_longer(plate_data, cols = c(2:97), names_to = "Well", values_to = "OD600") %>%
          arrange(Well) 
        
        # add metadata
        clean_data$Strain <- i # add info like strain, replicate, ps1/ps2
        clean_data$Library <- j
        clean_data$Replicate <- k
        clean_data$Exp <- ex
        
        # calculating auc values for all wells in that plate
        for (l in unique(clean_data$Well)){
          
          ## Get data for each well
          data_subset <- clean_data %>% filter(Well == l)
          Time <- data_subset$Time 
          OD600 <- data_subset$OD600
          
          # calculating auc values
          m <- AreaUnderSpline(data_subset$Time, data_subset$OD600, 1, 73)
          print(m)
          # giving auc of all wells from all plates (replica 1,2,3 & ps1, ps2) by adding values from loops
          auc_values[nrow(auc_values) + 1,] <- c(i,j,k,l,m,ex)
        }
        
        # giving growth curve data of all wells from all plates during 72 hrs by adding values from loops
        all_data <- bind_rows(all_data, clean_data)
      }
    }
  }
}


# read chemical data
ps1 <- read.csv("ps1.csv", header = TRUE)
ps2 <- read.csv("ps2.csv", header = TRUE)

# change the name from x1, x2... to 1, 2....
colnames(ps1) <- c("Row", seq(1, 12))
colnames(ps2) <- c("Row", seq(1, 12))

# Transform the data
ps1_long <- pivot_longer(ps1, cols = c(2:13), names_to = "Column", values_to = "Chemical") %>%
  arrange(Row)
ps2_long <- pivot_longer(ps2, cols = c(2:13), names_to = "Column", values_to = "Chemical") %>%
  arrange(Row)

# combine column 'row' and 'column' to 'Well'
ps1_long$Well <- paste(ps1_long$Row, ps1_long$Column, sep = "")
ps1_long = subset(ps1_long, select = -c(Row, Column) ) %>%
  arrange(Well)

ps2_long$Well <- paste(ps2_long$Row, ps2_long$Column, sep = "")
ps2_long = subset(ps2_long, select = -c(Row, Column) ) %>%
  arrange(Well)

# bind chemical data with corresponding plates & wells
auc_values_all <- data.frame()
for (i in unique(auc_values$Library)) {
  library_ps <- paste("ps",i,"_long",sep = "")
  auc_values_library <- auc_values[auc_values$Library == i,] %>%
    left_join(get(library_ps), by = "Well")
  auc_values_all <- bind_rows(auc_values_all, auc_values_library)
} 

#############################Section 3.5###################################
# correct the wells of plates of replica 2 and 3
# create data frame of well arrangements
row_alphabet <- LETTERS[1:8]
column_number <- c(seq(1,12))
well_arr <- expand.grid(row_alphabet, column_number)
colnames(well_arr) <- c("row","col")
well_arr$well <- paste(well_arr$row, well_arr$col, sep = ":")
well_arr = subset(well_arr) %>% 
  arrange(well)

# define functions to correct well arrangements
rep2_func_row <- function(data){
  data$well_row2 <- sapply(strsplit(data$well, ":"),
                           function(z){
                             if (which(z[1] == row_alphabet) > 6) {
                               new_row1 <- row_alphabet[which(z[1] == row_alphabet) - 6]
                               print(new_row1)
                             }
                             else{
                               new_row2 <- row_alphabet[which(z[1] == row_alphabet) + 2]
                               print(new_row2)
                             }
                           })
  return(data)
}

rep2_func_col <- function(data){
  data$well_col2 <- sapply(strsplit(data$well, ":"),
                           function(z){
                             if (which(as.numeric(z[2]) == column_number) > 10) {
                               new_col1 <- column_number[which(as.numeric(z[2]) == column_number) - 10]
                               print(new_col1)
                             }
                             else{
                               new_col2 <- column_number[which(as.numeric(z[2]) == column_number) + 2]
                               print(new_col2)
                             }
                           })
  return(data)
}

rep3_func_row <- function(data){
  data$well_row3 <- sapply(strsplit(data$well, ":"),
                           function(z){
                             if (which(z[1] == row_alphabet) > 4) {
                               new_row1 <- row_alphabet[which(z[1] == row_alphabet) - 4]
                               print(new_row1)
                             }
                             else{
                               new_row2 <- row_alphabet[which(z[1] == row_alphabet) + 4]
                               print(new_row2)
                             }
                           })
  return(data)
}

rep3_func_col <- function(data){
  data$well_col3 <- sapply(strsplit(data$well, ":"),
                           function(z){
                             if (which(as.numeric(z[2]) == column_number) > 8) {
                               new_col1 <- column_number[which(as.numeric(z[2]) == column_number) - 8]
                               print(new_col1)
                             }
                             else{
                               new_col2 <- column_number[which(as.numeric(z[2]) == column_number) + 4]
                               print(new_col2)
                             }
                           })
  return(data)
}

# get correct arrangements for replica 2 and 3, and paste them together.
well_arr <- rep2_func_row(well_arr)
well_arr <- rep2_func_col(well_arr)
well_arr <- rep3_func_row(well_arr)
well_arr <- rep3_func_col(well_arr)
well_arr$well_rep2 <- paste(well_arr$well_row2, well_arr$well_col2, sep = "")
well_arr$well_rep3 <- paste(well_arr$well_row3, well_arr$well_col3, sep = "")

# make the data with all auc values into correct wells
for (i in unique(auc_values_all$Replica)) {
  if (i == 1) {
    next
  }
  else if (i == 2) {
    ps1_r2 <- ps1_long
    ps2_r2 <- ps2_long
    ps1_r2$Well <- well_arr$well_rep2
    ps2_r2$Well <- well_arr$well_rep2
    ps1_r2 <- arrange(ps1_r2, Well)
    ps2_r2 <- arrange(ps2_r2, Well)
    auc_values_all$Chemical[auc_values_all$Replica == i & auc_values_all$Library == 1] <- ps1_r2$Chemical
    auc_values_all$Chemical[auc_values_all$Replica == i & auc_values_all$Library == 2] <- ps2_r2$Chemical
  }
  else {
    ps1_r3 <- ps1_long
    ps2_r3 <- ps2_long
    ps1_r3$Well <- well_arr$well_rep3
    ps2_r3$Well <- well_arr$well_rep3
    ps1_r3 <- arrange(ps1_r3, Well)
    ps2_r3 <- arrange(ps2_r3, Well)
    auc_values_all$Chemical[auc_values_all$Replica == i & auc_values_all$Library == 1] <- ps1_r3$Chemical
    auc_values_all$Chemical[auc_values_all$Replica == i & auc_values_all$Library == 2] <- ps2_r3$Chemical
  }
}

auc_values_all <- arrange(auc_values_all, Strain)


control_auc <- data.frame()
ratio <- data.frame()
ratio_strain <- data.frame()
for (i in unique(auc_values_all$Strain)) {
  DMSO <- mean(as.numeric(auc_values_all[auc_values_all$Chemical == "DMSO" &
                                           auc_values_all$Strain == i,]$AUC)) %>%
    as.data.frame()
  colnames(DMSO) <- "control_auc"
  DMSO$Strain <- i
  control_auc <- bind_rows(control_auc, DMSO)
  
  ratio_strain <- as.numeric(auc_values_all$AUC[auc_values_all$Strain == i])/control_auc$control_auc[control_auc$Strain == i]
  ratio_strain <- as.data.frame(ratio_strain)
  ratio_strain$Strain <- i
  ratio <- bind_rows(ratio, ratio_strain)
} 

# bonding ratio to plate data
auc_values_all$AUC_ratio <- ratio$ratio_strain

# delete data of wells with empty chemical
auc_values_all <- filter(auc_values_all, Chemical != "empty")


auc_values_all$ID <- paste("PS:", 
                           auc_values_all$Library, 
                           "_Rep:", 
                           auc_values_all$Replica, 
                           "_Well:", 
                           auc_values_all$Well, 
                           "_Exp:", 
                           auc_values_all$Exp, 
                           sep =  "")

all_data$ID <- paste("PS:", 
                     all_data$Library, 
                     "_Rep:", 
                     all_data$Replicate, 
                     "_Well:", 
                     all_data$Well, 
                     "_Exp:", 
                     all_data$Exp, 
                     sep =  "")

IDC <- data.frame()
IDC <- auc_values_all[auc_values_all$Chemical == "DMSO" & auc_values_all$Strain == "85",]
all_data_dmso <- all_data[all_data$ID %in% IDC$ID,]
auc_values_all_dmso <- auc_values_all[auc_values_all$Chemical == "DMSO",]



ggplot(all_data %>% filter(Strain == "85" & Replicate == "1" & Library == "2" & Exp == "_2"), aes(x = Time, y = OD600, col = Well)) +
  geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Well, nrow = 8)



ggplot(all_data_dmso %>% filter(Strain == "436" & Exp == "_2" & Library == "2" & Replicate == "3"), aes(x = Time, y = OD600, col = ID)) +
  geom_point() +
  facet_wrap(~Well, nrow = 8)
