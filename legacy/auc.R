####################Section 1####################################
## Cleaning commands and graphs & Setting working directory
cat("\014")
rm(list=ls())
if(!is.null(dev.list())) dev.off
setwd("D:/Masters_project/data/growth-curves")
getwd()

#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

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
#########################Section 2#############################

























####################################################################
## Change this to point to one of your files
filepath <- "333_p_1_r_1.txt"

# read in the data
plate_data <-  clean.time(remove.temp(read.table(filepath,  sep = '\t', header = TRUE, check.names = FALSE)))

# make the data from wide to long format
clean_data <- pivot_longer(plate_data, cols = c(2:97), names_to = "Well", values_to = "OD600") %>%
  arrange(Well) 

# try plotting the data
# for each well (LOTS OF GRAPHS!!! CAUTION!)
for (welln in unique(clean_data$Well)){
  dev.new()
  print( ggplot(clean_data[clean_data$Well == welln,], aes(x = Time, y = OD600, col = Well)) +
      geom_point() +
      theme(legend.position = "none") )
}

# All together
ggplot(clean_data, aes(x = Time, y = OD600, col = Well)) +
  geom_point() +
  theme(legend.position = "none")

# All together
ggplot(clean_data, aes(x = Time, y = OD600, col = Well)) +
  geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Well, nrow = 8)



###########################################################
# Area under the curve

# subset data to the first well
data_subset <- clean_data %>% filter(Well == "A1")

# pull out the time and optical densities for the curve
Time <- data_subset$Time 
OD600 <- data_subset$OD600

# spline function
spline.curve <- stats::splinefun(Time, OD600, method = "natural")
# calculate the area under the spline curve
spline.area <- (stats::integrate(spline.curve, lower = 1, upper = 73))$value

## Function to calculate AUC all in one go

AreaUnderSpline <- function(Time, OD600, minT, maxT){
  spline.area <- stats::integrate(stats::splinefun(Time, OD600, method = "natural"), 
                                  lower = minT, upper = maxT)$value
  return(spline.area)
}



###############################Section 3######################################
columns <- c("Well", "AUC", "strain", "library", "replica")
auc_values <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(auc_values) <- columns

#Calculate auc for each well
for (l in unique(all_data$Well)){
  
  ## Get data for each well
  data_subset <- all_data %>% filter(Well == l)
  Time <- data_subset$Time 
  OD600 <- data_subset$OD600
  
  ## Defining function to calculate auc
  AreaUnderSpline <- function(Time, OD600, minT, maxT){
    spline.area <- stats::integrate(stats::splinefun(Time, OD600, method = "natural"), 
                                    lower = minT, upper = maxT)$value
    return(spline.area)
  }
  
  m <- AreaUnderSpline(data_subset$Time, data_subset$OD600, 1, 73)
  print(m)
  auc_values[nrow(auc_values) + 1,] <- c(l,m)
  
}


##############################################################################
# read data of all plates and calculate auc for each well
strain <- c("100", "186")
library <- c("1", "2")
rep <- c("1", "2", "3")

  all_data <- data.frame()
  columns <- c("Strain", "Library", "Replica", "Well", "AUC")
  auc_values <- data.frame(matrix(nrow = 0, ncol = length(columns)))
    colnames(auc_values) <- columns

for(i in strain){
  for(j in library){
    for(k in rep){
      filepath <- paste(i, "_p_", j, "_r_", k, ".txt", sep = "")
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
      
      # calculating auc values for all wells in that plate
      for (l in unique(clean_data$Well)){
        
        ## Get data for each well
        data_subset <- clean_data %>% filter(Well == l)
        Time <- data_subset$Time 
        OD600 <- data_subset$OD600
        
        ## Defining function to calculate auc
        AreaUnderSpline <- function(Time, OD600, minT, maxT){
          spline.area <- stats::integrate(stats::splinefun(Time, OD600, method = "natural"), 
                                          lower = minT, upper = maxT)$value
          return(spline.area)
        }
        
        m <- AreaUnderSpline(data_subset$Time, data_subset$OD600, 1, 73)
        print(m)
        auc_values[nrow(auc_values) + 1,] <- c(i,j,k,l,m)
      }
      
      all_data <- bind_rows(all_data, clean_data)
    }
  }
}

ggplot(all_data %>% filter(Strain == "100"), aes(x = Time, y = OD600)) +
  geom_point() +
  facet_wrap(~Well, nrow = 8)

################################################################################
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
ps1_long = subset(ps1_long, select = -c(Row, Column) )
ps2_long$Well <- paste(ps2_long$Row, ps2_long$Column, sep = "")
ps2_long = subset(ps2_long, select = -c(Row, Column) )

# bind chemical data with corresponding plates & wells
auc_values_all <- data.frame()
for (i in c(1,2)) {
  library_ps <- paste("ps",i,"_long",sep = "")
  auc_values_library <- auc_values[auc_values$Library == i,] %>%
    left_join(get(library_ps), by = "Well")
  auc_values_all <- bind_rows(auc_values_all, auc_values_library)
} 





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

well_arr <- rep2_func_row(well_arr)
well_arr <- rep2_func_col(well_arr)



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

well_arr <- rep3_func_row(well_arr)
well_arr <- rep3_func_col(well_arr)

well_arr$well_ps2 <- paste(well_arr$well_row2, well_arr$well_col2, sep = "")
well_arr$well_ps3 <- paste(well_arr$well_row3, well_arr$well_col3, sep = "")



  
row_alphabet <- c("A","B","C","D","E","F","G","H")
column_number <- c(seq(1,12))
  well_arr <- expand.grid(row_alphabet, column_number)
    colnames(well_arr) <- c("row","col")
  well_arr$well <- paste(well_arr$row, well_arr$col, sep = ":")
  well_arr = subset(well_arr) %>% 
    arrange(well)
  
  
for (i in auc_values_all$Replica) {
  if (i == 1) {
    next
  }
  else if (i == 2) {
    auc_values_all$Well[auc_values_all$Replica == i] <- well_arr$well_ps2
  }
  else {
    auc_values_all$Well[auc_values_all$Replica == i] <- well_arr$well_ps3
  }
}




# calculate control AUC
control_auc <- mean(as.numeric(auc_values[auc_values$Chemical == "DMSO",]$AUC))

auc_values$ratio <- as.numeric(auc_values$AUC)/control_auc

ggplot(auc_values, aes(x = ratio)) + geom_histogram()