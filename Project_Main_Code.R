###########################Section 1###############################
## Cleaning commands and graphs & Setting working directory
cat("\014")
rm(list=ls())
setwd("D:/Masters_project/data/growth-curves")
getwd()

#load packages
for (i in "packages loaded") {
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
  suppressPackageStartupMessages(library(plot3D))
  suppressPackageStartupMessages(library(berryFunctions))
  suppressPackageStartupMessages(library(plotrix))
  suppressPackageStartupMessages(library(gridExtra))
  suppressPackageStartupMessages(library(MASS))
  suppressPackageStartupMessages(library(shape))
  suppressPackageStartupMessages(library(formattable))
  library(lattice)
  suppressPackageStartupMessages(library(vegan))
  library(ggrepel)
  library(ape)
  library(phytools)
  library(ggcor)
  library(glmm)
  library(wrapr)
  library(ggtree)
  print(i)
}

strain_name <- read.csv("aiden-strain-taxonomy.csv")
strain_name <- subset(strain_name, select = c("Species", "Strain.id")) %>% 
  arrange(as.character(Strain.id))


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
strain <- c("88", "100", "186", "322", "333", "350", "353", "374", "380", "390", "442", "448", "487", "527", "565") # unstable : 353 448 error: 390
library <- c("1", "2")
rep <- c("1", "2", "3")

# creating data frame for auc values
all_data <- data.frame()
columns <- c("Strain", "Library", "Replica", "Well", "AUC")
auc_values <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(auc_values) <- columns

# read data and calculate auc data for all wells
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
        
        # calculating auc values
        m <- AreaUnderSpline(data_subset$Time, data_subset$OD600, 1, 73)
        print(m)
        # giving auc of all wells from all plates (replica 1,2,3 & ps1, ps2) by adding values from loops
        auc_values[nrow(auc_values) + 1,] <- c(i,j,k,l,m)
      }
      
      # giving growth curve data of all wells from all plates during 72 hrs by adding values from loops
      all_data <- bind_rows(all_data, clean_data)
    }
  }
}


## load weird data from strain 436, 398, and 85
setwd("D:/Masters_project/data/weird")
filenames <- list.files(getwd(), pattern="*.txt", full.names=TRUE)
ldf <- lapply(filenames, read.table, sep = '\t', header = TRUE, check.names = FALSE)
setwd("D:/Masters_project/data/growth-curves")

testsp <- strsplit(filenames, "_")
ldf_data_all <- data.frame()
ldf_auc_sum <- data.frame()
ldf_auc <- data.frame("Strain", "Library", "Replica", "Well", "AUC")
colnames(ldf_auc) <- c("Strain", "Library", "Replica", "Well", "AUC")

for (i in 1:length(ldf)) {
  j <- testsp[[i]][5]
  k <- testsp[[i]][7]
  l <- testsp[[i]][3]
  ldf_data <- as.data.frame(ldf[[i]])
  ldf_data <- clean.time(remove.temp(ldf_data))
  ldf_data <- pivot_longer(ldf_data, cols = c(2:97), names_to = "Well", values_to = "OD600") %>%
    arrange(Well)
  ldf_data$Strain <- as.character(l)
  ldf_data$Library <- as.character(j)
  ldf_data$Replicate <- as.character(k)
  
  for (p in unique(ldf_data$Well)){
    
    ## Get data for each well
    ldf_data_subset <- ldf_data %>% filter(Well == p)
    Time <- ldf_data_subset$Time 
    OD600 <- ldf_data_subset$OD600
    
    # calculating auc values
    m <- AreaUnderSpline(ldf_data_subset$Time, ldf_data_subset$OD600, 1, 73)
    print(m)
    # giving auc of all wells from all plates (replica 1,2,3 & ps1, ps2) by adding values from loops
    ldf_auc$Strain <- as.character(l)
    ldf_auc$Library <- as.character(j)
    ldf_auc$Replica <- as.character(k)
    ldf_auc$Well <- p
    ldf_auc$AUC <- as.character(m)
    
    ldf_auc_sum <- bind_rows(ldf_auc_sum, ldf_auc)
  }
  
  ldf_data_all <- bind_rows(ldf_data_all, ldf_data)
}

all_data <- bind_rows(all_data, ldf_data_all)
auc_values <- bind_rows(auc_values, ldf_auc_sum)


## binding extra data from Tom S.
extra <- read.csv("spline-fits.csv")
extra <- subset(extra, 
                Strain %in% c("331", "371", "74"), # 302 306 not sure if should be added.
                c("Strain", "Library", "Replica", "Well", "AUC")) %>%
  arrange(Well) %>% arrange(Replica) %>% arrange(Library) %>% arrange(Strain)

extra[extra$Library == "PS1",]$Library <- 1
extra[extra$Library == "PS2",]$Library <- 2
extra$Replica <- as.character(extra$Replica)
extra$AUC <- as.character(extra$AUC)

auc_values <- bind_rows(auc_values, extra)

strain <- NULL
library <- NULL
rep <- NULL

################################Section 2.5###################################
# try to plot graphs
# for all strains and replicates
ggplot(all_data, aes(x = Time, y = OD600, col = Strain)) +
  geom_point()

# separately
ggplot(all_data, aes(x = Time, y = OD600, col = Library)) +
  geom_point() +
  facet_wrap(~Strain, nrow = 8)

# for specific strain (or library/replica)
# all together
ggplot(all_data %>% filter(Strain == "100" & Replicate == "1" & Library == "1"), aes(x = Time, y = OD600, col = Well)) +
  geom_point() +
  theme(legend.position = "none")

# separately
ggplot(all_data %>% filter(Strain == "100" & Replicate == "1" & Library == "1"), aes(x = Time, y = OD600)) +
  geom_point() +
  theme(legend.position = "none") +
  facet_wrap(~Well, nrow = 8)


################################Section 3###################################
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


#############################Section 4#####################################
# getting auc ratio values for each well
# delete data from 8th column of strain 100 library 1 replica 3, due to experiment error
well_delete <- expand.grid(LETTERS[1:8], 8)
well_delete$WD <- paste(well_delete$Var1, well_delete$Var2, sep = "")
auc_values_all <- auc_values_all[!(auc_values_all$Library %in% 1 & auc_values_all$Replica %in% 3 &
                                   auc_values_all$Strain %in% 100 &
                                   auc_values_all$Well %in% well_delete$WD),]

# delete Well B3 from strain 306 library 1 replica 3
#auc_values_all <- auc_values_all[!(auc_values_all$Strain == 306 & 
#                                     auc_values_all$Library == 1 & 
#                                     auc_values_all$Replica == 3 & 
#                                     auc_values_all$Well == "B3"),]

# calculate control AUC and AUC ratio
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


##############################Section 5#####################################
# get graph overview of auc ratio of each strain
# histogram of auc ratio
histogram_AUC_ratio <- ggplot(auc_values_all, aes(x = AUC_ratio, fill = Chemical, binwidth = 10)) +
  geom_histogram(bins = 100) +
  facet_wrap(~Strain, nrow = 8) +
  theme(legend.position = "none")
histogram_AUC_ratio

# get a heatmap of auc ratio of different strains and chemicals
# get means of auc ratios of all replicas
auc3 <- aggregate(auc_values_all$AUC_ratio,
                  by = list(Strain = auc_values_all$Strain,
                          Chemical = auc_values_all$Chemical),
                  data = auc_values_all,
                  FUN = mean) %>%
  arrange(Strain)

auc3$x <- as.numeric(auc3$x)
auc_wide <- pivot_wider(auc3, id_cols = Chemical, names_from = Strain, values_from = x)
auc_wide <- column_to_rownames(auc_wide, "Chemical")

auc_wide_hmp <- auc_wide
colnames(auc_wide_hmp) <- strain_name$Species

col_fun <- colorRamp2(c(-1, -0.5, -0.2, 0, -0.2, 0.5, 1), c("brown", "red", "pink", "white", "slateblue", "blue", "purple"))
Heatmap(log10(as.matrix(auc_wide_hmp)), col = col_fun, border = "grey", show_row_names = FALSE, column_names_rot = 80)

# do dunnett tests as statistical tests. This may take hours.
dta <- data.frame()
for (i in unique(auc_values_all$Strain)) {
  dt1 <- DunnettTest(as.numeric(AUC) ~ Chemical, 
                     data = auc_values_all[auc_values_all$Strain %in% i,], 
                     control = "DMSO")
  test2 <- as.data.frame(dt1$DMSO)
  test2$strain <- i
  dta <- bind_rows(dta, test2)
}

dta <- rownames_to_column(dta) # tried %>% head()
for (i in unique(dta$strain)){
  dta$rowname[dta$strain == i] <- auc3[auc3$Strain == i & auc3$Chemical != "DMSO",]$Chemical
}

# show results of dunnett tests in graphs
ggplot(dta, aes(x = rowname, y = log10(pval), col = strain)) +
  geom_point()

dta$AUC_Sig <- 0
if (nrow(dta[dta$pval < 0.05 & dta$diff > 0,]) == 0) {
  print("No significantly larger AUC")
} else {
  dta[dta$pval < 0.05 & dta$diff > 0,]$AUC_Sig <- "Significantly larger"
}
dta[dta$AUC_Sig != "Significantly larger" & dta$pval < 0.05,]$AUC_Sig <- "Significantly smaller"
dta[dta$AUC_Sig != "Significantly larger" & dta$AUC_Sig != "Significantly smaller",]$AUC_Sig <- "Nonsignificant"

dta_wide <- pivot_wider(dta, id_cols = rowname, names_from = strain, values_from = AUC_Sig)
dta_wide <- column_to_rownames(dta_wide, "rowname")
dta_wide_hmp <- dta_wide
colnames(dta_wide_hmp) <- strain_name$Species

Heatmap(as.matrix(dta_wide_hmp), col = c("white","slateblue","firebrick"), border = "grey", rect_gp = gpar(col = "grey"), show_row_names = FALSE, column_names_rot = 93)

#############################Section 5.5#######################################
# get rid of chemical data with no significant effect on any of the strains
dta_wide_clean <- dta_wide %>% 
  filter_all(any_vars(. %in% c("Significantly smaller") | . %in% c("Significantly larger")))

# draw the heatmap of auc values with cleaned data
clean_chemical <- row.names(dta_wide_clean)
auc_wide_clean <- filter(auc_wide, row.names(auc_wide) %in% clean_chemical)
auc_wide_clean_hmp <- auc_wide_clean
colnames(auc_wide_clean_hmp) <- strain_name$Species
Heatmap(log10(as.matrix(auc_wide_clean_hmp)), col = col_fun, border = "grey", rect_gp = gpar(col = "grey"), row_names_side = "left", show_row_dend = FALSE, column_names_rot = 93)

# draw the heatmap of the result of dunnett's test for the cleaned data
Heatmap(as.matrix(dta_wide_clean), col = c("white","black","gray"), border = "grey", rect_gp = gpar(col = "grey"), column_names_rot = 85)


#############################Section 6########################################
# clustering analyses
## Hierarchical clustering using hclust without redundant chemicals
# read chemical grouping details by their targets
chemical_d <- read.csv("chemicaldetails.csv")
chemical_t <- subset(chemical_d, select = c("Chemical","Target"))
chemical_t <- filter(chemical_t, Chemical %in% clean_chemical) %>% arrange(Chemical)

# make hierarchical clusters using hclust
clust_dist_2 <- dist(log10(auc_wide_clean), method = "euclidean") # measure the distance between the data (euclidean distance)
clust_output_2 <- hclust(d = clust_dist_2, method = "average") # do clustering using average linkage as fusion method 

# draw a scree plot to see the 'elbow point'
ggplot(clust_output_2$height %>%
                     as_tibble() %>%
                     add_column(groups = length(clust_output_2$height):1) %>%
                     rename(height = value),
               aes(x = groups, y = height)) +
       geom_point() +
       geom_line()

# Silhouette Widths 

# draw the dendrogram of clustering
fviz_dend(clust_output_2, 
          cex = 0.7, 
          k = 7, 
          palette = "jco", 
          rect = TRUE, 
          rect_border = "green", 
          rect_lty = 7, 
          horiz = FALSE)

# cut the groups and see the comparison between the clusters and original groups
test_tree <- cutree(clust_output_2, k = 7)
plot(test_tree)
table(chemical_t$Target, test_tree)

## data performing of hclust
# plot the table of comparison
test_tree_df <- as.data.frame(test_tree) %>% rownames_to_column(var = "Chemical")
test_tree_df <- left_join(test_tree_df, chemical_t, by = "Chemical")
ggplot(test_tree_df, aes(x = Target, y = test_tree)) +
  geom_count()

# plot a 3D graph of the table
# data preparation
testqqn <- as.data.frame(table(chemical_t$Target, test_tree))
testqqn <- pivot_wider(testqqn, 
                       id_cols = Var1, 
                       names_from = test_tree, 
                       values_from = Freq) %>% 
  column_to_rownames(var = "Var1")

# add buffer rows and columns
testqqn_b <- testqqn
for (i in 1:(ncol(testqqn_b) - 1)) {
  testqqn_b <- add_column(testqqn_b, buffer = 0, .after = as.character(i))}
testqqn_b <- rownames_to_column(testqqn_b)
testqqn_rn <- testqqn_b[testqqn_b$rowname != testqqn_b[nrow(testqqn_b),]$rowname,]$rowname
for (i in testqqn_rn) {
  testqqn_b <- insertRows(testqqn_b,
                          (which(testqqn_b$rowname %in% i) + 1),
                          new = 0)
  testqqn_b$rowname[testqqn_b$rowname %in% 0] <- as.character(paste(i, "buffer", sep = "."))}
testqqn_b <- column_to_rownames(testqqn_b, var = "rowname")

# plot the 3D histogram
hist3D(z = log(as.matrix(testqqn_b) + 1), phi = 30, theta = 50,
       shade = 0.1, 
       border = "black",
       main = "Count Index", clab = c("log(count) + 1"), 
       breaks = seq(0,3, by = 0.05))

# Define a function to add 3D bars
add_3Dbar <- function(p, x,y,z, width = 0.4) {
  w <- width
  add_trace(p, type="mesh3d",
            x = c(x-w, x-w, x+w, x+w, x-w, x-w, x+w, x+w),
            y = c(y-w, y+w, y+w, y-w, y-w, y+w, y+w, y-w),
            z = c(0, 0, 0, 0, z, z, z, z),
            i = c(7, 0, 0, 0, 4, 4, 2, 6, 4, 0, 3, 7),
            j = c(3, 4, 1, 2, 5, 6, 5, 5, 0, 1, 2, 2),
            k = c(0, 7, 2, 3, 6, 7, 1, 2, 5, 5, 7, 6),
            facecolor = rep(toRGB(viridisLite::plasma(6)), each = 2)) 
}

# Draw the 3D histogram
fig <- plot_ly()
for (k1 in 1:nrow(as.matrix(testqqn))) {
  for (k2 in 1:ncol(as.matrix(testqqn))) {
    fig <- fig %>% add_3Dbar(k1,k2,as.matrix(testqqn)[k1,k2])
  }
}
  fig


##############################Section 6.1####################################
## Expectation-Maximization (EM) algorithm for model-based clustering using mclust (soft assignment)
## Gaussian mixture distribution assumed, choosing from 14 models
# prepare the data
auc_wide_clean_m <- rownames_to_column(auc_wide_clean)
auc6 <- pivot_longer(auc_wide_clean_m, cols = c(unique(colnames(auc_wide_clean_m[-1]))), names_to = "strain", values_to = "auc_ratio")
  auc6 <- auc6 %>% relocate(rowname, .after = strain) %>% arrange(strain)
  auc6$Target <- ""
    for (i in unique(auc6$strain)) {
      auc6[auc6$strain == i,]$Target <- chemical_t$Target
    }
auc7 <- auc6[,c(1,2,3)]
auc7_wide <- pivot_wider(auc7, id_cols = rowname, names_from = strain, values_from = auc_ratio)
auc7_wide <- column_to_rownames(auc7_wide, var = "rowname")
Xcl5 <- data.matrix(auc7_wide)
Xcl5log <- log10(Xcl5)

# build the clustering model
cl5 <- Mclust(Xcl5log, G = 1:7) # larger number of components increases BIC while causes 'overclustering' problem. Clusters of single (or very few) clustering object may be taken into account and scatter around.

# plot the BIC curves
plot.Mclust(cl5, what = "BIC")

# show properties of the model
cl5$modelName # the optimal model with highest BIC
cl5$n # number of chemicals
cl5$G # number of groups
cl5$bic # the highest BIC

# plot the results of clustering
plot(cl5, what = "classification") # clustering with each pair of variables (strains)
table(chemical_t$Target, cl5$classification) # comparison between clustering and groups by targets
# get adjusted Rand index to evaluate the comparison of the two clusterings (ARI; Hubert and Arabie, 1985)
adjustedRandIndex(chemical_t$Target, cl5$classification)

# dimension reduction into 2 dimensions
cl5_dr <- MclustDR(cl5, lambda = 1)
# plot the dimension reduction graph
plot(cl5_dr, what = "classification", colors = c("red", "green", "skyblue"))
summary(cl5_dr)

# plot the clustering (DR) in ggplot2
# preparation of data
mclust_dr_table <- as.data.frame(cl5_dr$dir) # take the values of dir 1 and 2.
mlabels <- rownames(Xcl5log)
mclust_dr_table$Chemical <- mlabels
cl5cl <- as.data.frame(cl5_dr$classification) %>% rownames_to_column(var = "Chemical") # take the clusters of each chemical
mclust_dr_table <- mclust_dr_table %>% left_join(cl5cl, by = "Chemical") # assign the clusters to the coordinates of chemical on dir 1&2
colnames(mclust_dr_table)[colnames(mclust_dr_table) == "cl5_dr$classification"] <- "Cluster"
mclust_dr_table <- mclust_dr_table[c("Dir1","Dir2","Chemical","Cluster")]

# plot the MclustDR graph in ggplot2 with density plot of dir 1 and dir 2
cldr_plot <- ggplot(mclust_dr_table, aes(x = Dir1, y = Dir2, label = Chemical, color = Cluster)) +
  geom_point(size = 3) +
  xlim(-0.1, 0.2) +
  geom_label_repel(aes(fill = factor(Cluster), hjust = -0.1), colour = "white", fontface = "bold", max.overlaps = 50, label.padding = 0.12, label.size = 0.17) + 
  scale_color_manual(values = c('red', 'green', 'skyblue')) + 
  theme(legend.position=c(0,1), legend.justification=c(0,1))

xdensity <- ggplot(mclust_dr_table, aes(x = Dir1, fill = Cluster, after_stat(scaled))) +
  geom_density(alpha = .5, kernel = "gaussian") + 
  xlim(-0.1, 0.2) +
  scale_fill_manual(values = c('red', 'green', 'skyblue')) + 
  theme(legend.position = "none")

ydensity <- ggplot(mclust_dr_table, aes(y = Dir2, fill = Cluster, after_stat(scaled))) + 
  geom_density(alpha = .5, kernel = "gaussian") + 
  scale_fill_manual(values = c('red', 'green', 'skyblue')) + 
  theme(legend.position = "none")

blankPlot <- ggplot() + geom_blank(aes(1, 1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

grid.arrange(xdensity, blankPlot, cldr_plot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

# plot Gaussian distributed density plot fitted by dir 1 and 2
plot(cl5_dr, what = "density", dimens = 1)
plot(cl5_dr, what = "density", dimens = 2)

# plot a 3D density graph
den3d <- kde2d(mclust_dr_table$Dir1, mclust_dr_table$Dir2)
persp(den3d, col = drapecol(den3d$z), theta = 30, phi = 30)

# get a final result table for chemical with both target groups and clusters
chemical_table_final <- left_join(mclust_dr_table, chemical_t, by = "Chemical") %>% 
  subset(select = c("Chemical", "Cluster", "Target", "Dir1", "Dir2")) %>%
  arrange(Chemical)

formattable(chemical_table_final)
write.csv(chemical_table_final, "chemical_table_final.csv")


############################Section 7######################################
## Ordination - PCA
# prepare the data for pca
auc_values_pca <- auc_values_all
#auc_values_pca <- auc_values_pca[auc_values_pca$Strain != "306" & auc_values_pca$Strain != "302",]
auc_values_pca$StrainRep <- paste(auc_values_pca$Strain, auc_values_pca$Replica, sep = "_")
auc_values_pca <- subset(auc_values_pca, Chemical != "DMSO", select = c("StrainRep", "Chemical", "AUC_ratio"))
auc_values_pca <- pivot_wider(auc_values_pca, id_cols = StrainRep, names_from = Chemical, values_from = AUC_ratio)
auc_values_pca <- subset(auc_values_pca, select = c("StrainRep", clean_chemical)) %>% 
  column_to_rownames(var = "StrainRep")
auc_values_pca <- auc_values_pca %>% 
  filter_all(all_vars(. %nin% NA))

# run the pca using rda
auc_pca <- prcomp(log10(auc_values_pca))
PCAloadings <- data.frame(Variables = rownames(auc_pca$rotation), auc_pca$rotation)

# show scree plot
screeplot(auc_pca)
abline(a=1,b=0)

# show the biplot
kolr <- 1
chemical_table_final$Color <- kolr
for (i in unique(chemical_table_final$Target)) {
kolr <- kolr + 1
chemical_table_final[chemical_table_final$Target %in% i,]$Color <- kolr
}
chemical_table_final$Color2 <- as.numeric(chemical_table_final$Cluster) + 1

suppressWarnings(ggplot2:::print.ggplot(
ggplot(as.data.frame(auc_pca[["x"]][,1:2]), aes(PC1, PC2, label = rownames(auc_pca[["x"]][,1:2]))) + 
  geom_point(size = 2) +
  geom_segment(data = as.data.frame(auc_pca[["rotation"]])[,1:2], 
               aes(0, 0, xend = PC1*3, yend =PC2*3, label = NULL),
               arrow = arrow(length = unit(1/2, "picas")),
               col = (chemical_table_final$Color2), 
               linewidth = 2) +
  geom_label_repel(data = PCAloadings, 
                  x = (PCAloadings$PC1) * 3.1, 
                  y = (PCAloadings$PC2) * 3.1, 
                  label = PCAloadings$Variables, 
                  col = chemical_table_final$Color2, 
                  hjust = 1, vjust = 0.4,
                  max.overlaps = 100,
                  alpha = 0.6) +
  theme_gray() + 
  xlim(-2,2) +
  ylim(-2,2)))


suppressWarnings(ggplot2:::print.ggplot(
ggplot(as.data.frame(auc_pca[["x"]][,1:2]), aes(PC1, PC2*0, label = rownames(auc_pca[["x"]][,1:2]))) + 
  geom_point(size = 2) +
  geom_segment(data = as.data.frame(auc_pca[["rotation"]])[,1:2], 
               aes(0, 0, xend = PC1*3, yend =PC2*0, label = NULL),
               arrow = arrow(length = unit(0.8, "picas")),
               col = (chemical_table_final$Color2), 
               linewidth = 2) +
  theme_gray() + 
  xlim(-2,2) +
  ylim(-0.5,0.5)))


###################################Section 8###################################
# mantel test of phylogenetic and phenotypic (reaction to chemicals) grouping
# extract the distances from the phylogenetic tree
pltree <- read.tree("phylo.io_n.nwk")
plotTree(pltree)
new_pltree <- reroot(pltree, node = 25)
new_pltree <- ladderize(new_pltree)
write.tree(new_pltree, file = "phylo.io_new.nwk")
new_pltree <- read.tree("phylo.io_new.nwk")
plotTree(new_pltree)
phylo.dist <- cophenetic.phylo(new_pltree)
col.order <- rownames(phylo.dist)

# Use "upper", to only get the upper half of the matrix (i.e. 1 number for each pair)
# get phenotypic matrix distance
auc_mt <- log10(t(auc_wide_clean[col.order]))
pheno.dist <- as.matrix(dist(as.matrix(auc_mt), diag = TRUE, upper = TRUE, method = "euclidean"))

## MAKE SURE THE COLUMNS ARE IN THE SAME ORDER FOR BOTH DISTANCE MATRICES
#pheno.dist <- pheno.dist[col.order,col.order]
#phylo.dist <- phylo.dist[col.order,col.order]

# do mantel test with all chemical data
auc_mt_all <- log10(t(auc_wide[col.order]))
pheno.dist_all <- as.matrix(dist(as.matrix(auc_mt_all), diag = TRUE, upper = TRUE, method = "euclidean"))
pheno.dist_all <- pheno.dist_all[col.order,col.order]
plmt_all <- mantel(phylo.dist, pheno.dist_all, method = "kendall", permutations = 9999)
plmt_all

# then do the mantel test:
plmt <- mantel(phylo.dist, pheno.dist, method = "kendall", permutations = 9999)
plmt
plmt <- mantel_test(phylo.dist, pheno.dist, 
                    method = "kendall", 
                    permutations = 9999, 
                    use = "pairwise")

# plot the results of mantel test
plmt <- plmt %>% 
  mutate(lty = cut(r, breaks = c(-Inf, 0, Inf), 
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 0.1, 1),
                   labels = c("< 0.01", "< 0.05", "<0.10", ">= 0.10"),
                   right = FALSE, include.lowest = TRUE))

quickcor(pheno.dist, type = "upper") +
  geom_square() +
  anno_link(plmt, mapping = aes(colour = col,
                                size = r^3,
                                linetype = lty)) +
  scale_fill_gradient2n() +
  scale_size_area(max_size = 3) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "Mantel's p", order = 2),
    size = guide_legend(title = "Mantel's r", order = 3),
    linetype = "none"
  )


mtt_data_pheno <- rownames_to_column(as.data.frame(pheno.dist), var = "Initial.Strain")
mtt_data_pheno <- pivot_longer(mtt_data_pheno, 
                        cols = 2:ncol(mtt_data_pheno), 
                        names_to = "Paired.Strain", 
                        values_to = "dist.phenotypic")

mtt_data_phylo <- rownames_to_column(as.data.frame(phylo.dist), var = "Initial.Strain")
mtt_data_phylo <- pivot_longer(mtt_data_phylo, 
                        cols = 2:ncol(mtt_data_phylo), 
                        names_to = "Paired.Strain", 
                        values_to = "dist.phylogenetic")

mtt_data <- bind_cols(mtt_data_pheno, mtt_data_phylo[3])
par(mfrow = c(1,1))
plot(mtt_data$dist.phenotypic, mtt_data_phylo$dist.phylogenetic)

# density plot of distances from both matrices
ggplot(data = mtt_data[!(mtt_data$dist.phenotypic == 0),], aes(x = scale(dist.phenotypic))) +
  geom_density(fill = "grey", col = "grey") +
  geom_density(aes(x = scale(dist.phylogenetic), col = "red", fill = "red", alpha = 1/15))


############################Section 9##########################################
# get phylogenetic results
# get auc data align with strain order in phylo tree
auc_wide_t <- log10(as.data.frame(t(auc_wide_clean))) %>% rownames_to_column(var = "Strain")
strain_order <- match_order(auc_wide_t$Strain, new_pltree$tip.label)
auc_wide_t <- auc_wide_t[strain_order,]

auc_wide_t_nl <- as.data.frame(t(auc_wide_clean)) %>% rownames_to_column(var = "Strain")
auc_wide_t_nl <- auc_wide_t_nl[strain_order,]

# set data and names for single chemicals from 3 different clusters
svl1 <- setNames(auc_wide_t$ABAMECTIN, auc_wide_t$Strain)
svl2 <- setNames(auc_wide_t$FENPYROXIMATE, auc_wide_t$Strain)
svl3 <- setNames(auc_wide_t$DICHLOROPHEN, auc_wide_t$Strain)

svl1nl <- setNames(auc_wide_t_nl$ABAMECTIN, auc_wide_t_nl$Strain)
svl2nl <- setNames(auc_wide_t_nl$FENPYROXIMATE, auc_wide_t_nl$Strain)
svl3nl <- setNames(auc_wide_t_nl$DICHLOROPHEN, auc_wide_t_nl$Strain)

# test phylogenetic signals
phylosig(new_pltree,svl1nl,method="lambda",test=TRUE)
phylosig(new_pltree,svl2nl,method="lambda",test=TRUE)
phylosig(new_pltree,svl3nl,method="lambda",test=TRUE)

# plot phylo tree with color for signal
fit_tree <- fastAnc(new_pltree, svl1, vars=TRUE, CI=TRUE)
print(fit_tree, printlen = 10)
trpoj1 <- contMap(new_pltree, svl1, plot=FALSE)
par(mfrow = c(1,1))
plot(trpoj1, type="fan", legend = 0.7*max(nodeHeights(new_pltree)),
     sig = 1, fsize = c(0.9, 0.9))

fit_tree <- fastAnc(new_pltree, svl2, vars=TRUE, CI=TRUE)
print(fit_tree, printlen = 10)
trpoj2 <- contMap(new_pltree, svl2, plot=FALSE)
par(mfrow = c(1,1))
plot(trpoj2, type="fan", legend = 0.7*max(nodeHeights(new_pltree)),
     sig = 1, fsize = c(0.9, 0.9))

fit_tree <- fastAnc(new_pltree, svl3, vars=TRUE, CI=TRUE)
print(fit_tree, printlen = 10)
trpoj3 <- contMap(new_pltree, svl3, plot=FALSE)
par(mfrow = c(1,1))
plot(trpoj3, type="fan", legend = 0.7*max(nodeHeights(new_pltree)),
     sig = 1, fsize = c(0.9, 0.9))


# plot phylomorphospace
rownames(auc_wide_t) <- c(1:21)
auc_wide_t <- column_to_rownames(auc_wide_t, var = "Strain")
obj<-fancyTree(new_pltree,type="scattergram", 
               X = auc_wide_t[,c("ABAMECTIN", "FENPYROXIMATE", "DICHLOROPHEN")])

