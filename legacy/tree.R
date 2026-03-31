#comp <- data.frame(chemical = rownames(as.data.frame(test_tree)))
test_tree_df <- as.data.frame(test_tree) %>% rownames_to_column(var = "Chemical")
test_tree_df <- left_join(test_tree_df, chemical_t, by = "Chemical")
#comp <- left_join(comp, test_tree_df, by = "chemical")
ggplot(test_tree_df, aes(x = Target, y = test_tree)) +
  geom_count()

library(plot3D)
library(berryFunctions)
library(plotrix)
library(gridExtra)
library(MASS)
library(shape)

testqqn <- as.data.frame(table(chemical_t$Target, test_tree))
testqqn <- pivot_wider(testqqn, id_cols = Var1, names_from = test_tree, values_from = Freq) %>% column_to_rownames(var = "Var1")


for (i in 1:3) {
  testqqn <- add_column(testqqn, buffer = 0, .after = as.character(i))
}

testqqn <- rownames_to_column(testqqn)
testqqn_rn <- testqqn[testqqn$rowname != testqqn[nrow(testqqn),]$rowname,]$rowname

for (i in testqqn_rn) {
  testqqn <- insertRows(testqqn, 
      (which(testqqn$rowname %in% i) + 1),
      new = 0)
  testqqn$rowname[testqqn$rowname %in% 0] <- as.character(paste(i, "buffer", sep = "."))
}

testqqn <- column_to_rownames(testqqn, var = "rowname")



hist3D(z = log(as.matrix(testqqn) + 1), phi = 30, theta = 50,
       shade = 0.1, 
       border = "black",
       main = "count", clab = c("counts"), 
       breaks = seq(0,3, by = 0.05))



fig <- plot_ly(z = ~ 3 * (log(as.matrix(testqqn) + 1)), colors = c("slateblue", "darkgreen", "yellow", "orange"))
fig <- fig %>% add_surface()
fig



# 3d heatmap

Heatmap3D(as.matrix(auc_wide_clean), col = col_fun,)


# mclust

plot(cl5_dr, what = "classification")

mclust_dr_table <- as.data.frame(cl5_dr$dir)


mlabels <- rownames(Xcl5)
mclust_dr_table$Chemical <- mlabels

cl5cl <- as.data.frame(cl5_dr$classification) %>% rownames_to_column(var = "Chemical")
mclust_dr_table <- mclust_dr_table %>% left_join(cl5cl, by = "Chemical")
colnames(mclust_dr_table)[colnames(mclust_dr_table) == "cl5_dr$classification"] <- "Cluster"


cldr_plot <- ggplot(mclust_dr_table, aes(x = Dir1, y = Dir2, label = Chemical, color = Cluster)) +
  xlim(-0.25, 0.5) +
  geom_point(size = 3) +
  geom_label(aes(fill = factor(Cluster), hjust = -0.1), colour = "white", fontface = "bold") + 
  scale_color_manual(values = c('red', 'green', 'skyblue')) + 
  theme(legend.position=c(0,1), legend.justification=c(0,1))

xdensity <- ggplot(mclust_dr_table, aes(x = Dir1, fill = Cluster)) +
  xlim(-0.25, 0.5) +
  geom_density(alpha = .5, kernel = "gaussian") + 
  scale_fill_manual(values = c('red', 'green', 'skyblue')) + 
  theme(legend.position = "none")

ydensity <- ggplot(mclust_dr_table, aes(y = Dir2, fill = Cluster)) + 
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

plot(cl5_dr, what = "density")

den3d <- kde2d(mclust_dr_table$Dir1, mclust_dr_table$Dir2)
persp(den3d, col = drapecol(den3d$z), theta = 30, phi = 30)
