# Mantel test to compare phylogenetic distance with phenotypic distance (i.e. do closely related strains have similar responses to the chemicals?)

library(vegan) # good package for doing the mantel test
library(ape) # for phylogenetic stuff
library(phytools) # for more phylogenetic stuff
library(ggcor)

# You need a distance matrix of the responses of the species to all the chemicals
# Use "upper", to only get the upper half of the matrix (i.e. 1 number for each pair)
auc_mt <- log10(t(auc_wide_clean))
pheno.dist <- as.matrix(dist(as.matrix(auc_mt), diag = TRUE, upper = TRUE, method = "euclidean"))

# You also need to extract the distances from the phylogenetic tree
# cophenetic.phylo in the ape package can do this
pltree <- read.tree("phylo.io_n.nwk")
phylo.dist <- cophenetic.phylo(pltree)

## VERY IMPORTANT:
## YOU WILL NEED TO MAKE SURE THE COLUMNS
## ARE IN THE SAME ORDER FOR BOTH DISTANCE MATRICES!!!

# example:
col.order <- colnames(pheno.dist)
pheno.dist <- pheno.dist[col.order,col.order]
phylo.dist <- phylo.dist[col.order,col.order]

# then do the mantel test:
plmt <- mantel_test(phylo.dist, pheno.dist, 
                    method = "kendall", 
                    permutations = 9999, 
                    use = "pairwise")


plmt <- plmt %>% 
  mutate(lty = cut(r, breaks = c(-Inf, 0, Inf), 
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))

quickcor(pheno.dist, type = "upper") +
  geom_square() +
  anno_link(plmt, mapping = aes(colour = col,
                             size = r,
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






# suggest kendall as its non parametric and your data are probably not normally distributed
# but you could also try pearson or spearman if you think they suit your data better

library("wrapr")
# get auc data align with strain order in phylo tree
auc_wide_t <- as.data.frame(t(auc_wide_clean)) %>% rownames_to_column(var = "Strain")
strain_order <- match_order(auc_wide_t$Strain, pltree$tip.label)
auc_wide_t <- auc_wide_t[strain_order,]

# set data and names for single chemicals from 3 different clusters
svl1 <- setNames(auc_wide_t$FLUBENDIAMID, auc_wide_t$Strain)
svl2 <- setNames(auc_wide_t$BIXAFEN, auc_wide_t$Strain)
svl3 <- setNames(auc_wide_t$DICHLOROPHEN, auc_wide_t$Strain)

# plot phylo tree with color for signal
fit_tree <- fastAnc(pltree, svl1, vars=TRUE, CI=TRUE)
print(fit_tree, printlen = 10)
trpoj <- contMap(pltree, svl1, plot=FALSE)
par(mfrow = c(1,1))
plot(trpoj, type="fan", legend = 0.7*max(nodeHeights(pltree)),
     sig = 2, fsize = c(0.9, 0.9))

fit_tree <- fastAnc(pltree, svl2, vars=TRUE, CI=TRUE)
print(fit_tree, printlen = 10)
trpoj <- contMap(pltree, svl2, plot=FALSE)
par(mfrow = c(1,1))
plot(trpoj, type="fan", legend = 0.7*max(nodeHeights(pltree)),
     sig = 2, fsize = c(0.9, 0.9))

fit_tree <- fastAnc(pltree, svl2, vars=TRUE, CI=TRUE)
print(fit_tree, printlen = 10)
trpoj <- contMap(pltree, svl2, plot=FALSE)
par(mfrow = c(1,1))
plot(trpoj, type="fan", legend = 0.7*max(nodeHeights(pltree)),
     sig = 2, fsize = c(0.9, 0.9))


# plot phylomorphospace
rownames(auc_wide_t) <- c(1:21)
auc_wide_t <- column_to_rownames(auc_wide_t, var = "Strain")
obj<-fancyTree(pltree,type="scattergram", X = auc_wide_t[,c("BIXAFEN", "FLUBENDIAMID", "DICHLOROPHEN")])


auc_strain_m <- as.data.frame(t(auc_wide_clean)) %>% rownames_to_column(var = "Strain")
strain_mclust <- Mclust(log10(as.matrix(auc_strain_m[,-1])), G = 2)
auc_strain_m$Cluster <- as.data.frame(strain_mclust$classification)





