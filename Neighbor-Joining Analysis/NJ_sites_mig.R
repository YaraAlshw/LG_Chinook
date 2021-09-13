library(vegan)
library(ape)
library(ggplot2)
library(ggtree)

library(readr)
EnvSite1 = read_csv("env_Site.csv")
EnvMig1 = read_csv("env_migration.csv")


EnvSite1 = data.frame(EnvSite1,row.names="Population")
EnvMig1 = data.frame(EnvMig1,row.names="Population")


##-----------Migration paths: Environmental NJ-----------------##

#not scaled
Migmatrix1 = vegdist(EnvMig1, method = "euclidean") #distance matrix based on euclidean distances
tre <- nj(Migmatrix1) # Neighbor joining analysis
class(tre1)
treladder1 <- ladderize(tre1) # set as ladder
resolved.phy1 <- multi2di(treladder1) #resolve basal polytomy

#scaled
EnvMig1_sc= scale(EnvMig1)
MigMatrix_sc=vegdist(EnvMig1_sc, method = "euclidean")
bionjmg_sc=bionj(MigMatrix_sc)
bionjmglad_sc <- ladderize(bionjmg_sc) # set as ladder
bionj1mgla_sc.phy <- multi2di(bionjmglad_sc) #resolve basal polytomy
plot(bionj1mgla_sc.phy)
axisPhylo()
nodelabels() # add node numbers
tiplabels()  # add tip numbers

##--------------Locations: Environmental NJ -------------------##
EnvSite

#Migration paths relationships
Migmatrix_site1 = vegdist(EnvSite1, method = "euclidean") #distance matrix based on euclidean distances
tre_site1 <- nj(Migmatrix_site1) # Neighbor joining analysis
class(tre_site1)
treladder_site1 <- ladderize(tre_site1) # set as ladder
resolved_site.phy1 <- multi2di(treladder_site1) #resolve basal polytomy


#scaled
EnvSite1_sc= scale(EnvSite1)
SiteMatrix_sc=vegdist(EnvSite1_sc, method = "euclidean")
bionjst_sc=bionj(SiteMatrix_sc)
bionjstlad_sc <- ladderize(bionjst_sc) # set as ladder
bionj1stla_sc.phy <- multi2di(bionjstlad_sc) #resolve basal polytomy
plot(bionj1stla_sc.phy)
axisPhylo()
nodelabels() # add node numbers
tiplabels()  # add tip numbers



#Final Figures------------------ Based on scaled data and bioNJ



#site- colors - small scale
ggtree(bionj1stla_sc.phy) %>% rotate(8)%>% rotate(12)+ geom_tiplab()+ theme_tree(plot.margin = margin(20, 140, 20, 10))+coord_cartesian(clip = 'off') + geom_cladelabel(node=12, label="Summer run",align=TRUE, offset=1.22, color="black", offset.text=.02, fontface="bold") + geom_strip(1,2, barsize=3, color="blue", offset=1.2, offset.text=.02) + geom_cladelabel(node=5, label="Fall run",align=TRUE, offset=1.3, color="black", offset.text=.01, fontface="bold") + geom_strip(7,6, barsize=3, color='red', offset=1.2)+ geom_treescale(x = NULL, y = NULL, width = 0.5, offset = 0.2, color = "black", linesize = 1, fontsize = 3.88, family = "sans")

#migration corridor- colors - Small scale
ggtree(bionj1mgla_sc.phy) %>% rotate(8)%>% rotate(12)+ geom_tiplab()+ theme_tree(plot.margin = margin(20, 140, 20, 10))+coord_cartesian(clip = 'off') + geom_cladelabel(node=11, label="Summer run",align=TRUE, offset=1.22, color="black", offset.text=.02, fontface="bold") + geom_strip(1,2, barsize=3, color="blue", offset=1.2, offset.text=.02) + geom_cladelabel(node=7, label="Fall run",align=TRUE, offset=1.3, color="black", offset.text=.01, fontface="bold") + geom_strip(4,6, barsize=3, color='red', offset=1.22)+ geom_treescale(x = NULL, y = NULL, width = 0.5, offset = 0.2, color = "black", linesize = 1, fontsize = 3.88, family = "sans")



