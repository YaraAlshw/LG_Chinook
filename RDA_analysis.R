#RDA analysis
#Written by Yara Alshwairikh based on script from http://popgen.nescent.org/2018-03-27_RDA_GEA.html
#Script consists of 3 parts
#Note that 
#freqs.RData is available on Dryad (doi to be added later)

library("tidyverse")
library("dplyr")
library("psych")
library("vegan")
library("miceadds")


###---Start script part 1---###

load.Rdata(filename = "freqs.RData", "maf") #Load the allele frequency matrix

maf_subset <- maf[, sample(ncol(maf), 100000)] #This RDA script is computationally heavy and will crash on a regular computer. To try it out first, you can use this line of code to susbet the allele frequency matrix and run the full script to make sure everything works. When ready to run the full analysis, run it on the cluster

env_site_five <- read.csv("env_site_five.csv") #In addition to subsetting the allele frequency matric, here we will only use FIVE environmental variables to run the analysis and make sure everything works before moving to the cluster

#Run RDA using vegan package
ots.rda <- vegan::rda(maf_subset ~ ., data=env_site_five, Scale=T)
ots.rda #summary of RDA

vegan::RsquareAdj(ots.rda) #R squared values

summary(ots.rda)$concont #proportion of variance explaines by each axis
screeplot(ots.rda) #visualize the canconical eignevalues

vegan::vif.cca(ots.rda) #Check VIF for the variables

#Check the FULL RDA model for significance
signif.full <- anova.cca(ots.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

#Check each constrained axis for significance
signif.axis <- anova.cca(ots.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis #Check to see which axes have p = 0.05. Ideally this shoud match the results in the screeplots

###---End script part 1---###


###---Start script part 2---###

#Find candidate SNPs: this part relies on knowledge of which axes are significant. 
#Currently, the script is written assuming axes 1, 2, and 3 are significant

load.rda <- summary(ots.rda)$species[,1:3]
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
} #function to find the outlier SNPs

#apply the outliers() function to each axis
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
cand3 <- outliers(load.rda[,3],3) 

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand #total number of candidate SNPs

#Organzie the results into a datadram with the axis, SNP, laoding, a correlation with each environmental variable
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=5)  # create 5 columns for 5 predictors
colnames(foo) <- c("Wsp_range", "HLI_Site", "Bio15_site", "Bio3_site", "Bio5_site")
                   
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- maf[,nam]
  foo[i,] <- apply(env_site_five,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  

length(cand$snp[duplicated(cand$snp)]) #check for duplicate SNPs

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) #check for duplicates on axis 1
table(foo[foo[,1]==2,2]) #check for duplicates on axis 2
table(foo[foo[,1]==3,2]) #check for duplicates on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

#To find which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,9] <- names(which.max(abs(bar[4:8]))) # gives the variable
  cand[i,10] <- max(abs(bar[4:8]))              # gives the correlation
}
#note: the 4:8 are the columns where the environmental variables are, must manually adjust this for every dataset
#can confirm through str(cand). MUST change cand[i,#] to specify the column right after col 8, 
#which here is column 9, or else you get an error 
#do the same thing for the correlation cand[ ]

#assign column names
colnames(cand)[9] <- "predictor"
colnames(cand)[10] <- "correlation"

table(cand$predictor) #lists top associations 

#Write the output into a .csv
write.csv(cand, "RDA_cand_site.csv")

###---End script part 2---###

###---Start script part 3---###

###---Plots---###
#Just like part 2, this part 3 depends on how many axes are significant and we want to plot
#Currently, this script assumes the first 3 axes are significant

#A) Triplots
#read a "env" dataframe containing 1 column for "population", and 1 column for "type"
env <- read.csv("env.csv")
#set levels
levels(env$population) <- c("Yakima","Wenatchee","Lyons Ferry","Methow","Deschutes","Clearwater", "Priest")
levels(env$type) <- c("Fall","Summer")

population <- env$population
type <- env$type

#set colors for each population
bg <- c("#A50026","#DC3D2D", "#F57D4A", "#354B99", "#FDB366", "#6DA5CC", "#FED98B")

# axes 1 & 2
png(file="RDA_plot1.png") 
plot(ots.rda, type="n", scaling=3)
points(ots.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(ots.rda, display="sites", pch=c(21, 24) [type], cex=1.3, col="gray32", scaling=3, bg=bg[population]) # the populations (7 pops)
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors (5)
legend("bottomright", legend=levels(population), bty="n", col="gray32", pch=c(21, 24) [type], cex=1, pt.bg=bg)
legend("bottomleft", legend=levels(type), bty="n", col="gray32", pch=c(1, 2), cex=1, pt.bg=type) 
dev.off()

