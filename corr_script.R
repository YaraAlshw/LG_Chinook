#Script written by Travis Seaborn
#This script detaisl the steps that were udnertaken to create the 3 environmental datasets, "Site", "Migration", and "Combination" from an initial list of 140 variables
#The ful lset of 140 variables along with their correlation matrix can be found on Dryad (link to be added later)

library(tidyverse)
library(caret)
library(corrplot)
library(psych)

env_full <- read.csv("Steve_Files/env_final.csv") #this is the 44 after the three removal steps
pop_name <- env_full[,1] #moving out pop names
env_full <- env_full[,-1]

# main correlation plot
cor_full <- corr.test(env_full, method = "spearman", adjust = "none")
cor_full$r[cor_full$p > 0.05] <- 0
corrplot(cor_full$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))

# manual stepwise removal
#with me selecting until threshold
#source("http://www.sthda.com/upload/rquery_cormat.r")
step1 <- rquery.cormat(env_full, type="flatten", graph=FALSE)
#again removing all not significant at the start
step1 <- step1$r
step1$cor[step1$p > 0.05] <- 0
step1o <- step1[order(step1$cor),]
#bio 1, 2 are very problematic, removing from the 44
step1_df <- env_full %>% select(-c(Bio1_mean, Bio1_site, Bio2_mean, Bio2_range, Bio2_site, Bio16_mean, Bio16_rang, Bio16_site))


step2 <- rquery.cormat(step1_df, type="flatten", graph=FALSE)
#again removing all not significant at the start
step2 <- step2$r
step2$cor[step2$p > 0.05] <- 0
step2o <- step2[order(step2$cor),]
#dropping bio17 and DEM
step2_df <- step1_df %>% select(-c(Bio17_mean, Bio17_rang, Bio17_site, DEM_site, DEMmean, DEMrange))

step3 <- rquery.cormat(step2_df, type="flatten", graph=FALSE)
#again removing all not significant at the start
step3 <- step3$r
step3$cor[step3$p > 0.05] <- 0
step3o <- step3[order(step3$cor),]
##dropping roughness and tree cover
step3_df <- step2_df %>% select(-c(Rough_mean, Rough_Site, Treecov_me, Treecov_site))

#checking where we are at:
cor_s3 <- corr.test(step3_df, method = "spearman", adjust = "none")
cor_s3$r[cor_s3$p > 0.05] <- 0
highlyCorrelated <- findCorrelation(cor_s3$r, cutoff=0.8, exact = TRUE)
highlyCorCol <- colnames(step3_df)[highlyCorrelated]
s3_var <- step3_df[, -which(colnames(step3_df) %in% highlyCorCol)]

s3_final <- corr.test(s3_var, method = "spearman")
corrplot(s3_final$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0), insig = "blank")

#going to 1) manually add in some variables of interest and 2) do site and corridor separate, and keep them that way
#load csv of full list of mig and site variables so manual add back in can occur for bio3 and 5. these were created prior to the 44 final list decided upon after the a priori 3 step removal
mig <- read.csv("Steve_Files/env_final_mig.csv")
site <- read.csv("Steve_Files/env_final_site.csv")
mig <- mig[,-1]
site <- site[,-1]


#Now dealing with the adding back in variables from other studies
s3_var_add <- cbind(s3_var, env_full$Bio3_range, env_full$Bio3_site, site$Bio5_site, mig$Bio5_range)
s3_var_new <- corr.test(s3_var_add, method = "spearman")
#check correlation
corrplot(s3_var_new$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0), insig = "blank")

#bio12 is really problematic with 3 and 5
s3_var_add$Bio12_rang <- NULL
s3_var_add$Bio12_site <- NULL
s3_var_add$'env_full$Bio3_range' <- NULL
s3_var_new <- corr.test(s3_var_add, method = "spearman")

corrplot(s3_var_new$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0), insig = "blank")

#this gives the final "full" data set

#now to repeat for site and one for mig, going full automated to start
#may end up with different results here each time run due to automated removal
cor_site <- corr.test(site, method = "spearman", adjust = "none")
cor_site$r[cor_site$p > 0.05] <- 0
highlyCorrelated <- findCorrelation(cor_site$r, cutoff=0.8, exact = TRUE)
highlyCorCol <- colnames(site)[highlyCorrelated]
site_var <- site[, -which(colnames(site) %in% highlyCorCol)]

site_var <- cbind(site_var, site$Bio3_site, site$Bio5_site)
site_var$Bio8_site <- NULL
site_var$Bio14_site <- NULL
site_var$Treecov_site <- NULL

site_final <- corr.test(site_var, method = "spearman")
corrplot(site_final$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0), insig = "blank")

cor_mig <- corr.test(mig, method = "spearman", adjust = "none")
cor_mig$r[cor_mig$p > 0.05] <- 0
highlyCorrelated <- findCorrelation(cor_mig$r, cutoff=0.8, exact = TRUE)
highlyCorCol <- colnames(mig)[highlyCorrelated]
mig_var <- mig[, -which(colnames(mig) %in% highlyCorCol)]
#manual add back in of 7, 11, and dist
mig_var <- cbind(mig_var, mig$Distance_.m., mig$Bio5_range)
mig_var$Bio11_rang <- NULL
mig_var$Bio7_range <- NULL
mig_var$SRad_mean <- NULL

mig_final <- corr.test(mig_var, method = "spearman")
corrplot(mig_final$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0), insig = "blank")

write.csv(mig_var, "env_migration.csv")
write.csv(site_var, "env_site.csv")
write.csv(s3_var_add, "env_combination.csv")