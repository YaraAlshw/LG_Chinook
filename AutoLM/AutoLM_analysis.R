#AutoLM analysis
#Originall written by Steven Micheletti (Original https://github.com/StevenMicheletti/autoLM)
#Adapted by Yara Alshwairikh


#-------Part 1: autoLM analysis------#


print("R ALERT: Checking for R dependencies. Will attempt to install automatically if not present.")

#Install packages if they don't exist already
list.of.packages <- c("data.table", "nlme","MuMIn","lme4")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) print("R ALERT: Installing dependencies for first time use....")
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

suppressMessages(require(lme4))
suppressMessages(require(nlme))
suppressMessages(require(MuMIn))
suppressMessages(require(data.table))

args <- commandArgs()

fnameo = args[6]
env.input = args[7]
freq.input = args[8]
standardize = args[9]
coords.input = args[10]


if (is.na(fnameo) ==  T){
  print("Not running via shell wrapper, variables must be included in Rscript")
  fnameo = 'autolm_results.out'
  env.input = 'env.txt'
  freq.input = 'maf.txt'#maf.txt available on Dryad (DOI to be added later)
  standardize = TRUE
  coords.input = "nope" # set to "nope" if not interested
}


autolm <- function(fnameo,
                   env.input ,
                   freq.input,
                   standardize = TRUE,
                   coords.input = NULL) {
  
  #standardized env pop variables 
  e.data <-  fread(env.input, sep='\t', header=TRUE)
  e.data <- as.data.frame(e.data)
  if (standardize == TRUE | standardize == "TRUE"){
    e.data <- data.frame(scale(e.data))
  }
  
  # Population, clust, east, north 
  if (file.exists(coords.input)) {
    c.data <- fread(coords.input, header=TRUE)
    if (ncol(c.data) < 4) stop ("Coordinate file has too many columns. 
              Make sure it is tab delimited and includes: Pop Name, genetic cluster, east (m), north (m)")
    c.data <- as.data.frame(c.data)
    colnames(c.data) = c("ID", "clust", "east", "north")
    if (min(c.data$east) < 0 | (min(c.data$north) < 0) ) stop ("Coordinates must be in meters!")   
    east  = c.data$east
    north = c.data$north
    clust = c.data$clust
    v= (as.numeric(length(e.data)))
    Var <- list()
    for (i in 1:v) {
      Var[[i]] <-e.data[, i]
    }
    names(Var) <- colnames(e.data)
    if (length(Var) < 1) stop ("No environmental variables loaded. Check input file")
  }
  
  # population, minor allele fz with row and column headings
  f.data <-  data.frame(fread(freq.input, header=TRUE), row.names=1)
  
  f.data <- as.data.frame(t(f.data))
  f.data[1,] <-as.character(f.data[1,])
  #f.data <- sapply( f.data, as.numeric )
  
  vpop = as.numeric(nrow(f.data))
  epop = as.numeric(ncol(f.data))
  env = as.numeric(ncol(e.data))
  
  Snp <- list()
  for (j in 1:epop) {
    Snp[[j]] <-f.data[, j]
  }
  
  
  names(Snp) <-colnames(f.data)
  cfile = paste0(getwd(),"/",fnameo)
  
  if (file.exists(fnameo)){
    print(" Warning: Output file exists and is being overwritten!")
  }
  
  out.headerz <-  noquote(paste("snp", "env_var", "pval", "r2", "max_freq_diff", "mean_freq", "score", sep= '\t'))
  write.table(out.headerz, file= fnameo, sep='\t', col.names = F, row.names = F, quote = F)
  
  
  if (file.exists(coords.input)) {            
    print ("Coordinates and clusters provided as random effects: Running model")
    for (i in 1:env) { 
      test <- as.vector(Var[[i]])
      ci <- names(Var[i])
      for (j in 1:vpop) { 
        tryCatch({
          cn <- row.names(f.data[j,])
          freq <- as.numeric(f.data[j,])
          #model
          mod.gau <- lme(fixed = freq ~ test , random = ~ 1 | clust, method = "ML", 
                         correlation = corGaus (1, form = ~ east + north), control = lmeControl(returnObject = TRUE))
          #Print
          me <-mean(freq)
          mxmi <-  abs((max(freq) - min(freq)))
          mo <- median(max(freq), min(freq))
          mdiff <- mean(abs(diff(freq)))
          spr = mean(as.numeric(diff(quantile(freq))))
          s1 <- summary(mod.gau)$tTable[2,5]
          s3 <- r.squaredLR(mod.gau, null.RE = TRUE)[1]
          scorez <- ((1-s1) + (mxmi) + 1-abs(0.5 - me - 0.5) + (s3*7) + (spr * 4)) / 10.5
          s4 <- cbind (cn,ci,s1,s3, mxmi, me, scorez)
          s5 <-as.data.frame(s4)
          write.table(s5, file=fnameo, sep='\t',
                      col.names = F, row.names = F, quote=FALSE, append=TRUE)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      } 
    }
  }
  
  
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  if (is.null(coords.input) | (!file.exists(coords.input))) {
    print ("No coordinates provided, there will be no random effects: Running model")
    for (i in 1:env) { 
      test <- as.vector(e.data[[i]])
      ci <- colnames(e.data[i])
      for (j in 1:vpop) { 
        tryCatch({
          cn <- row.names(f.data[j,])
          freq <- as.numeric(f.data[j,])
          #model
          mod.gau <- lm(freq ~ test, method = "qr")
          #Print
          me <-mean(freq)
          mxmi <-  abs((max(freq) - min(freq)))
          mo <- median(max(freq), min(freq))
          mdiff <- mean(abs(diff(freq)))
          spr = mean(as.numeric(diff(quantile(freq))))
          s1 <- lmp(mod.gau)
          s3 <- summary(mod.gau)$r.squared
          scorez <- ((1-s1) + (mxmi) + 1-abs(0.5 - me - 0.5) + (s3*7) + (spr * 4)) / 10.5
          s4 <- cbind (cn,ci,s1,s3, mxmi, me, scorez)
          s5 <-as.data.frame(s4)
          write.table(s5, file= fnameo, sep='\t',
                      col.names = F, row.names = F, quote=FALSE, append=TRUE)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      } 
    }
  }
}

#Run the autoLM function and save the output as an object
autolm_out <- autolm(fnameo, env.input, freq.input, standardize, coords.input)

#-------Part 2: calculate Bonferoni and BH corrections

#The following code for calcualting BH thresholds was writter by Nikos Ignatiadis and Wolfgang Huber (Source code: https://rdrr.io/bioc/IHW/src/R/helpers.R)

get_bh_threshold <- function(pvals, alpha, mtests = length(pvals)){
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m)/mtests*alpha)
  ifelse(length(prejected)==0, 0, pvals[prejected[which.max(prejected)]])
}

criticalValue_0.05 <- get_bh_threshold(autolm_out$pval, 0.05)
criticalValue_0.05 #print the threshold at alpha = 0.05


#-------Part 3: Create Manhattan plots analysis------#

#Load libraries
library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)


autolm_out <- auto_lm_out[ grep("NW", auto_lm_out$snp, invert = TRUE) , ] #remove all NW chromosomes

print("removed all NW chromosomes")

#separate the snp column into "chr_number" and "snp_position" columns
newColNames <- c("chr","BP")
newCols <- colsplit(auto_lm_out$snp, "_", newColNames)
newColNames2 <- c("junk","chr")
newCols2 <- colsplit(newCols$ch, "Chr", newColNames2)
newCols$chr <- NULL
newCols2$junk <- NULL
auto_lm_out <- cbind(auto_lm_out, newCols, newCols2) 

auto_lm_out$BP <- as.numeric(auto_lm_out$BP) #convert int to numeric
auto_lm_out$chr <- as.numeric(auto_lm_out$chr) #convert int to numeric

print("created chr column and BP column")


auto_lm_out_ggplot <- auto_lm_out %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(auto_lm_out, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate( BPcum=BP+tot)

#Then we need to prepare the X axis. Want to display the cumulative position of SNP in bp, but just show the chromosome name instead

axisdf = auto_lm_out_ggplot %>% 
  group_by(chr) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#next we split the column values to different dataframes, one for each env. variable and add names to the list

auto_lm_out_env <- auto_lm_out_ggplot %>%
  group_split(env_var)

names(auto_lm_out_env) <- c(levels(auto_lm_out$env_var))

#function to create plots and save each one as a .png
plotstuff <- function(dataset, name) {
  plot_title <- paste("Plot of ", name, sep=" ")
  ggplot(dataset, aes(x=BPcum, y=-log10(pval))) + #change pval to score if desired
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    #add in correction lines
    geom_hline(yintercept=-log10(1.187e-8), color = "darkred") + #line for Bonferroni correction, calculated manually by dividing 0.05 by the total number of SNPs (0.05/3212126 = 1.187049e-08)
    geom_hline(yintercept=-log10(criticalValue), color = "blue") + #line for BH correction as calculated in part 2
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    ggtitle(plot_title) +
    theme_classic() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  ggsave(filename = paste(plot_title,".png", sep = ""), plot = last_plot())
}

#run function, which will display each plot, but will also create an object so individual plots can be called again

mapply(plotstuff, dataset = auto_lm_out_env, name = names(auto_lm_out_env))

