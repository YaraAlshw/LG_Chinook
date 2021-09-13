#Full lfmm2 code
#Written by Yara Alshwairikh based on source code from https://rdrr.io/bioc/LEA/man/lfmm2.html
#this script needs to be run for each individual environmental variable

library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)

print("lfmm2 analysis using Bio11_site var") ##EDIT 


#Load the allele frequency matrix with 20 simulated individuals per population
rbeta_data <- readRDS("rbeta_data.rds") #this file is available on Dryad (DOI to be added later)
rbeta_data$.id <- NULL #remove .id column that gets generated when reading the data

#Load the individual environemntal csv containing 20 repeated rows pers population (total 140 rows)
Bio11_site <- read.csv("Bio11_site.csv") ###EDIT ME###

print("Loaded rbeta_data and environmental data")

# Create the lfmm2 function
setClass("lfmm2Class",
         slots = c(K = "integer", 
                   lambda = "numeric",
                   U = "matrix",
                   V = "matrix"
         )
)

lfmm2 <- function(input,
                  env, 
                  K, 
                  lambda = 1e-5){
  
  ## Check response input matrix 
  ## LEA  
  if (is.character(input)){
    Y <- read.lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## Y is an R object       
    if (is.null(input)){
      stop("NULL value for argument 'input'.")
    }
    Y <- as.matrix(input)
    Y[Y == 9] <- NA
    Y[Y == -9] <- NA
    if (anyNA(Y)) {
      stop("The input matrix contains missing values: NA, 9 or -9 not allowed.")
    }
  }
  
  ## Check independent/covariate env matrix  
  ## LEA 
  if (is.character(env)){
    X <- read.env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }
  
  if (length(K) > 1){
    stop("Multiple values of K not allowed.")
  }
  if (lambda <= 0){
    stop("The ridge regularization parameter must be positive.")
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  # run SVD of X: X = Q Sigma R
  
  svx <- svd(x = scale(X, scale = FALSE), nu = n)
  Q <- svx$u
  
  d_lambda <- c(sqrt(lambda/(lambda + svx$d)), rep(1, n-d))
  d_lambda_inv <- c(sqrt((lambda + svx$d)/lambda), rep(1, n-d))
  D_inv <- diag(d_lambda_inv)
  D  <- diag(d_lambda)
  
  # run SVD of modified Y    
  svk <- svd(D %*% t(Q) %*% scale(Y, scale = FALSE), nu = K)
  
  if (K > 1) {
    Sigma_k <- diag(svk$d[1:K])
  } else {
    Sigma_k <- as.matrix(svk$d[1])
  }
  
  # compute the latent matrix W
  W <- Q %*% D_inv %*% tcrossprod(svk$u %*% Sigma_k, svk$v[,1:K])
  
  # compute LFMM factors U and loadings V
  # Non orthogonal factors
  U <- crossprod(t(Q %*% D_inv), svk$u %*% Sigma_k)
  #U <- Q %*% D_inv %*% svk$u %*% Sigma_k
  V <- svk$v[,1:K]
  
  obj <- new("lfmm2Class")
  obj@K <- as.integer(K)
  obj@lambda <- as.numeric(lambda)
  obj@U <- as.matrix(U)
  obj@V <- as.matrix(V)
  
  ## LEA 
  return(obj)
}


setGeneric("lfmm2.test", function(object, input, env, 
                                  genomic.control = TRUE, 
                                  linear = TRUE, 
                                  family  = binomial(link = "logit")) matrix);
setMethod("lfmm2.test", "lfmm2Class",
          function(object, 
                   input,
                   env,
                   genomic.control, 
                   linear,
                   family
          ) {
            
            ## Check input matrix   
            ## LEA  
            if (is.character(input)){
              Y <- read.lfmm(input)
              lst.unique <- unique(as.numeric(Y))
              if (9 %in% lst.unique){
                stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
              }
              if (-9 %in% lst.unique){
                stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
              }
            } else {
              ## Y is an R object       
              if (is.null(input)){
                stop("NULL value for argument 'input'.")
              }
              Y <- as.matrix(input)
              Y[Y == 9] <- NA
              Y[Y == -9] <- NA
              if (anyNA(Y)) {
                stop("The input matrix contains missing values (NA or 9).")
              }
            }
            
            ## Check independent/covariate matrix  
            ## LEA 
            if (is.character(env)){
              X <- read.env(env)
              if (anyNA(X)){
                stop("'env' file contains missing data (NA).")
              }
            } else {
              if (is.null(env)){
                stop("NULL value for argument 'env'.")
              }
              X <- as.matrix(env)
              if (anyNA(X)) {
                stop("The environmental matrix contains NA.")
              }
            }
            
            d <-  ncol(X) #number of environmental variables
            n <-  nrow(X) #number of individuals
            
            if (nrow(Y) != n){
              stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
            }
            
            if (n < d) {
              stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
            }
            
            p <- ncol(Y)
            p_value <- NULL
            z_score <- NULL
            
            if (linear){
              mod_lm <- lm(Y ~ ., data = data.frame(X, object@U)) 
              sm <- summary(mod_lm)
              p_value <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 4])
              z_score <- as.matrix(sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 3]))
            } else {
              for (j in 1:p) {
                mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object@U), family = family)
                sm <- summary(mod_glm)
                p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
                z_score <- rbind(z_score, sm$coeff[2:(d + 1), 3])
              }
            }
            if (genomic.control){
              gif <- apply(z_score^2, 2, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
              p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)
            } else {
              gif <- NULL
            }
            res <- list(pval = p_value, zscores = z_score, gif = gif)
            return(res)
          }
)

print("finished creating lfmm2 function")


######################################
# Fitting an LFMM with K = 2 factors
######################################

#run the model
model <- lfmm2(input = rbeta_data, env = Bio11_site, K = 2) ##Set value of K that is appropriate for your data
print("finished running lfmm2 model")


#Code below edited on Nov 11, 2020. Changes: save the RAW output without transposing 
#to check that it's happening ok. Also print the GIF value


# Computing P-values
#save p values, z-scores, and gif into a list of 3 (pv)
pv <- lfmm2.test(object = model, input = rbeta_data, 
                 env = Bio11_site, linear = TRUE)
#df <- pv[[1]] #extarct ONLY the pval into an object

gif <- pv[[3]] #extarct and save the gif value into gif object
print("gif value")
gif #prints the gif value

pv <- as.data.frame(pv)

pv$X <- rownames(pv)
pv$gif <- NULL

write.csv(pv, file = "allraw_output_Bio11_site_K2.csv", row.names=FALSE) 

####Clean and process data for Manhattan plots

df <- pv

#rename columns into "X" for chromsoomes, and "pval"
#colnames(df) <- c('X', 'pval')

#print("renames columns into X for chromsoomes, and pval")

##calculate BH threshold

#The following code for calcualting BH thresholds was writter by Nikos Ignatiadis and Wolfgang Huber (Source code: https://rdrr.io/bioc/IHW/src/R/helpers.R)
get_bh_threshold <- function(pval, alpha, mtests = length(pval)){
  m <- length(pval)
  pval <- sort(pval)
  prejected <- which(pval <= (1:m)/mtests*alpha)
  ifelse(length(prejected)==0, 0, pval[prejected[which.max(prejected)]])
}

criticalValue <- get_bh_threshold(df$pval, 0.05)
print("BH criticla value at 0.05")
criticalValue

df <- df[ grep("NW", df$X, invert = TRUE) , ] #remove all NW chromosomes
print("removed all NW chromosomes")

#Add column env_var, and manually edit the variable name 
df['env_var']='Bio11_site' ##EDIT ME#####
print("added column env_var with manual env var name")

#separate the snp column (here, it is X) into "chr_number" and "snp_position" columns
newColNames <- c("chr","BP")
newCols <- colsplit(df$X, "_", newColNames)
newColNames2 <- c("junk","chr")
newCols2 <- colsplit(newCols$ch, "Chr", newColNames2)
newCols$chr <- NULL
newCols2$junk <- NULL
df <- cbind(df, newCols, newCols2) 

df$BP <- as.numeric(df$BP) #convert int to numeric
df$chr <- as.numeric(df$chr) #convert int to numeric

print("created chr column and BP column")

#write the fully processed data into a .csv, this format is ready for the Manhattan plot script
#write.csv(df, "final_output_Bio11_site_K2.csv") ###EDIT ME###

cand_SNP_Bio11_site_K2 <- subset(df, pval <= criticalValue)
write.csv(cand_SNP_Bio11_site_K2, "cand_SNP_Bio11_site_K2.csv", row.names=FALSE) ###EDIT ME###


df_ggplot <- df %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate( BPcum=BP+tot)


#Then we need to prepare the X axis. Want to display the cumulative position of SNP in bp, but just show the chromosome name instead

axisdf = df_ggplot %>% 
  group_by(chr) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


Bio11_site <- "Bio11_site"  ###EDIT ME###

#function to create plots and save each one as a .png
plotstuff <- function(dataset) {
  plot_title <- paste("Plot of Bio11_site, K=2") ##EDIT ME###
  ggplot(dataset, aes(x=BPcum, y=-log10(pval))) + #change pval to score if desired
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    #add corrections lines
    geom_hline(yintercept=-log10(1.187e-8), color = "darkred") + #line for Bonferroni correction, calculated manually by dividing 0.05 by the total number of SNPs (0.05/3212126 = 1.187049e-08)
    geom_hline(yintercept=-log10(criticalValue), color = "blue") + #line for BH correction as calculated in part 2
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
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


#Create the Manhattan plot 
plotstuff(dataset = df_ggplot)

#Save the plot as an object
p1 <- plotstuff(dataset = df1_ggplot)

#Save the plot as a .png file
png("p1.png")
p1
dev.off()