#Simulating genotypes using the rbeta function 
#Written by Alexis Garretson and Yara Alshwairikh 
#Last edited July 3rd, 2020

library("tidyverse")
library("dplyr")
library("miceadds")

load.Rdata(filename = "freqs.RData", "maf") #load the allele frequency matrix

# sample size n for each population
#1	lower Yakima R fall-run	46
#2	Wenatchee R summer-run	61
#3	Lyons Ferry weir fall-run	92
#4	Methow R summer-run	68
#5	upper Deschutes River fall-run	48
#6	Clearwater River fall-run	96
#7	Priest Rapids fall-run	46

maf2 <- as.data.frame(maf) #convert maf into a data frame for slice() to work

#Subset the maf2 dataframe based on population
Yakima <- maf2 %>% slice(1) 
Wenatchee <- maf2 %>% slice(2)
LyonsFerry <- maf2 %>% slice(3) 
Methow <- maf2 %>% slice(4) 
Deschutes <- maf2 %>% slice(5) 
Clearwater <- maf2 %>% slice(6) 
Priest <- maf2 %>% slice(7) 

#run every function separatly on each population

#-----#Yakima, n=46
number_of_columns <- ncol(Yakima)
number_of_rows <- nrow(Yakima)
output_lfmm_Yakima <- data.frame(rep(NA, 46))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(46, 
                           shape1= ((46 * as.numeric(Yakima[1, column])) + 1), 
                           shape2= (46 * (1 - as.numeric(Yakima[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 46)
    output_lfmm_Yakima <- cbind(output_lfmm_Yakima, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_Yakima)[column] <- names(Yakima)[column]
    next
  }
  output_lfmm_Yakima <- cbind(output_lfmm_Yakima, rbeta_value)##put the value in the right place
  names(output_lfmm_Yakima)[column] <- names(Yakima)[column]
} 
output_lfmm_Yakima <- output_lfmm_Yakima[,-1] 

#add a new column to indicate sample number
output_lfmm_Yakima <- output_lfmm_Yakima %>% mutate(sample = c(1:46)) 
write.csv(output_lfmm_Yakima, "output_lfmm_Yakima.csv") #write and save the data into .csv file

#-----#Wenatchee, n=61
number_of_columns <- ncol(Wenatchee)
number_of_rows <- nrow(Wenatchee)
output_lfmm_Wenatchee <- data.frame(rep(NA, 61))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(61, 
                           shape1= ((61 * as.numeric(Wenatchee[1, column])) + 1), 
                           shape2= (61 * (1 - as.numeric(Wenatchee[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 61)
    output_lfmm_Wenatchee <- cbind(output_lfmm_Wenatchee, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_Wenatchee)[column] <- names(Wenatchee)[column]
    next
  }
  output_lfmm_Wenatchee <- cbind(output_lfmm_Wenatchee, rbeta_value)##put the value in the right place
  names(output_lfmm_Wenatchee)[column] <- names(Wenatchee)[column]
} 
output_lfmm_Wenatchee <- output_lfmm_Wenatchee[,-1]

#add a new column to indicate sample number
output_lfmm_Wenatchee <- output_lfmm_Wenatchee %>% mutate(sample = c(47:107)) 


#-----#LyonsFerry, n=92
number_of_columns <- ncol(LyonsFerry)
number_of_rows <- nrow(LyonsFerry)
output_lfmm_LyonsFerry <- data.frame(rep(NA, 92))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(92, 
                           shape1= ((92 * as.numeric(LyonsFerry[1, column])) + 1), 
                           shape2= (92 * (1 - as.numeric(LyonsFerry[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 92)
    output_lfmm_LyonsFerry <- cbind(output_lfmm_LyonsFerry, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_LyonsFerry)[column] <- names(LyonsFerry)[column]
    next
  }
  output_lfmm_LyonsFerry <- cbind(output_lfmm_LyonsFerry, rbeta_value)##put the value in the right place
  names(output_lfmm_LyonsFerry)[column] <- names(LyonsFerry)[column]
} 
output_lfmm_LyonsFerry <- output_lfmm_LyonsFerry[,-1]

#add a new column to indicate sample number
output_lfmm_LyonsFerry <- output_lfmm_LyonsFerry %>% mutate(sample = c(108:199)) 
write.csv(output_lfmm_LyonsFerry, "output_lfmm_LyonsFerry.csv") #write and save the data into .csv file

#-----#Methow, n=68
number_of_columns <- ncol(Methow)
number_of_rows <- nrow(Methow)
output_lfmm_Methow <- data.frame(rep(NA, 68))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(68, 
                           shape1= ((68 * as.numeric(Methow[1, column])) + 1), 
                           shape2= (68 * (1 - as.numeric(Methow[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 68)
    output_lfmm_Methow <- cbind(output_lfmm_Methow, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_Methow)[column] <- names(Methow)[column]
    next
  }
  output_lfmm_Methow <- cbind(output_lfmm_Methow, rbeta_value)##put the value in the right place
  names(output_lfmm_Methow)[column] <- names(Methow)[column]
} 
output_lfmm_Methow <- output_lfmm_Methow[,-1]
#add a new column to indicate sample number
output_lfmm_Methow <- output_lfmm_Methow %>% mutate(sample = c(200:267)) 
write.csv(output_lfmm_Methow, "output_lfmm_Methow.csv") #write and save the data into .csv file

#-----#Deschutes, n=48
number_of_columns <- ncol(Deschutes)
number_of_rows <- nrow(Deschutes)
output_lfmm_Deschutes <- data.frame(rep(NA, 48))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(48, 
                           shape1= ((48 * as.numeric(Deschutes[1, column])) + 1), 
                           shape2= (48 * (1 - as.numeric(Deschutes[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 48)
    output_lfmm_Deschutes <- cbind(output_lfmm_Deschutes, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_Deschutes)[column] <- names(Deschutes)[column]
    next
  }
  output_lfmm_Deschutes <- cbind(output_lfmm_Deschutes, rbeta_value)##put the value in the right place
  names(output_lfmm_Deschutes)[column] <- names(Deschutes)[column]
} 
output_lfmm_Deschutes <- output_lfmm_Deschutes[,-1]

#add a new column to indicate sample number
output_lfmm_Deschutes <- output_lfmm_Deschutes %>% mutate(sample = c(268:315)) 
write.csv(output_lfmm_Deschutes, "output_lfmm_Deschutes.csv") #write and save the data into .csv file

#-----#Clearwater, n=96
number_of_columns <- ncol(Clearwater)
number_of_rows <- nrow(Clearwater)
output_lfmm_Clearwater <- data.frame(rep(NA, 96))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(96, 
                           shape1= ((96 * as.numeric(Clearwater[1, column])) + 1), 
                           shape2= (96 * (1 - as.numeric(Clearwater[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 96)
    output_lfmm_Clearwater <- cbind(output_lfmm_Clearwater, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_Clearwater)[column] <- names(Clearwater)[column]
    next
  }
  output_lfmm_Clearwater <- cbind(output_lfmm_Clearwater, rbeta_value)##put the value in the right place
  names(output_lfmm_Clearwater)[column] <- names(Clearwater)[column]
} 
output_lfmm_Clearwater <- output_lfmm_Clearwater[,-1]

#add a new column to indicate sample number
output_lfmm_Clearwater <- output_lfmm_Clearwater %>% mutate(sample = c(316:411)) 
write.csv(output_lfmm_Clearwater, "output_lfmm_Clearwater.csv") #write and save the data into .csv file

#-----#Priest, n=46
number_of_columns <- ncol(Priest)
number_of_rows <- nrow(Priest)
output_lfmm_Priest <- data.frame(rep(NA, 46))
for (column in 1:number_of_columns){          #Go through all columns
  rbeta_value <- try(rbeta(46, 
                           shape1= ((46 * as.numeric(Priest[1, column])) + 1), 
                           shape2= (46 * (1 - as.numeric(Priest[1, column])) + 1)))
  ##get the rbeta value
  if(inherits(rbeta_value, "try-error")){
    rbeta_value <-rep(NA, 46)
    output_lfmm_Priest <- cbind(output_lfmm_Priest, rbeta_value) ##If the rbeta function fails, for some reason, make the value NA
    names(output_lfmm_Priest)[column] <- names(Priest)[column]
    next
  }
  output_lfmm_Priest <- cbind(output_lfmm_Priest, rbeta_value) ##put the value in the right place
  names(output_lfmm_Priest)[column] <- names(Priest)[column]
} 
output_lfmm_Priest <- output_lfmm_Priest[,-1]

#add a new column to indicate sample number
output_lfmm_Priest <- output_lfmm_Priest %>% mutate(sample = c(412:457)) 
write.csv(output_lfmm_Priest, "output_lfmm_Priest.csv") #write and save the data into .csv file

#bind all data into one master file
all_genotypes <- rbind(output_lfmm_Yakima, output_lfmm_Wenatchee, output_lfmm_LyonsFerry,
                       output_lfmm_Methow, output_lfmm_Deschutes, output_lfmm_Clearwater, output_lfmm_Priest)


write.csv(all_genotypes, "all_genotypes.csv") #write and save the data into .csv file


