#Code for making figures of DayMet data for Alaska and Georgia
#August 29, 2021
#Written by Yara Alshwairikh

TEST <- transform(alaska_daymet, Min = pmin(prcp), 
                  Max = pmax(prcp))


library(ggplot2)
ggplot(TEST) +
  geom_line(aes(yday, prcp), group = 1) +
  geom_ribbon(aes(x = prcp, ymax = Max, ymin = Min), alpha = 0.6, fill = "skyblue")


plotLineWithRange <- function(x, yVal, yMin, yMax,
                              lineColor="Black", rangeColor="LightBlue",
                              main="", xlab="X", ylab="Y"){
  if(missing(x)){
    x <- 1:length(yVal)
  }
  stopifnot(length(yVal) == length(yMin) && length(yVal) == length(yMax))
  
  plot(x=c(min(x),max(x)),y=c(min(yMin),max(yMax)),type="n", main=main,xlab=xlab,ylab=ylab)
  polygon(x=c(x,rev(x)),y=c(yMax,rev(yVal)),col=rangeColor,border=NA)
  polygon(x=c(x,rev(x)),y=c(yMin,rev(yVal)),col=rangeColor,border=NA)
  lines(x=x,y=yVal,col=lineColor)
}


# usage example:
d <- plotLineWithRange(alaska_daymet$yday, 
                  yVal=alaska_daymet$prcp,yMin=alaska_daymet$prcp,
                  yMax=alaska_daymet$prcp,main="Average")

a <- ggplot(alaska_daymet, aes(yday, prcp))

a + geom_ribbon(aes(ymin = prcp, ymax = prcp)) - x, 
ymax, ymin, alpha = 0.5, color = blue, fill = "skyblue")

###

df2 <- alaska_daymet %>% group_by(yday) %>%   # the grouping variable
  summarise(mean_min = mean(tmin), min = min(tmin), max = max(tmin))  # calculates the mean of each group
df2$state <- "Alaska" #add state column

df3 <- georgia_daymet %>% group_by(yday) %>%   # the grouping variable
  summarise(mean_min = mean(tmin), min = min(tmin), max = max(tmin))  # calculates the mean of each group
df3$state <- "Georgia" #add state column

dd <- bind_rows(df2, df3) #bind into one dataframe



##THIS WORKS
minplot <- ggplot(data=dd, aes(y=mean_min, x=yday, group=state, colour=state,
                    fill=state)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=min, ymax=max), alpha=.3, linetype=0) +
  theme_bw()+
  labs(x = "x label",
       y = "y label",
       title = "Title")


df4 <- alaska_daymet %>% group_by(yday) %>%   # the grouping variable
  summarise(mean_max = mean(tmax), min = min(tmax), max = max(tmax))  # calculates the mean of each group
df4$state <- "Alaska" #add state column

df5 <- georgia_daymet %>% group_by(yday) %>%   # the grouping variable
  summarise(mean_max = mean(tmax), min = min(tmax), max = max(tmax))  # calculates the mean of each group
df5$state <- "Georgia"
ff <- bind_rows(df4, df5)

maxplot <- ggplot(data=ff, aes(y=mean_max, x=yday, group=state, colour=state,
                    fill=state)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=min, ymax=max), alpha=.3, linetype=0) 
  theme_bw()+
  labs(x = "x label",
       y = "y label",
       title = "Title") 
  
  #+ scale_fill_manual(values=c("blue", "red"))
  
  
  require(gridExtra)
  plot1 <- qplot(1)
  plot2 <- qplot(1)
  grid.arrange(maxplot, minplot, ncol=2)

  
  par(mfrow = c(1,2))

  maxplot
minplot  
    
  