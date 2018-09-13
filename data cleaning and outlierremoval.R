
#load libraries used
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.us.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("ggplot2", 
              "plyr", 
              "reshape2", 
              "RColorBrewer", 
              "scales", 
              "grid",
              "VennDiagram",
              "vegan",
              "RAM",
              "tidyr",
              "Rmisc",
              "reshape","reshape2",
              "gapminder","magrittr","dplyr","ggpubr","gridExtra","
              patternplot","tibble","gplots","broom","data.table","nlme",
              "devtools","bibtex","tidyverse","multcomp","pander","xtable",
              "XLconnect","ggsignif","ggpubr","ggpmisc","packagemetrics",
              "flextable","multcompView","cowplot","lubridate","lme4")
ipak(packages)


#devtools::install_github("const-ae/ggsignif")

####Data import####

cv.data=read.csv("compiled_cv.csv")
cv.data<-cv.data %>% 
  gather(key=Date,value=OD_cm,-Substrate,-Rep)
cv.data$Date<-cv.data$Date %>% 
  gsub(pattern="X",replacement="")

metabolism.data=read.csv("compiled metabolism_longform.csv")
chl.data=read.csv("compiled_chl.csv")
field.data=read.csv("fielddata.csv")
eea.df=read.csv("Compiled_EEA.csv")
par.df=read.csv("par.csv")
usgs.df=read.csv("usgs_niles.csv")
water=read.csv("waterchem.wide.csv")




###reformat dates####
metabolism.data$Date.Sampled<-mdy(metabolism.data$Date.Sampled)
water$Date<-mdy(water$Date)
field.data$Date<-mdy(field.data$Date)
eea.df$Date.Sampled<-mdy(eea.df$Date.Sampled)
cv.data$Date.Sampled<-mdy(cv.data$Date)
par.df$Date=dmy(par.df$Date) 
usgs.df$Datetime=mdy_hm(usgs.df$timestamp)



metabolism.data$time=metabolism.data$doy-227
metabolism.data$id=seq(1:15)
eea.df$time=metabolism.data$time
cv.data$time=metabolism.data$time
combo.df=cbind(metabolism.data,eea.df[,7:12])
chl.data=droplevels(chl.data[-which(chl.data$Sample=='Cr'),])
chl.data$time=metabolism.data$time
eea.df$id=seq(1:15)
chl.data$id=seq(1:15)
cv.data$id=seq(1:15)
#s.numeric(field.data$xvelocity)





####Data cleaning-Metabolism####
metabolism.data$GPP<-ifelse(metabolism.data$GPP<0,0.000001,metabolism.data$GPP)
metabolism.data$Resp<-ifelse(metabolism.data$Resp>0,0.000001,metabolism.data$Resp)
#pc.df=cbind(metabolism.data[1],metabolism.data[3],metabolism.data[5:7],eea.df[8],eea.df[10],eea.df[12])

cv.resp.df=cbind(cv.data,metabolism.data[,5:7])
chl.data$chla_cm=286*(((abs(chl.data$X665B)-abs(chl.data$X750B))-(abs(chl.data$X665A)-abs(chl.data$X750A))*5)/chl.data$area)
chl.data$chl_abs=abs(chl.data$chla_cm)

pc.df=cbind(metabolism.data[1],metabolism.data[3],metabolism.data[5:7],eea.df[8],
            eea.df[10],eea.df[12],chl.data[13],cv.data[4],metabolism.data[8])


water.Date=list(unique(water$Date))
no3.sum=aggregate(no3~Date,water,mean)
nh4.sum=aggregate(nh4~Date,water,mean)
srp.sum=aggregate(srp~Date,water,mean)
water.sum=full_join(no3.sum,nh4.sum,by="Date")
water.sum=full_join(srp.sum,water.sum,by="Date")
pc.df=full_join(pc.df,water.sum,by=c("Date.Sampled"="Date"))
pc.df=full_join(pc.df,field.data,by=c("Date.Sampled"="Date"))

write.csv(pc.df, file="env.df")
pc.df=na.omit(pc.df)
metabolism.data=full_join(metabolism.data,water.sum,by=c("Date.Sampled"="Date"))

eea.df=full_join(eea.df,water.sum,by=c("Date.Sampled"="Date"))



####Remove outliers####
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}

outlierKD(metabolism.data,Resp)
yes
outlierKD(metabolism.data,GPP)
yes
outlierKD(metabolism.data,NEP)
yes
outlierKD(metabolism.data,Resp_OD)
no
outlierKD(metabolism.data,GPP_OD)
no
outlierKD(metabolism.data,NEP_OD)
no
outlierKD(eea.df,BGT2)
yes
outlierKD(eea.df,NAGT2)
yes
outlierKD(eea.df,PT2)
yes
outlierKD(eea.df,BGT2_OD)
no
outlierKD(eea.df,NAGT2_OD)
no
outlierKD(eea.df,PT2_OD)
no
outlierKD(eea.df,srp)
yes
outlierKD(eea.df,no3)
no
outlierKD(eea.df,nh4)
no
outlierKD(cv.data,OD_cm)
yes
outlierKD(chl.data,chl_abs_OD)
no


save.image(file="assays.Rdata")
