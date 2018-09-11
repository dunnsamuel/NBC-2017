#rm(list = ls(all = TRUE))



####Data import####


cv.data=read.csv("compiled_cv.csv")
metabolism.data=read.csv("compiled metabolism_longform.csv")
chl.data=read.csv("compiled_chl.csv")
field.data=read.csv("fielddata.csv")
eea.df=read.csv("Compiled_EEA.csv")
par.df=read.csv("par.csv")
usgs.df=read.csv("usgs_niles.csv")
water=read.csv("waterchem.wide.csv")

###reformat dates####
metabolism.data$Date.Sampled=strptime(as.character(metabolism.data$Date.Sampled),"%m/%d/%Y")
metabolism.data$Date.Sampled=format(metabolism.data$Date.Sampled, "%m/%d/%y")
water$Date=strptime(as.character(water$Date),"%m/%d/%Y")
water$Date=format(water$Date,"%m/%d/%y")
field.data$Date=strptime(as.character(field.data$Date),"%m/%d/%Y")
field.data$Date=format(field.data$Date,"%m/%d/%y")
#field.data$sDate=format(field.data$Date,"%m/%d/%")
eea.df$Date.Sampled=strptime(as.character(eea.df$Date.Sampled),"%m/%d/%Y")
eea.df$Date.Sampled=format(eea.df$Date.Sampled, "%m/%d/%y")
cv.data$Date.Sampled=strptime(as.character(cv.data$Date.Sampled),"%m/%d/%Y")
cv.data$Date.Sampled=format(cv.data$Date.Sampled, "%m/%d/%y")
par.df$Date=strptime(as.character(par.df$Date),"%d/%m/%Y")
par.df$Date=format(par.df$Date, "%m/%d/%y")
usgs.df$Date=as.character(lapply(strsplit(as.character(usgs.df$timestamp), split=" "), "[", 1))
usgs.df$time=as.character(lapply(strsplit(as.character(usgs.df$timestamp), split=" "), "[", 2))
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
cv.data$OD_cm=cv.data$OD/cv.data$area
cv.resp.df=cbind(cv.data,metabolism.data[,5:7])
chl.data$chla_cm=286*(((abs(chl.data$X665B)-abs(chl.data$X750B))-(abs(chl.data$X665A)-abs(chl.data$X750A))*5)/chl.data$area)
chl.data$chl_abs=abs(chl.data$chla_cm)

pc.df=cbind(metabolism.data[1],metabolism.data[3],metabolism.data[5:7],eea.df[8],
            eea.df[10],eea.df[12],chl.data[13],cv.data[9],metabolism.data[8])


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

#normalize per OD
metabolism.data$GPP_OD=metabolism.data$GPP/cv.data$OD_cm
metabolism.data$NEP_OD=metabolism.data$NEP/cv.data$OD_cm
metabolism.data$Resp_OD=metabolism.data$Resp/cv.data$OD_cm
chl.data$chl_abs_OD=chl.data$chl_abs/cv.data$OD_cm
eea.df$BGT2_OD=eea.df$BGT2/cv.data$OD_cm
eea.df$NAGT2_OD=eea.df$NAGT2/cv.data$OD_cm
eea.df$PT2_OD=eea.df$PT2/cv.data$OD_cm


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
y
outlierKD(metabolism.data,GPP)
y
outlierKD(metabolism.data,NEP)
y
outlierKD(metabolism.data,Resp_OD)
n
outlierKD(metabolism.data,GPP_OD)
n
outlierKD(metabolism.data,NEP_OD)
n
outlierKD(eea.df,BGT2)
y
outlierKD(eea.df,NAGT2)
y
outlierKD(eea.df,PT2)
y
outlierKD(eea.df,BGT2_OD)
n#
outlierKD(eea.df,NAGT2_OD)
n#
outlierKD(eea.df,PT2_OD)
n#
outlierKD(eea.df,srp)
y
outlierKD(eea.df,no3)
n
outlierKD(eea.df,nh4)
n
outlierKD(cv.data,OD_cm)
y
outlierKD(chl.data,chl_abs_OD)
n

#rename substrates
metabolism.data$Sample=factor(metabolism.data$Sample,levels=c("H","So","Sh","F","T"),labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))
eea.df$Sample=factor(eea.df$Sample,levels=c("H","So","Sh","F","T"),labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))
cv.data$Substrate=factor(cv.data$Substrate,levels=c("H","So","Sh","F","T"),labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))
chl.data$Sample=factor(chl.data$Sample,levels=c("H","So","Sh","F","T"),labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))
pc.df$Sample=factor(pc.df$Sample,levels=c("H","So","Sh","F","T"),labels=c("PVC","Soft PE","Sheet PE","PS","Tile"))

###

library(extrafont)
#font_import()
y

#### Global Fonts####
theme_min = function (size=14, font=NA, face='plain', 
                      panelColor=backgroundColor, axisColor='#999999', 
                      gridColor=gridLinesColor, textColor='black') 
{
  theme_text = function(...)
    ggplot2::theme_text(family=font, face=face, colour=textColor, 
                        size=size, ...)
  
  opts(
    axis.text.x = theme_text(),
    axis.text.y = theme_text(),
    axis.line = theme_blank(),
    axis.ticks = theme_segment(colour=axisColor, size=0.25),
    panel.border = theme_rect(colour=backgroundColor),
    legend.background = theme_blank(),
    legend.key = theme_blank(),
    legend.key.size = unit(1.5, 'lines'),
    legend.text = theme_text(hjust=0),
    legend.title = theme_text(hjust=0),
    panel.background = theme_rect(fill=panelColor, colour=NA),
    panel.grid.major = theme_line(colour=gridColor, size=0.33),
    panel.grid.minor = theme_blank(),
    strip.background = theme_rect(fill=NA, colour=NA),
    strip.text.x = theme_text(hjust=0),
    strip.text.y = theme_text(angle=-90),
    plot.title = theme_text(hjust=0),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'lines'))
}

##Create a custom font type. Could be 'F', 'TEST', whatever
#windowsFonts(F = windowsFont('Wide Latin'))



####integrate PAR####
library(plyr)
int_par=aggregate(par.df$microMole.s.m2,by=list((substr(par.df$Date,1,385))),mean)
int_par=rename(int_par,c("Group.1"="Date.Sampled","x"="Summed.PAR"))


###Time Series rmANOVAS####
####two-way anovas, no error terms####
require(nlme)
require(multcomp)
require(knitr)

####Resp####
resp.aov=aov(Resp~(Sample*Date.Sampled),metabolism.data)
resp.sum=summary(resp.aov)
tuk.resp.sample=TukeyHSD(resp.aov,"Sample")
tuk.resp.date=TukeyHSD(resp.aov,"Date.Sampled")
par(mfrow=c(2,2))
plot(resp.aov)#> in residual plot, falls of 0 9due to negative values)45,1110,1218 are outliers to be removed
resp.results=list(resp.sum,tuk.resp.sample,tuk.resp.date)
capture.output(resp.results,file="resp.twoway.doc")

####GPP####
GPP.aov=aov(GPP~Sample*Date.Sampled,metabolism.data)
gpp.sum=summary(GPP.aov)
tuk.gpp.sample=TukeyHSD(GPP.aov,"Sample")
tuk.gpp.date=TukeyHSD(GPP.aov,"Date.Sampled")
par(mfrow=c(2,2))
plot(GPP.aov)#< in residual plot, centered on zero
gpp.results=list(gpp.sum,tuk.gpp.sample,tuk.gpp.date)
capture.output(gpp.results,file="gpp.twoway.doc")

####NEP####
NEP.aov=aov(NEP~Sample*Date.Sampled,metabolism.data)
nep.sum=summary(NEP.aov)
tuk.nep.sample=TukeyHSD(NEP.aov,"Sample")
tuk.nep.date=TukeyHSD(NEP.aov,"Date.Sampled")
par(mfrow=c(2,2))
plot(NEP.aov) #no problem
nep.results=list(nep.sum,tuk.nep.sample,tuk.nep.date)
capture.output(nep.results,file="nep.twoway.doc")

####2-way RM anova via LME model approach#####
library(pander)
library(xtable)
require(XLConnect)
require(Rmisc)
require(nlme)
library(broom)
panderOptions('digits',4)
panderOptions('round',4)
panderOptions('keep.trailing.zeros',TRUE)
#Resp
lme.resp=anova(lme(Resp~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=metabolism.data))
lme.resp
summary(lme.resp)
#summary(lme.resp)

#GPP
lme.gpp=anova(lme(GPP~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=metabolism.data))
lme.gpp

#NEP
lme.nep=anova(lme(NEP~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=metabolism.data))
lme.nep

#BG
lme.bg=anova(lme(BGT2~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=eea.df))
lme.bg

#NAG
lme.nag=anova(lme(NAGT2~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=eea.df))
lme.nag

#P
lme.p=anova(lme(PT2~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=eea.df))
lme.p

#Chl
lme.chl=anova(lme(chl_abs~Date.Sampled*Sample,random=list(id=pdBlocked(list(~1,pdIdent(~Sample-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=chl.data))
lme.chl
l=xtable(lme.chl)
print(l)
#CV
require(lme4)
cv.data.drop=droplevels(cv.data[-which(cv.data$Date.Sampled=='08/15/17'),])
lme.cv=anova(lmer(OD_cm~Date.Sampled*Substrate+(Date.Sampled|Substrate),na.action=na.omit,data=cv.data.drop))
lme.cv
#model fails to converge, basically there is a date effect but no effect of sample
lme.cv=anova(lm(OD_cm~Date.Sampled*Substrate,random=list(id=pdBlocked(list(~1,pdIdent(~Substrate-1),pdIdent(~Date.Sampled-1)))),method="ML",na.action=na.omit,data=cv.data))
lme.cv
cv.anova=aov(data=cv.data, OD_cm~Date.Sampled)
TukeyHSD(cv.anova)
TUKEY.cv=TukeyHSD(cv.anova,conf.level=.95,'Date.Sampled',na.rm=T)
generate_label_df<-function(TUKEY.cv,Date.Sampled){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY.cv[[Date.Sampled]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
LABELS.cv=generate_label_df(TUKEY.cv,"Date.Sampled")

####Table 1####
writeWorksheetToFile("lmeoutput.xlsx",
                     data=list(lme.resp,lme.gpp,lme.nep,lme.bg,lme.nag,lme.p,lme.chl,lme.cv),
                     sheet=c("resp","gpp","nep","bg","nag","p","chl","cv"),
                     header=T,clearSheets = T,rownames="var")

table1=list(lme.resp,lme.gpp,lme.nep,lme.bg,lme.nag,lme.p,lme.chl)

names(table1)=c("resp","gpp","nep","bg","nag","p","chl")
table1.df=do.call(rbind.data.frame,table1)
table1.df=cbind(rownames(table1.df),table1.df)
rownames(table1.df)=NULL
colnames(table1.df)[1]="assay"
table1.df%>%separate(table1.df,assay,into=c("Assay","Factor"),sep=".")

#table1.assay=data.frame(do.call('rbind',strsplit(as.character(table1.df$assay),'.',fixed=TRUE)))


###one way anovas by factor level Time####
r.o=lapply(split(x=metabolism.data,f=metabolism.data$Date.Sampled),aov,formula=Resp~Sample)
r.o.sum=lapply(r.o,summary)
g.o=lapply(split(x=metabolism.data,f=metabolism.data$Date.Sampled),aov,formula=GPP~Sample)
g.o.sum=lapply(g.o,summary)
n.o=lapply(split(x=metabolism.data,f=metabolism.data$Date.Sampled),aov,formula=NEP~Sample)
n.o.sum=lapply(n.o,summary)
oneway.time=list(r.o.sum,g.o.sum,n.o.sum)
capture.output(oneway.time,file="oneway.time.doc")
bt.o=lapply(split(x=eea.df,f=eea.df$Date.Sampled),aov,formula=BGT2~Sample)
bt.o.sum=lapply(bt.o,summary) #p correction a=0.005, non sig with correction a


####Tukey for loop, TUKEY HSD on date by date basis####
require(multcomp)
metabolism.data$Date.Sampled=as.factor(metabolism.data$Date.Sampled)
####Resp####
dates.resp=lapply(levels(metabolism.data$Date.Sampled), function(i){
  dat=metabolism.data[metabolism.data$Date.Sampled==i,] #subset for date i
  part1<-with(dat, aov(Resp~Sample))
  print(part1)
  list(
    part1=part1,
    part2=TukeyHSD(x=part1,which="Sample", conf.level=0.95)
  )
})
names(dates.resp)<-levels(metabolism.data$Date.Sampled)
dates.resp
capture.output(dates.resp,file="tuk.resp.bydate.doc")
####NEP####
dates.nep=lapply(levels(metabolism.data$Date.Sampled), function(i){
  dat=metabolism.data[metabolism.data$Date.Sampled==i,] #subset for date i
  part1<-with(dat, aov(NEP~Sample))
  print(part1)
  list(
    part1=part1,
    part2=TukeyHSD(x=part1,which="Sample", conf.level=0.95)
  )
})
names(dates.nep)<-levels(metabolism.data$Date.Sampled)
dates.nep
capture.output(dates.nep,file="tuk.nep.bydate.doc")
####GPP####

dates.gpp=lapply(levels(metabolism.data$Date.Sampled), function(i){
  dat=metabolism.data[metabolism.data$Date.Sampled==i,] #subset for date i
  part1<-with(dat, aov(GPP~Sample))
  print(part1)
  list(
    part1=part1,
    part2=TukeyHSD(x=part1,which="Sample", conf.level=0.95)
  )
})
names(dates.gpp)<-levels(metabolism.data$Date.Sampled)
dates.gpp
capture.output(dates.gpp,file="tuk.gpp.bydate.doc")

eea.df$Date.Sampled=as.factor(eea.df$Date.Sampled)
###BG#### 
dates.bg=lapply(levels(eea.df$Date.Sampled), function(i){
  dat=eea.df[eea.df$Date.Sampled==i,] #subset for date i
  part1<-with(dat, aov(BGT2~Sample))
  print(summary(part1))
  list(
    part1=part1,
    part2=TukeyHSD(x=part1,which="Sample", conf.level=0.95)
  )
})
names(dates.bg)<-levels(eea.df$Date.Sampled)
dates.bg
capture.output(dates.bg,file="tuk.bg.bydate.doc")

cv.data$Date.Sampled=as.factor(cv.data$Date.Sampled)
###CV#### 
dates.cv=lapply(levels(cv.data$Date.Sampled), function(i){
  dat=cv.data[cv.data$Date.Sampled==i,] #subset for date i
  part1<-with(dat, aov(OD_cm~Substrate))
  print(summary(part1))
  list(
    part1=part1,
    part2=TukeyHSD(x=part1,which="Substrate", conf.level=0.95)
  )
})
names(dates.cv)<-levels(cv.data$Date.Sampled)
dates.cv
capture.output(dates.cv,file="tuk.cv.bydate.doc")

####metabolism first cut####
library(ggplot2)
library(gridExtra)

tabletheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = .5)),
  rowhead = list(fg_params=list(cex = .5)))


library(pander)
panderOptions('digits',4)
panderOptions('round',4)
panderOptions('keep.trailing.zeros',TRUE)

raw.plot.g=ggplot(data=metabolism.data,aes(x=time,y=GPP,shape=Sample,group=Sample))+geom_point()
raw.plot.g+geom_smooth(inherit.aes=TRUE)+
  theme(panel.background=element_blank())+
  ylab(bquote('GPP (mg'~cm^-2~h^-1~')'))+xlab("Days of Incubation")+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.position=c(0.15,0.85),
        legend.title = element_blank())+
  annotation_custom(tableGrob(round(as.matrix(lme.gpp),digits=4),theme=tabletheme),xmin=0,xmax=20,ymin=0.75)
#ggsave("GPP.jpeg",device="jpeg", width=6.5,height=6.2, units="in")


raw.plot.r=ggplot(data=metabolism.data,aes(x=time,y=Resp,color=Sample,group=Sample))+geom_point()
raw.plot.r+geom_smooth(inherit.aes=TRUE)+
  theme(panel.background=element_blank())+
  ylab(bquote('Respiration ('*mu~'g'~cm^-2~h^-1~')'))+xlab("Day of Incubation")+
  theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
  annotation_custom(tableGrob(round(as.matrix(lme.resp),digits=4),theme=tabletheme),xmin=0,xmax=20,ymax=-0.8)
#ggsave("Resp.jpeg",device="jpeg", width=6.5,height=6.2, units="in")

raw.plot.n=ggplot(data=metabolism.data,aes(x=time,y=NEP,group=Sample,color=Sample))+geom_point()
raw.plot.n+geom_smooth(inherit.aes=TRUE)+
  theme(panel.background=element_blank())+
  geom_hline(yintercept=0)+
  ylab(bquote('NEP ('*mu~'g'~cm^-2~h^-1~')'))+xlab("Day of Incubation")+
  theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
  annotation_custom(tableGrob(round(as.matrix(lme.nep),digits=4),theme=tabletheme),xmin=20,xmax=30,ymax=-0.6)

#ggsave("NEP.jpeg",device="jpeg", width=6.5,height=6.2, units="in")

####crystal violet first cut####
devtools::install_github("const-ae/ggsignif")
library(ggsignif)
library(ggplot2)
cv.data=na.omit(cv.data)
theme_set(theme_gray(base_size = 18))
cv.sum=summarySE(cv.data,measurevar="OD_cm",groupvars=c("time","Substrate"))
substrate.plot.cv=ggplot(cv.sum,aes(time,group=Substrate,na.rm=T,y=OD_cm))+
  geom_point(data=cv.sum,aes(y=OD_cm,shape=Substrate),stat="identity",color="black",size=3)+
  geom_line(aes(y=OD_cm,linetype=Substrate),width=.5)+
  theme_classic()+
  geom_errorbar(aes(ymin=OD_cm-se,ymax=OD_cm+se,fill=Substrate),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote('Optical Density  ('~cm^-2~')'))+
  theme(panel.background=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.15,0.85),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14))+
  xlab("Days of Incubation")+
  scale_x_continuous(breaks=seq(0,30,5))
substrate.plot.cv
ggsave("cv_od.jpeg",device="jpeg", width=10.52 ,height=6.5, units="in")





####chla first cu####
require(ggplot2)
require(gridExtra)
require(ggpubr)

#chl.data$chla_abs=abs(chl.data$chla_cm)
chl.data=na.omit(chl.data)
chl.sum=summarySE(chl.data,measurevar="chl_abs",groupvars=c("time","Sample"))
chl.sum.od=summarySE(chl.data,measurevar="chl_abs_OD",groupvars=c("time","Sample"))
raw.plot.chlv=ggplot(chl.sum,aes(time,group=Sample,na.rm=T))+
  geom_point(data=chl.sum,aes(y=chl_abs, shape=Sample),stat="identity",color="black",size=3)+
  geom_line(aes(y=chl_abs, linetype=Sample),width=.5)+theme_classic()+
  geom_errorbar(aes(ymin=chl_abs-se,ymax=chl_abs+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote('Chlorohyll-a ( '*mu~'g'~cm^-2~')'))+xlab("Days of Incubation")+
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.position=c(0.15,0.85),
        legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14))+
  scale_y_continuous(breaks=c(0,1,2,3))+scale_x_continuous(breaks=seq(0,30,5))
raw.plot.chlv
#ggarrange(substrate.plot.cv+theme(axis.text.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank()),
#raw.plot.chlv+theme(axis.text.x=element_text(angle=35,hjust=1)),common.legend = T,nrow=2,legend="right",align="v")
ggsave("chl.jpeg",device="jpeg", width=6.5,height=6.5, units="in")
#ggsave(ggarrange(substrate.plot.cv+theme(axis.text.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x = element_blank()),
#raw.plot.chlv+theme(axis.text.x=element_text(angle=35,hjust=1)),common.legend = T,nrow=2,legend="right",align="v"),
#device="jpeg",file="cv_od.jpeg",width=6.5,height=6.2, units="in")

#chl_od

raw.plot.chl.od=ggplot(chl.sum.od,aes(time,group=Sample,na.rm=T))+
  geom_point(data=chl.sum.od,aes(y=chl_abs_OD, fill=Sample),stat="identity",color="black",size=3,pch=21)+
  geom_line(aes(y=chl_abs_OD, color=Sample),width=.5)+theme_classic()+
  geom_errorbar(aes(ymin=chl_abs_OD-se,ymax=chl_abs_OD+se,fill=Sample,color=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote('Chlorohyll-a ( '*mu~'g'~cm^-2~')'))+xlab("Days of Incubation")+
  theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+theme_classic()
raw.plot.chl.od
ggsave("chl.od.jpeg",device="jpeg", width=6.5,height=6.2, units="in")

####Tukey plots####
####split df by date####
metabolism.data$dates=as.Date(metabolism.data$Date.Sampled, "%m/%d/%Y")
metabolism.data$dates=as.factor(metabolism.data$Date.Sampled)
metabolism.data$Date.Sampled2=as.character(metabolism.data$Date.Sampled)
list.df<-split(metabolism.data, metabolism.data$dates)
list2env(list.df,envir=.GlobalEnv)

#for(i in 1:length(list.df)){
#  assign(new.names[i],list.df[i])
#}
####ADDD IN SIG LETTERS?####
#####plots####
require(ggplot2)
tukey.plots.gpp<-function(metabolism.data,na.rm=TRUE,...){
  dates_list=unique(metabolism.data$dates)
  for (i in seq_along(dates_list)){
    plot<-ggplot(metabolism.data[metabolism.data$dates==levels(metabolism.data$dates)[i],],
                 aes(x=Sample,y=GPP,color=Sample))+geom_boxplot()+ggtitle(paste(dates_list[i]))
    print(plot)
  }
}

#tukey.plots.gpp(metabolism.data)


tukey.plots.nep<-function(metabolism.data,na.rm=TRUE,...){
  dates_list=unique(metabolism.data$dates)
  for (i in seq_along(dates_list)){
    plot<-ggplot(metabolism.data[metabolism.data$dates==levels(metabolism.data$dates)[i],],
                 aes(x=Sample,y=NEP,color=Sample))+geom_boxplot()+ggtitle(paste(dates_list[i]))
    print(plot)
  }
}

#tukey.plots.nep(metabolism.data)

tukey.plots.Resp<-function(metabolism.data,na.rm=TRUE,...){
  dates_list=unique(metabolism.data$dates)
  for (i in seq_along(dates_list)){
    plot<-ggplot(metabolism.data[metabolism.data$dates==levels(metabolism.data$dates)[i],],
                 aes(x=Sample,y=Resp,color=Sample))+geom_boxplot()+ggtitle(paste(dates_list[i]))
    print(plot)
  }
}

#tukey.plots.Resp(metabolism.data)


####EEA ####
require(ggplot2)
require(ggpubr)
eea.df$dBG=eea.df$BGT2-eea.df$BGT1
eea.df$dNAG=eea.df$NAGT2-eea.df$NAGT1
eea.df$dP=eea.df$PT2-eea.df$PT1
eea.df=na.omit(eea.df)
eea.df=droplevels(eea.df[-which(eea.df$Date.Sampled=='09/05/17'),])  ####PLATE DROPPED THIS DATE!!!!!!!

bg.aov=aov(BGT2~Date.Sampled,data=eea.df)
TUKEY.bg=TukeyHSD(bg.aov,conf.level=.95,'Date.Sampled',na.rm=T)
plot(TUKEY.bg,las=1,col="brown")
generate_label_df<-function(TUKEY.bg,Date.Sampled){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY.bg[[Date.Sampled]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
LABELS.BG=generate_label_df(TUKEY.bg,"Date.Sampled")
capture.output(TUKEY.bg,file="tukeybg.docx")


nag.aov=aov(NAGT2~Date.Sampled,data=eea.df)
TUKEY.nag=TukeyHSD(nag.aov,conf.level=.95,'Date.Sampled',na.rm=T)
plot(TUKEY.nag,las=1,col="brown")
generate_label_df<-function(TUKEY.nag,Date.Sampled){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY.bg[[Date.Sampled]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
LABELS.NAG=generate_label_df(TUKEY.nag,"Date.Sampled")

p.aov=aov(PT2~Date.Sampled,data=eea.df)
TUKEY.p=TukeyHSD(p.aov,conf.level=.95,'Date.Sampled',na.rm=T)
plot(TUKEY.p,las=1,col="brown")
generate_label_df<-function(TUKEY.p,Date.Sampled){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY.bg[[Date.Sampled]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
LABELS.P=generate_label_df(TUKEY.p,"Date.Sampled")

capture.output(LABELS.BG,LABELS.NAG,LABELS.P,file="enz_tukey.docx")


bg.sum=summarySE(eea.df,measurevar="BGT2",groupvars=c("time","Sample"))
nag.sum=summarySE(eea.df,measurevar="NAGT2",groupvars=c("time","Sample"))
p.sum=summarySE(eea.df,measurevar="PT2",groupvars=c("time","Sample"))

exp.plot.bg=ggplot(bg.sum,aes(time,group=Sample,na.rm=T))+
  geom_point(data=bg.sum,aes(y=BGT2, shape=Sample),stat="identity",color="black",size=3)+
  geom_line(aes(y=BGT2, linetype=Sample),width=.5)+
  geom_errorbar(aes(ymin=BGT2-se,ymax=BGT2+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote(atop(beta~'glucosidase', 
                   '( '*mu~'g'~cm^-2~h^-1~')')))+
  theme(legend.position=c(0.05,0.65),legend.title=element_blank())+
  scale_x_continuous(breaks=seq(0,30,5))+
  xlab("Days of Incubation")+annotate("text",x=0,y=200, label="paste(bold(a))",parse=T)


exp.plot.bg


exp.plot.nag=ggplot(nag.sum,aes(time,group=Sample,na.rm=T))+
  geom_point(data=nag.sum,aes(y=NAGT2, shape=Sample),stat="identity",color="black",size=3)+
  geom_line(aes(y=NAGT2, linetype=Sample),width=.5)+
  geom_errorbar(aes(ymin=NAGT2-se,ymax=NAGT2+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote(atop('N-acetyl- '*beta~'-D glucosaminidase', '( '*mu~'g'~cm^-2~h^-1~')')))+
  theme(legend.position=c(0.05,0.85),legend.title=element_blank())+
  scale_x_continuous(breaks=seq(0,30,5))+
  xlab("Days of Incubation")+annotate("text",x=0,y=200, label="paste(bold(b))",parse=T,hjust=0,size=4)
exp.plot.nag


exp.plot.p=ggplot(p.sum,aes(time,group=Sample,na.rm=T))+
  geom_point(data=p.sum,aes(y=PT2, shape=Sample),stat="identity",color="black",size=3)+
  geom_line(aes(y=PT2, linetype=Sample),width=.5)+
  geom_errorbar(aes(ymin=PT2-se,ymax=PT2+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote(atop('Phosphatase', '( '*mu~'g'~cm^-2~h^-1~')')))+xlab("Date")+scale_y_continuous(limits=c(0,200))+
  scale_x_continuous(breaks=seq(0,30,5))+
  theme(legend.position=c(0.05,0.85),legend.title=element_blank())+
  scale_x_continuous(breaks=seq(0,30,5))+
  xlab("Days of Incubation")+annotate("text",x=0,y=200, label="paste(bold(c))",parse=T,hjust=0,size=4)

exp.plot.p


enz.plots=ggarrange(exp.plot.bg+theme(axis.text.x=element_blank(),axis.title.x = element_blank()),
                    exp.plot.nag+theme(legend.position = "none",axis.text.x=element_blank(),axis.title.x = element_blank()),
                    exp.plot.p+theme(legend.position = "none"),nrow=3,align="v",common.legend = FALSE,heights=c(1,1,1.25) )
ggsave("eea_compiled.jpeg",device="jpeg", width=8.5,height=11, units="in")
cols=c(13,15,17)
iv.string=names(combo.df)[cols]


#Enz_OD plots
bg.sum.od=summarySE(eea.df,measurevar="BGT2_OD",groupvars=c("time","Sample"))
nag.sum.od=summarySE(eea.df,measurevar="NAGT2_OD",groupvars=c("time","Sample"))
p.sum.od=summarySE(eea.df,measurevar="PT2_OD",groupvars=c("time","Sample"))

exp.plot.bg.od=ggplot(bg.sum.od,aes(time,group=Sample,na.rm=T))+
  geom_point(data=bg.sum.od,aes(y=BGT2_OD, fill=Sample),stat="identity",color="black",size=3,pch=21)+
  geom_line(aes(y=BGT2_OD, color=Sample),width=.5)+
  #geom_hline(yintercept=0)+theme_classic()+
  geom_errorbar(aes(ymin=BGT2_OD-se,ymax=BGT2_OD+se,fill=Sample,color=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote(atop(beta~'glucosidase', 
                   '( '*mu~'g'~cm^-2~h^-1~')')))+
  scale_x_continuous(breaks=seq(0,30,5))
#geom_boxplot(aes(y=BGT2,group=time),alpha=0.01,position=position_dodge())

exp.plot.bg.od
#ggsave("bgt2.od.jpeg",device="jpeg", width=6.5,height=6.2, units="in")

exp.plot.nag.od=ggplot(nag.sum.od,aes(time,group=Sample,na.rm=T))+
  geom_point(data=nag.sum.od,aes(y=NAGT2_OD, fill=Sample),stat="identity",color="black",size=3,pch=21)+
  geom_line(aes(y=NAGT2_OD, color=Sample),width=.5)+
  #geom_hline(yintercept=0)+theme_classic()+
  geom_errorbar(aes(ymin=NAGT2_OD-se,ymax=NAGT2_OD+se,fill=Sample,color=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote(atop('N-acetyl- '*beta~'-D glucosaminidase', '( '*mu~'g'~cm^-2~h^-1~')')))+
  scale_x_continuous(breaks=seq(0,30,5))
#geom_boxplot(aes(y=NAGT2,group=time),alpha=0.01,position=position_dodge())
exp.plot.nag.od
#ggsave("nagt2.od.jpeg",device="jpeg", width=6.5,height=6.2, units="in")


exp.plot.p.od=ggplot(p.sum.od,aes(time,group=Sample,na.rm=T))+
  geom_point(data=p.sum.od,aes(y=PT2_OD, fill=Sample),stat="identity",color="black",size=3,pch=21)+
  geom_line(aes(y=PT2_OD, color=Sample),width=.5)+
  geom_errorbar(aes(ymin=PT2_OD-se,ymax=PT2_OD+se,fill=Sample,color=Sample),
                stat="identity",width=.15,na.rm=T)+
  ylab(bquote(atop('Phosphatase', '( '*mu~'g'~cm^-2~h^-1~')')))+xlab("Date")+
  scale_x_continuous(breaks=seq(0,30,5))+
  xlab("Days of Incubation")
#geom_boxplot(aes(y=PT2,group=time),alpha=0.01,position=position_dodge())
exp.plot.p.od
#ggsave("pt2.od.jpeg",device="jpeg", width=6.5,height=6.2, units="in")

####Enz-metabolism regression plots####
write.csv(combo.df,file="combodf.csv")
require(ggpmisc)
reg.plots<-function(combo.df,na.rm=TRUE,...){
  sample_list=unique(combo.df$Sample)
  rv.string=names(combo.df)[5:7]
  cols=c(11,13,15)
  iv.string=names(combo.df)[cols]
  plots=list()
  
  for (i in seq_along(sample_list)){
    for(j in rv.string){
      for(k in iv.string){
        
        plots<-ggplot(combo.df[combo.df$Sample==levels(combo.df$Sample)[i],],
                      aes_string(x=k,y=j))+
          geom_point(aes(colour=Date.Sampled))+
          geom_smooth(inherit.aes = TRUE, method="lm")+
          ggtitle(paste(sample_list[i]))+
          stat_poly_eq(formula=y~x,aes(label=paste(..eq.label..,..rr.label..,sep="~~~")),parse=TRUE)+
          stat_fit_glance(method='lm',method.args=list(formula=y~x),geom='text',aes(label=paste("P-value = ",signif(..p.value..,digits=4),sep="")),label.x.npc = 'right', label.y.npc = 0.35, size = 3)
        png(paste(i,j,k,".png",sep=""),width=600, height=500,res=120,units="px")
        print(plots)
        dev.off()
      }
    }
  }
}

#reg.plots(combo.df)



####PCA####
library(tidyr)
library(devtools)
install_github("vqv/ggbiplot")
install_github("ropenscilabs/packagemetrics")
require(packagemetrics)
require(ggbiplot)
library(pander)
devtools::install_github("davidgohel/flextable")
require(flextable)
library(data.table)
panderOptions('digits',4)
panderOptions('round',4)
panderOptions('keep.trailing.zeros',TRUE)
cols=c(13,15,17)
iv.string=names(combo.df)[cols]

#+OD_cm

pca=prcomp(~Resp+GPP+NEP+BGT2+NAGT2+PT2+chl_abs,data=pc.df,scale=TRUE,center=T)
#srp+nh4+no3++xvelocity+width+xdepth+xtemp+xspc+Q
pca1=summary(pca,scale=TRUE)
capture.output(pca1,file="pcasum.doc")

pca$sdev^2
capture.output(pca$sdev^2,file="eigenvalues.doc") #eigenvalues  #retain 1 and 2, 3 by kaiser criterios
predict(pca,newdata=tail(pc.df,2))
pca.group=pc.df[2]
pca.group2=pc.df[11]
print(pca)
plot(pca,type="l")

pca.tab=as.data.frame(pca$rotation[,1:2])
pca.tab=setattr(pca.tab,"row.names", 
                c("Respiration","GPP","NEP","BG","NAG",
                  "Phosphatase","Chlorophyll","Biomass"))




pca.vars=pca.tab$sdev^2

writeWorksheetToFile("pcaoutput.xlsx",
                     data=list(pca.tab),
                     sheet=c("pca"),
                     header=T,clearSheets = T,rownames="var")


p=ggbiplot(pca,obs.scale=1,
           var.scale=1,
           groups=pca.group2$time,ellipse=TRUE,var.axes=F,
           circle=FALSE)
p=p+theme_classic()+geom_point(pch=21,color="black",aes(color=pca.group2$time))+
  scale_color_gradient(low="black")+
  theme(legend.position=c(0.85,0.25))+
  scale_x_continuous(limits=c(-4,7))+
  scale_y_continuous(limits=c(-4,4))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14))

p
#ggsave("pca.date.jpeg",device="jpeg", width=10.5,height=6, units="in")


q=ggbiplot(pca,obs.scale=1,
           var.scale=1,var.axes=F,
           groups=pca.group$Sample,ellipse=TRUE,circle=FALSE)
q=q+theme_classic()+
  geom_point(aes(shape=pca.group$Sample,color=pca.group$Sample,fill=pca.group$Sample),size=4,color="black")+
  scale_shape_manual(name="Substrate",values=c(21,22,23,24,25))+scale_color_discrete(name="Substrate")+
  scale_x_continuous(limits=c(-4,7))+
  scale_y_continuous(limits=c(-4,4))+
  theme(legend.position=c(0.85,0.85))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_blank())
q  
#ggsave("pca.sample.jpeg",device="jpeg", width=8.5,height=6, units="in")

# 
pca.bw=ggarrange(q,p, nrow=2)
#ggsave("pca_plots_bw.jpeg",device="jpeg",width=8.5, height=11,units="in")

# Data for the supplementary quantitative variables
quanti.sup <- pc.df[1:21, 12:21, drop = FALSE]
head(quanti.sup)

# Calculate the correlations between supplementary variables
# and the principal components
ind.coord <- pca$x
quanti.coord <- cor(quanti.sup, ind.coord)
head(quanti.coord[, 1:2])




# Helper function : 
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

# Variable correlation/coordinates
loadings <- pca$rotation
sdev <- pca$sdev

var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord[, 1:4])

#ggarrange(p,q,align="v",nrow=2)
#ggsave("compiledpca.jpeg",device="jpeg", width=6.5,height=8.2, units="in")
#####STacked bar plots####
require(reshape)
require(dplyr)#,position=position_dodge(.9)
require(Rmisc)
require(ggplot2)
require(multcompView)
require(cowplot)
require(gridExtra)
require(ggpubr)





bar.df=melt(metabolism.data,id.vars=c("time","Sample"),measure.vars = c("GPP","Resp","NEP"),na.rm=TRUE)
head(bar.df)
bar.sum=summarySE(bar.df,measurevar="value",groupvars = c("time","variable","Sample"),na.rm=TRUE)
resp.gpp.bar.plot=ggplot(bar.sum,aes(time,group=Sample,shape=Sample))+
  geom_point(data=subset(bar.sum,variable=="GPP"), aes(y=value,shape=Sample),stat="identity",size=3)+
  geom_point(data=subset(bar.sum,variable=="Resp"),aes(y=value,shape=Sample),stat="identity",color="black",size=3)+
  geom_line(data=subset(bar.sum,variable=="GPP"), aes(y=value, linetype=Sample),width=.5)+
  geom_line(data=subset(bar.sum,variable=="Resp"), aes(y=value, linetype=Sample),width=.5)+
  geom_hline(yintercept=0)+theme_classic()+
  geom_errorbar(data=subset(bar.sum,variable=="GPP"),aes(ymin=value-se,ymax=value+se,fill=Sample),stat="identity",width=.15,na.rm=T)+
  geom_errorbar(data=subset(bar.sum,variable=="Resp"),aes(ymin=value-se,ymax=value+se,fill=Sample),stat="identity",width=.15,na.rm=T)+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = c(0.1,0.85),
        legend.title=element_blank())+
  scale_y_continuous(limits=c(-1.15,1.15))+scale_x_continuous(breaks=seq(0,30,5))+
  ylab(bquote('mg O'[2]~''~cm^-2~''~hr^-1))+xlab(NULL)+annotate("text",x=30,y=1, label="GPP",hjust=0)+annotate("text",x=30,y=-1,label="Respiration",hjust=0)


nep.bar.plot=ggplot(bar.sum,aes(time,group=Sample,shape=variable))+
  geom_point(data=subset(bar.sum,variable=="NEP"),
             aes(y=value, shape=Sample),stat="identity",
             color="black",size=3)+
  geom_line(data=subset(bar.sum,variable=="NEP"), 
            aes(y=value, linetype=Sample),width=.5)+
  geom_hline(yintercept=0)+theme_classic()+
  geom_errorbar(data=subset(bar.sum,variable=="NEP"),
                aes(ymin=value-se,ymax=value+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  scale_y_continuous(limits=c(-1.15,1.15))+
  scale_x_continuous(breaks=seq(0,30,5))+
  ylab(bquote('mg O'[2]~''~cm^-2~''~hr^-1))+
  xlab("Days of Incubation")+
  annotate("text",x=30,y=1, label="NEP",hjust=0)+
  theme(legend.position = c(0.15,0.85),legend.title=element_blank())



#stacked.bar=ggarrange(resp.gpp.bar.plot+theme(axis.line.x=element_blank(),axis.text.x=element_blank(),
# axis.ticks.x=element_blank()),nep.bar.plot+theme(legend.position = "none"),nrow=2,common.legend = F)
#stacked.bar
#("met.jpeg",device="jpeg", width=8.5,height=6.2, units="in")
#
#breakup GPP RESP

bar.sum=summarySE(bar.df,measurevar="value",groupvars = c("time","variable","Sample"),na.rm=TRUE)
na.omit(bar.sum)
resp.bar.plot=ggplot(bar.sum,aes(time,group=Sample,shape=Sample))+
  geom_point(data=subset(bar.sum,variable=="Resp"),
             aes(y=value,shape=Sample),stat="identity",
             color="black",
             size=3)+
  geom_line(data=subset(bar.sum,variable=="Resp"), 
            aes(y=value, linetype=Sample),
            width=.5)+
  theme_classic()+
  geom_errorbar(data=subset(bar.sum,variable=="Resp"),
                aes(ymin=value-se,ymax=value+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  scale_x_continuous(breaks=seq(0,30,5),position="top")+
  ylab(bquote(atop('Respiration','mg O'[2]~''~cm^-2~''~hr^-1)))+
  xlab(NULL)+theme(legend.position = c(0.15,0.85),
                   legend.title=element_blank(),
                   axis.text=element_text(size=16),
                   axis.title=element_text(size=16),
                   legend.text=element_text(size=14))

gpp.bar.plot=ggplot(bar.sum,aes(time,group=Sample,shape=Sample))+
  geom_point(data=subset(bar.sum,variable=="GPP"), 
             aes(y=value,shape=Sample),stat="identity",size=3)+
  geom_line(data=subset(bar.sum,variable=="GPP"), 
            aes(y=value, linetype=Sample),width=.5)+
  theme_classic()+
  geom_errorbar(data=subset(bar.sum,variable=="GPP"),
                aes(ymin=value-se,ymax=value+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  scale_x_continuous(breaks=seq(0,30,5))+
  ylab(bquote(atop('GPP', 'mg O'[2]~''~cm^-2~''~hr^-1)))+
  xlab(NULL)+theme(legend.position = c(0.15,0.85),
                   legend.title=element_blank(),
                   axis.text=element_text(size=16),
                   axis.title=element_text(size=16),
                   legend.text=element_text(size=14))

nep.bar.plot=ggplot(bar.sum,aes(time,shape=Sample,group=Sample))+
  geom_point(data=subset(bar.sum,variable=="NEP"),na.rm=T,
             aes(y=value, shape=Sample),stat="identity",color="black",size=3)+
  geom_line(data=subset(bar.sum,variable=="NEP"),na.rm=T, 
            aes(y=value, linetype=Sample),width=.5)+
  geom_hline(yintercept=0)+theme_classic()+
  geom_errorbar(data=subset(bar.sum,variable=="NEP"),
                aes(ymin=value-se,ymax=value+se,fill=Sample),
                stat="identity",width=.15,na.rm=T)+
  scale_x_continuous(breaks=seq(0,30,5))+
  ylab(bquote(atop('NEP','mg O'[2]~''~cm^-2~''~hr^-1)))+
  xlab("Days of Incubation")+
  scale_y_continuous(breaks=seq(-1,1,0.5))+
  theme(legend.position = c(0.15,0.85),
        legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=14))

stacked.bar=ggarrange(gpp.bar.plot+theme(axis.text.x=element_blank()),
                      resp.bar.plot+theme(axis.text.x=element_blank()),
                      nep.bar.plot,nrow=3,common.legend = TRUE,legend="bottom",align="v",heights=c(1,1,1.25))
ggsave("stacked_met.jpeg",device="jpeg",width=8.5, height=8.5,units="in")


####PAR and Chl####
require(ggplot2)
par.lm=lm(int_par$Summed.PAR~int_par$Date.Sampled)
psummary(par.lm)
partable=
  int_par$Date.Sampled=as.Date(int_par$Date.Sampled,format='%m/%d/%y')
par.plot=ggplot(data=int_par,aes(x=Date.Sampled,y=Summed.PAR))
par.plot+geom_smooth(method="auto")+geom_point()+scale_color_gradient(low="blue",high="red")#+geom_text(x=as.Date("2017-08-29"),y=5000,label=lm_eqn(int_par),parse=T)


####field data####
library(lubridate)

hobo.a=read.csv("NBC-A.csv",header=T)
hobo.b=read.csv("NBC-B.csv",header=T)
hobo=cbind(hobo.a,hobo.b$Temp,hobo.b$Intensity)
hobo$datetime=as.POSIXct(paste(hobo$Date,hobo$Time), format="%m/%d/%Y %H:%M:%S")
usgs.df$datetime=mdy_hm(usgs.df$timestamp)
field.data$disdate=as.POSIXct(paste(field.data$Date),format="%m/%d/%Y")
#field.data$disdate=mdy(field.data$Date)#
par.df$date.time=as.POSIXct(paste(par.df$Date,par.df$Time),format="%m/%d/%y %H:%M:%S")
par.df$datetime=ymd_hms(par.df$date.time)
water$Date=mdy(water$Date) #formatting issue

water$no310=water$no3/10
water.long=melt(water,id.vars="Date",measure.vars = c("srp","nh4","no310"),variable_name="compound")
water.sum=summarySE(water.long,measurevar="value",groupvars = c("Date","compound"),na.rm=TRUE)

water.sum$Ions=mapvalues(water.sum$compound,from=c("srp","nh4","no310"),to=c("SRP","NH4","NO3"))

####water chem by date####
srp.lm=lm(srp~Date,data=water)
summary(srp.lm)

no3.lm=lm(no3~Date,data=water)
summary(no3.lm)

nh4.lm=lm(nh4~Date,data=water)
summary(nh4.lm)


###water chem plots####
q.density.df=read.csv("usgs_seasonal.csv")

dil.fac=10
chem.plot=ggplot(data=water,aes(x=Date,y=srp))+geom_point()+geom_smooth(data=water,aes(x=Date,y=srp),color="black")+
  geom_point(data=water,aes(x=Date,y=no3/dil.fac),color="Red")+geom_smooth(data=water,aes(x=Date,y=no3/dil.fac),color="red")+
  geom_point(data=water,aes(x=Date,y=nh4),color="Blue")+geom_smooth(data=water,aes(x=Date,y=nh4),color="blue")
chem.plot+scale_y_continuous(sec.axis=sec_axis(~.*dil.fac,name="NO3"))+ylab(bquote("SRP, NH4"))

chem.plot=ggplot(data=water.sum,aes(x=Date,y=value,group=Ions,color=Ions))+geom_point()+
  geom_errorbar(data=water.sum, aes(ymin=value-se,ymax=value+se,fill=Ions,color=Ions),
                stat="identity",width=.15,na.rm=T)+
  geom_line()+scale_y_continuous(sec.axis=sec_axis(~.*dil.fac,name=bquote(atop('NO'[3]~'','('*mu~'g  l'^-1~')'))))+
  ylab(bquote(atop('SRP, NH'[4]~'','('*mu~'g  l'^-1~')')))+geom_smooth(method="lm",se=F,color="black")+
  scale_x_date(breaks=pretty_breaks(6))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chem.plot
ggarrange(chem.plot,legend="bottom")
ggsave("waterchem.jpeg",device="jpeg",width=6.5,height=6.5,units="in")


field.plot=ggplot(data=usgs.df,aes(x=datetime,y=Q*.0283168))
field.plot=field.plot+geom_point()+geom_point(data=field.data, aes(x=disdate,y=Q),color="red")+geom_line()+
  ylab(bquote('(Q) '~m^-3~''~s^-1))+xlab(bquote("Date"))+scale_x_datetime(breaks=pretty_breaks(9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_y_continuous(breaks=seq(0,4,0.5))
ggsave("discharge.jpeg",device="jpeg",width=6.5,height=6.5,units="in")


td.plot=ggplot(hobo,aes(x=datetime,y=Temp))+geom_point()
td.plot=td.plot+theme_classic()+geom_line()+
  geom_point(data=par.df,aes(x=date.time,y=microMole.s.m2/4000),color="red")+
  geom_line(data=par.df, aes(x=date.time,y=microMole.s.m2/4000),color="red")+
  scale_y_continuous(sec.axis=sec_axis(~.*4000,name=bquote(atop('PAR','(mol photons m'^-2~' s'^-1~')'))))+
  xlab(bquote("Date"))+
  ylab(bquote(atop('Water Temperature','(C'*degree~')')))
ggsave("temp_depth.jpeg",device="jpeg",width=6.5,height=6.5,units="in")

ggarrange(td.plot,field.plot,nrow=2)
ggsave("Qtemppar.jpeg",device="jpeg",width=6.5,height=6.5, units="in")

####water chem and enzymes####
p.enz.plot=ggplot(eea.df,aes(x=srp,y=PT2,color=Sample))+geom_point()
p.enz.plot+geom_smooth(method="lm",se=F)

nh4.enz.plot=ggplot(eea.df,aes(x=nh4,y=NAGT2,color=Sample))+geom_point()
nh4.enz.plot+geom_smooth(method="lm",se=F)

no3.enz.plot=ggplot(eea.df,aes(x=no3,y=NAGT2,color=Sample))+geom_point()
no3.enz.plot+geom_smooth(method="lm",se=F)

eea.df=na.omit(eea.df)
p.lme=lmList(PT2~srp|Sample,eea.df)
summary(p.lme)

no3.lme=lmList(NAGT2~no3|Sample,eea.df)
summary(no3.lme)

nh4.lme=lmList(NAGT2~nh4|Sample,eea.df)
summary(nh4.lme)

n.lme=lmList(NAGT2~no3+nh4|Sample,eea.df)
summary(n.lme)


####ciation####
ip = sessionInfo()
ip=ip$loadedOnly

