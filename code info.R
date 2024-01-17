
####thermal traits estimated by the extended BA model
#for verify rate:u>0

library(minpack.lm)
dat<-read.csv(file.choose())
str(dat)

Johnson<-function(Temp,mu0,Ea,dED,Topt,Tref,kb=8.62E-5){
  stopifnot(dED!=0)
  T0<-273.15
  Tref<-15
  ED<-Ea+dED
  a<-Ea/dED
  b<-exp(ED/kb*(1/(Topt+T0)-1/(Temp+T0)))
  b<-a*b
  a<-Ea/kb*(1/(T0+Tref)-1/(T0+Temp))
  h<-mu0*exp(a)/(1+b)
  return(h)}

u1<-dat[dat$Nutrient=="n1",]
u11<-nlsLM(formula=GR~Johnson(Temp,mu0,Ea,dED,Topt,Tref=15),data=u1,control=nls.lm.control(ftol = 1E-12,ptol = 1E-12,maxiter = 1024),start=list(mu0=1,Ea=0.6,dED=0.5,Topt=20),jac=NULL,algorithm = "LM",lower = c(mu0=1E-6,Ea=0.001,dED=1E-2,Topt=0),upper = c(mu0=10,Ea=10,dED=20,Topt=40))
summary(u11)

#####using GAM to draw the relationship between DIN concentration and thermal traits
### use Topt~s(DIN) as an example
library(ggplot2)
library(tidyverse)
library(mgcv)

dat<-read.csv(file.choose())
Topt1<-gam(Topt ~ s(Nutrient,bs="cs",k=4), data= dat %>% filter(Group == "diatom"))
summary(Topt1)
gam.check(Topt1)

dat %>% ggplot(aes(x = Nutrient, y = Topt, group = Group))+
  geom_point(aes(color = Group,shape = Group),alpha=1,size=3.7)+
     stat_smooth(method="gam",formula = y~s(x, bs="cs",k=4),
              size=1.6,
              aes(color = Group,fill=Group))+
     xlim(0,33)+ylim(26,31.1)+
     scale_colour_manual(values=c("black","blue"))+scale_fill_manual(values = c("grey","lightblue"))

###### using GAMMs to analyze historical data
###use cell density of total diatoms as the example
library(mgcv)
library(itsadug)

TD<-read.csv(file.choose())

gamm_TD<-gamm(logDiatom~s(Temperature,bs="cr")+s(logDIN,bs="cr")+s(logDIP,bs="cr"), data=TD,method='REML',random=list(id=~1))
summary(gamm_TD$lme)
summary(gamm_TD$gam)

gamm_TD_inter1<-gamm(logDiatom~te(Temperature,logDIN)+s(logDIP,bs="cr"), data=TD, method='REML',random=list(id=~1))
summary(gamm_TD_inter1$lme)
summary(gamm_TD_inter1$gam)

gamm_TD_inter1<-gamm(logDiatom~ti(Temperature,logDIN)+s(logDIN,bs="cr")+s(Temperature,bs="cr")+s(logDIP,bs="cr"), data=TD, method='REML',random=list(id=~1))
summary(gamm_TD_inter1$lme)
summary(gamm_TD_inter1$gam)

####Pairplot 

env<-read.csv(file.choose())

panel.cor<-function(x,y,digits=2, prefix="",cex.cor){
  usr<-par("usr")
  on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r<-round(cor(x,y),digits=2)
  txt<-format(c(r,0.123456789), digits=digits)[1]
  txt<-paste(prefix,txt,sep="")
  if(missing(cex.cor)) cex.cor<-0.8/strwidth(txt)
  text(0.5,0.5,txt,cex=cex.cor/1.8)
}


panel.hist <- function(x, ...){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}


pairs(env, upper.panel = panel.cor,
      diag.panel = panel.hist,
      lower.panel = panel.smooth)

#####VIF estimation
library(car)
library(olsrr)
env_diatom<-read.csv(file.choose())
model<-lm(Diatom~Temperature+DIN+DIP+Salinity+Silica, data=env_diatom)
vif(model)
ols_vif_tol(model)
ols_eigen_cindex(model)






