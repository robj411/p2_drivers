library(data.table)
library(dplyr)
library(lubridate)
library(OECD)
library(readODS)
library(readxl)
library(xlsx)
library(stringi)
library(tidyverse)
library(splines)
library(zoo)
library(splines)
library(BayesTools)
library(runjags)
library(rjags)
library(extraDistr)

## paths and functions ####################################

data_path <- '../../data'
figure_path <- '../../README_files/figure-gfm/'

epiwrap <- function(x){
  out <- c()
  for(i in 1:length(x)) out[i] <- ceiling(difftime(x[i],'2019-12-28',units='weeks')[[1]])
  out
}

get_month_day <- function(x){
  out <- c()
  for(i in 1:length(x)) out[i] <- ceiling(difftime(paste0(x[i],'-01'),'2019-12-31',units='days')[[1]])
  out
}
get_day_from_jan1 <- function(x){
  out <- c()
  for(i in 1:length(x)) out[i] <- ceiling(difftime(x[i],'2019-12-31',units='days')[[1]])
  out
} 

get_day_from_week <- function(x){
  as.Date('2019-12-31') + weeks(x)
}

readxlsx <- function(...) suppressMessages(read_xlsx(...))

smooth_counts <- function(x,window=3){
  y <- x
  for(i in 1:length(x)){
    ind <- max(1,i-window) : min(length(x),i+window)
    y[i] <- mean(x[ind],na.rm=T)
  }
  y
}



## data #######################################

owid <- read.csv('https://covid.ourworldindata.org/data/owid-covid-data.csv',stringsAsFactors = F)

mobilityfull <- setDT(read.csv(file.path(data_path,'12.mobility.csv')))
response <- setDT(read.csv(file.path(data_path,'11.response.csv')))
deathcsv <- read.csv(file.path(data_path,'12.deaths.csv'))
# excessdeathcsv <- read.csv(file.path(data_path,'12.excess_deaths.csv'),allowEscapes=F)

## process mobility ##############################

mobilityfull[,date:=as.Date(Day,"%d/%m/%Y")]
mobilityfull$Day <- NULL
mobilityfull[,Entity:=gsub(' ','.',Entity)]
mobilityfull[,Entity:=gsub("'",'.',Entity)]
mobilityfull[,Entity:=gsub('_','',Entity)]
mobilityfull$parks <- NULL
mobilityfull$residential <- NULL
mobilityfull[,mob:=mean(c(retail_and_recreation,transit_stations,workplaces),na.rm=T),by=.(Code,date)]
mobilityfull$grocery_and_pharmacy <- NULL
mobilityfull$retail_and_recreation <- NULL
mobilityfull$transit_stations <- NULL
mobilityfull$workplaces <- NULL

setorder(mobilityfull,'Code','date')

## smooth mobility ##########################

kn <- 80
mobilityfull[,smoothmob:=lapply(.SD,function(z){
  fn<-smooth.spline(x=as.numeric(date),y=z,nknots=kn);
  predict(fn,as.numeric(date))$y
}),by=Code,.SDcols='mob']
mobilityfull[,deriv:=lapply(.SD,function(z){
  fn<-smooth.spline(x=as.numeric(date),y=z,nknots=kn);
  andt <- as.numeric(date)
  predt <- sort(c(andt,andt+diff(andt)[1]/10))
  diff(predict(fn,predt)$y)[seq(1,length(predt),by=2)]
}),by=Code,.SDcols='smoothmob']
mobilityfull[,minimum:=lapply(.SD,function(z){
  c(z[1:(length(z)-1)] <= 0 & z[2:(length(z))] >= 0,F)
}),by=Code,.SDcols='deriv']
mobilityfull[,maximum:=lapply(.SD,function(z){
  c(z[1:(length(z)-1)] >= 0 & z[2:(length(z))] <= 0,F)
}),by=Code,.SDcols='deriv']
mobilityfull[,maximum:=lapply(.SD,function(z){
  maxes <- c(z[1:(length(z)-1)] >= 0 & z[2:(length(z))] <= 0,F)
  mins <- c(z[1:(length(z)-1)] <= 0 & z[2:(length(z))] >= 0,F)
  if(which(maxes)[1] > which(mins)[1]) maxes[1] <- T
  maxes
}),by=Code,.SDcols='deriv']

## extract mins and maxes ########################

minmax <- subset(mobilityfull,(maximum|minimum)&date<'2020-12-31')
minmax[,prev:=c(0,smoothmob[-length(smoothmob)]),by=Code]
minmax[,loss:=100*(smoothmob-prev)/(100+prev)]
minmax[,minval:=min(loss),by=Code]
minmob <- minmax[minval==loss,.(Code,Entity,date,minval,smoothmob)]
ggplot(mobilityfull) + 
  geom_line(aes(x=date,y=smoothmob,colour=Code),show.legend = F) +
  geom_point(data=subset(mobilityfull,minimum==T),aes(x=date,y=smoothmob,colour=Code),show.legend = F) 

## find largest drop

mobility <- mobilityfull[date<as.Date("2020-10-01"),]
mobility[,minval:=min(smoothmob),by=Code]
mobility[,minmin:=minval==smoothmob,by=Code]

ggplot(subset(mobilityfull,date<'2020-10-01')) + 
  geom_line(aes(x=date,y=smoothmob,colour=Code),show.legend = F) +
  geom_point(data=minmob,aes(x=date,y=smoothmob,colour=Code),show.legend = F) 
# geom_point(data=subset(mobility,minmin==T),aes(x=date,y=smoothmob,colour=Code),show.legend = F) 

ggplot(subset(mobility,Code%in%subset(mobility,minmin&date>"2020-06-01")$Code)) + 
  geom_line(aes(x=date,y=smoothmob,colour=Code)) +
  geom_point(data=subset(minmob,date>"2020-06-01"),aes(x=date,y=smoothmob,colour=Code),show.legend = F) 
# geom_point(data=subset(mobility,minmin==T&Code%in%subset(mobility,minmin&date>"2020-06-01")$Code),aes(x=date,y=smoothmob,colour=Code),show.legend = F) 
ggplot(subset(mobility,Code%in%subset(mobility,minmin&date>"2020-06-01")$Code)) + 
  geom_line(aes(x=date,y=mob,colour=Code))

# minmob <- mobility[minmin==T,.(Code,Entity,date,minval)]

## process deaths data ################################

# excessdeathcsv$location_name <- gsub(' ','.',excessdeathcsv$location_name)
# excessdeathcsv$location_name <- gsub('\'','.',excessdeathcsv$location_name)
# excessdeathcsv <- subset(excessdeathcsv,location_name%in%minmob$Entity)
deathcsv <- deathcsv[,c(1,which(colnames(deathcsv)%in%minmob$Entity))]
deathcsv[is.na(deathcsv)] <- 0  

## subset for which countries we have data ##################

# minmob <- subset(minmob,Entity%in%colnames(deathcsv)&Entity%in%excessdeathcsv$location_name&Entity%in%colnames(response))
minmob <- subset(minmob,Entity%in%colnames(deathcsv)&Entity%in%colnames(response))

# combine datasets
minmob[,deaths:=mean(deathcsv[[Entity]][as.Date(deathcsv$date,"%d/%m/%y")%in%(date+c(-3:3))])
       # *excessdeathcsv$mean_value[excessdeathcsv$location_name==Entity]
       ,by=.(Entity,date)]
minmob[,response:=response[[Entity]][as.Date(response$country_name,"%d/%m/%Y")==date],by=.(Entity,date)]

## process owid data #################
# maybe only use for continent....

owid2020 <- setDT(subset(owid,date%in%minmob$date&iso_code%in%minmob$Code))
councon <- unique(owid2020[,.(continent,iso_code)])

minmob[,continent:=councon$continent[councon$iso_code==Code],by=.(Entity,date)]


# with(minmob, table(incomelevels$IncomeGroup[match(Code,incomelevels$Country.Code)]))

## plot #########################################

ggplot(subset(mobilityfull,Code%in%subset(minmob,response<40)$Code)) + 
  theme_bw(base_size = 15) +
  labs(x='',y='Smoothed average mobility',colour='') +
  geom_line(aes(x=date,y=smoothmob,colour=Entity)) +
  geom_point(data=subset(minmob,response<40),aes(x=date,y=minval,colour=Entity),show.legend = F)  +
  theme(legend.position = 'top')
p <- ggplot(subset(mobilityfull,date<'2020-12-31')) + 
  theme_bw(base_size = 15) +
  labs(x='',y='Smoothed average mobility',colour='') +
  geom_line(aes(x=date,y=smoothmob,colour=Entity),show.legend = F,alpha=.5) +
  geom_point(data=subset(minmob,date<'2020-12-01'),aes(x=date,y=smoothmob,colour=Entity),show.legend = F) 
ggsave(p,filename=file.path(figure_path,'smoothmobility.png'))

p <- ggplot(minmob,aes(x=response,y=minval)) +
  geom_text(aes(colour=continent,label=Code),size=4,show.legend = F) +
  # geom_smooth(data=subset(minmob,response>50),formula=y~ns(x,df=1),method=glm,se = F) +
  theme_bw(base_size = 15) +
  labs(x='Stringency at largest drop in mobility',y='Largest drop in mobility',colour='') +
  scale_y_continuous(limits = c(0,100)-100) +
  # geom_text(data=subset(minmob,Code%in%c('NER','BFA','YEM')),aes(label=paste0(Entity,' (',substr(date,1,4),')')),vjust=2) +
  theme(legend.position = 'top') 

ggsave(p,filename=file.path(figure_path,'mobilitydrop.png'),width=10)

## bayestools #################################

# format data
data <- list(
  mobility = (100+minmob$smoothmob)/100,
  deathspp = minmob$deaths,
  response = minmob$response/10,
  N = nrow(minmob)
)

## create and fit models
# define priors
p1 <- BayesTools::prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
p3 <- BayesTools::prior(distribution = "normal", parameters = list(mean = 0, sd = 1), truncation = list(0, Inf))
p4 <- BayesTools::prior(distribution = "normal", parameters = list(mean = 10, sd = 5), truncation = list(0, Inf))

priors_list1 <- list(mu1 = p3, 
                     sd1 = p3, 
                     mu2 = p3, 
                     sd2 = p3, 
                     mu3 = p1, 
                     sd3 = p3, 
                     mub = p1, 
                     sdb = p3, 
                     scale = p4)

# sketch the shape of the function
x <- 0:30
b <- 1/3
a <- -1
r <- .75
plot(x,1/((1+x+r^exp(a)))*(1-b)+b)

# define likelihood for the data
model_syntax <-
  "model{
    for(i in 1:N){
      # rate per country
      k1i[i] ~ dnorm(mu1, sd1)T(0,)
      k2i[i] ~ dnorm(mu2, sd2)T(0,)
      # lower bound per country
      bi[i] ~ dnorm(mub, sdb)
      bl[i] = exp(bi[i])/(1+exp(bi[i]))
      # response is on a scale 0 to 10
      r[i] ~ dbeta(1+response[i], 1+10-response[i])
      # expected mobility
      # divide deaths by 10 to get scale
      # multiply r by 0.9 because the mandate cannot be '100' as it is modelled in daedalus
      v1[i] = (k1i[i])*deathspp[i]
      v2[i] = (k2i[i])*(r[i]*0.9)
      mout[i] = 1/((1+v1[i]+v2[i]))*(1-bl[i])+bl[i]
      mobility[i] ~ dbeta(.1 + 5*mout[i], .1 + 5*(1-mout[i]))
    }
  }"

# fit the models
fit1 <- JAGS_fit(model_syntax, data, priors_list1, seed = 1,thin = 10, 
                 add_parameters = c('k1i','k2i','bl','mout'))

allsamples <- as.data.frame(do.call(rbind,fit1$mcmc))
plot(allsamples$mub,type='l')
plot((cbind(allsamples$mu1,allsamples$mub)))
plot((cbind(allsamples$mu1,allsamples$mu2)))

(smm <- summary(fit1))
indices <- grepl('mout',rownames(smm))
minmob$predmedian <- smm[indices,2]
minmob$predl95 <- smm[indices,1]
minmob$predu95 <- smm[indices,3]
# minmob$predmean <- smm[indices,4]
plot(data$mobility,minmob$predmedian)

(ggplot(allsamples) + geom_point(aes(mu1,(mub)),alpha=.2) +
  theme_bw(base_size=15) +
  labs(x='Deaths coefficient',y='Logit baseline') -> pmobpost)
# as rate gets larger, population gets more responsive


(ggplot(minmob) + geom_errorbar(aes(x=smoothmob/100+1,ymin=predl95,ymax=predu95),colour='grey') +
  geom_point(aes(x=smoothmob/100+1,y=predmedian))+#,colour=tdeathpp<1)) +
  theme_bw(base_size=15) +
  labs(x='Observed',y='Predicted') -> pmobfit)


# function definition
fx <- function(k1,k2,bl,x, r) {
  v1 <- (k1)*x
  v2 <- (k2)*(r*.9)
  1/(1+(v1+v2))*(1-bl)+bl
}
countrysamples <- setDT(data.frame( k1 = unlist(allsamples[,grepl('k1i',colnames(allsamples))]),
                                    k2 = unlist(allsamples[,grepl('k2i',colnames(allsamples))]),
                                    bl = unlist(allsamples[,grepl('bl',colnames(allsamples))])
))

# extract samples
countrysamples <- setDT(allsamples[,colnames(allsamples)%in%c('mu1','mu2','sd1','sd2','mub','sdb')])
countrysamples[,k1:=rtnorm(nrow(countrysamples),mu1,sd1,a=0)]
countrysamples[,k2:=rtnorm(nrow(countrysamples),mu2,sd2,a=0)]
countrysamples[,bi:=rnorm(nrow(countrysamples),mub,sdb)]
countrysamples[,bl:=exp(bi)/(1+exp(bi))]

# generate some curves for illustration
plot(countrysamples[sample(1:nrow(countrysamples),1000),])
nplot <- 10
sx <- sample(1:nrow(countrysamples),nplot)
k1 <- countrysamples$k1[sx]
k2 <- countrysamples$k2[sx]
b <- countrysamples$bl[sx]
plotsamples <- setDT(expand.grid(x=seq(0,30,by=.1),r=seq(0,1,.25),k1=k1))
nx <- length(unique(plotsamples$x))
plotsamples$b <- rep(b,each=nrow(plotsamples)/nplot)
plotsamples$k2 <- rep(k2,each=nrow(plotsamples)/nplot)
plotsamples$group <- rep(1:(nrow(plotsamples)/nx),each=nx)
plotsamples[,y:=fx(k1,k2,b,x,r),by=.(k1,k2,b,x,r)]
minmob$r <- cut(minmob$response,breaks=c(0,25,50,75,101),labels=c(25,50,75,100))

(ggplot(plotsamples) + geom_line(aes(x,y,colour=as.factor(r*100),group=group)) +
  theme_bw(base_size=15) +
  labs(x='Deaths per million',y='Mobility',colour='Mandate') +
  geom_point(data=minmob,aes(x=deaths,y=smoothmob/100+1,colour=as.factor(r))) -> pcurve)

if(file.exists('../../README_files/figure-gfm')){
  ggsave(pmobpost,filename='../../README_files/figure-gfm/mobilityposterior.png')
  ggsave(pmobfit,filename='../../README_files/figure-gfm/mobilityfitted.png')
  ggsave(pcurve,filename='../../README_files/figure-gfm/mobilitycurves.png')
}

# label and save
newnames <- c('deathcoef','mandatecoef','baseline')
colnames(countrysamples)[match(c('k1','k2','bl'),colnames(countrysamples))] <- newnames
write.csv(countrysamples[sample(1:nrow(countrysamples),2^13),..newnames],
          file.path(data_path,'utr_coefs.csv'),row.names = F)

