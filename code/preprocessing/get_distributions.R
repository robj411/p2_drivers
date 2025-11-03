## start #######################
library(MASS)
library(data.table)
library(ggplot2)
library(wbstats)
library(dplyr)
library(haven)
library(squire)
library(cowplot)
library(readODS)

datapath <- '../../data'
# yougovdata <- file.path(datapath,'covid-19-tracker/data/')
countrydatafile <- file.path(datapath,'country_data.csv')
incomelevels <- read.csv(file.path(datapath,'Metadata_Country_API_IT.NET.USER.ZS_DS2_en_csv_v2_5455054.csv'))
colnames(incomelevels)[1] <- 'Country.Code'

p2data <- read.csv(countrydatafile)
p2countries <- p2data$country


## internet coverage ##############
print('internet coverage')

internet <- read.csv(file.path(datapath,'API_IT.NET.USER.ZS_DS2_en_csv_v2_5455054.csv'),header = F) ## more recent
internet <- internet[-c(1:3),-c(3:63,65:nrow(internet))]
colnames(internet) <- c('country','Country.Code','internet')

incomeint <- left_join(internet,incomelevels,by='Country.Code')

hic <- fitdistr(subset(incomeint,!is.na(internet)&IncomeGroup=='High income')$internet/100,"beta",start=list(shape1=1,shape2=1))
umic <- fitdistr(subset(incomeint,!is.na(internet)&IncomeGroup%in%c('Upper middle income'))$internet/100,"beta",start=list(shape1=1,shape2=1))
llmic <- fitdistr(subset(incomeint,!is.na(internet)&IncomeGroup%in%c('Low income','Lower middle income'))$internet/100,"beta",start=list(shape1=1,shape2=1))

internet_distributions <- data.frame(parameter_name='internet_coverage',
                                     igroup=c('LLMIC','UMIC','HIC'),
                                     distribution='betainv',
                                     `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape1']]),
                                     `Parameter 2`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape2']]))

## labour share ####################
print('labour share')

labs <- setDT(read_dta(file.path(datapath,'lab_share_data.dta')))
labs[!is.na(labsh),maxyear:=max(year),by=countrycode]
labincome <- left_join(subset(labs,year==maxyear),incomelevels,by=c('countrycode'="Country.Code"))
labincome <- labincome[,.(labsh,IncomeGroup)]

hic <- fitdistr(subset(labincome,!is.na(labsh)&IncomeGroup=='High income')$labsh,"beta",start=list(shape1=1,shape2=1))
umic <- fitdistr(subset(labincome,!is.na(labsh)&IncomeGroup%in%c('Upper middle income'))$labsh,"beta",start=list(shape1=1,shape2=1))
llmic <- fitdistr(subset(labincome,!is.na(labsh)&IncomeGroup%in%c('Low income','Lower middle income'))$labsh,"beta",start=list(shape1=1,shape2=1))

labsh_distributions <- data.frame(parameter_name='labsh',
                                  igroup=c('LLMIC','UMIC','HIC'),
                                  distribution='betainv',
                                  `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape1']]),
                                  `Parameter 2`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape2']]))

## bmi ##########################
print('bmi')

bmi <- read.csv(file.path(datapath,'bmidata.csv'),stringsAs=F)
bmi <- subset(bmi,Dim1=='Both sexes'&!is.na(FactValueNumeric))
setDT(bmi)
bmi[,latest:=max(Period),by=Location]
bmi <- bmi[Period==latest,colnames(bmi)%in%c('SpatialDimValueCode','Location','FactValueNumeric'),with=F]

bmi <- left_join(bmi,incomelevels,by=c('SpatialDimValueCode'="Country.Code"))

hic <- fitdistr(subset(bmi,!is.na(FactValueNumeric)&IncomeGroup=='High income')$FactValueNumeric,"normal")
umic <- fitdistr(subset(bmi,!is.na(FactValueNumeric)&IncomeGroup%in%c('Upper middle income'))$FactValueNumeric,"normal")
llmic <- fitdistr(subset(bmi,!is.na(FactValueNumeric)&IncomeGroup%in%c('Lower middle income','Low income'))$FactValueNumeric,"normal")

bmi_distributions <- data.frame(parameter_name='bmi',
                                igroup=c('LLMIC','UMIC','HIC'),
                                distribution='norminv',
                                `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['mean']]),
                                `Parameter 2`=sapply(list(llmic,umic,hic),function(x)x$estimate[['sd']]))

## hosp capacity #####################
print('hosp capacity')

hic <- fitdistr(subset(p2data,!is.na(Hmax)&igroup=='HIC')$Hmax,"gamma")
umic <- fitdistr(subset(p2data,!is.na(Hmax)&igroup%in%c('UMIC'))$Hmax,"gamma")
llmic <- fitdistr(subset(p2data,!is.na(Hmax)&igroup%in%c('LMIC','LIC'))$Hmax,"gamma")

Hmax_distributions <- data.frame(parameter_name='Hmax',
                                 igroup=c('LLMIC','UMIC','HIC'),
                                 distribution='gaminv',
                                 `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape']]),
                                 `Parameter 2`=sapply(list(llmic,umic,hic),function(x)1/x$estimate[['rate']]))


## workforce in place ##########################
print('workforce in place')

workerNs <- grepl('NNs',colnames(p2data))
workerNpop <- colnames(p2data)%in%paste0('Npop',5:13)
fits <- list()
for(igp in 1:3){
  igroupdata <- subset(p2data,igroup%in%list(c('LMIC','LIC'),'UMIC','HIC')[[igp]])
  workersinplace <- apply(igroupdata[,workerNs],1,sum)
  workeragetotal <- apply(igroupdata[,workerNpop],1,sum)
  fits[[igp]] <- fitdistr(na.omit(workersinplace/workeragetotal),"beta",start=list(shape1=1,shape2=1))
}
workforce_in_place <- data.frame(parameter_name='workforce_in_place',
                                 igroup=c('LLMIC','UMIC','HIC'),
                                 distribution='betainv',
                                 `Parameter 1`=sapply(fits,function(x)x$estimate[['shape1']]),
                                 `Parameter 2`=sapply(fits,function(x)x$estimate[['shape2']]))

## public transport ###################
print('public transport')

# modeshare <- setDT(readxl::read_xlsx('~/projects/DAEDALUS/data/mode_shares.xlsx',sheet=2,skip = 1))
# modeshare2 <- setDT(readxl::read_xlsx('~/projects/DAEDALUS/data/mode_shares.xlsx',sheet=3,skip = 1))
# 
# modeshare[,pt:=(`Bus/Coach`+Rail)/(`Bus/Coach`+Rail+`Passenger Car`)]
# modeshare2[,pt:=(`Public transport (coach, bus, rail)`)/(`Public transport (coach, bus, rail)`+`Passenger Car`)]
# 
# ms <- rbind(modeshare[,c(1,5)],modeshare2[,c(1,4)])
# colnames(ms) <- c('Country','pt')
# 
# 
# ms$Country[ms$Country=="Lao People's Democratic Republic"] <- 'Laos'
# ms$Country[ms$Country=="Russian Federation"] <- 'Russia'
# ms$Country[ms$Country=="Venezuela (Bolivarian Republic of)"] <- 'Venezuela'
# ms$Country[ms$Country=="Bolivia (Plurinational State of)"] <- 'Bolivia'
# ms$Country[ms$Country=="Slovak Republic"] <- 'Slovakia'
# ms$Country[ms$Country=="United Republic of Tanzania"] <- 'Tanzania'
# ms$Country[ms$Country=="Côte d'Ivoire"] <- "Cote d'Ivoire"
# ms$Country[ms$Country=="Bosnia-Herzegovina"] <- "Bosnia and Herzegovina"
# ms$Country[ms$Country=="Czech Republic"] <- "Czechia"
# ms <- subset(ms,Country%in%p2countries)
# 
# p2data$public_transport_share <- NA
# p2data$public_transport_share[match(ms$Country,p2data$country)] <- ms$pt
# 
# hic <- fitdistr(subset(p2data,!is.na(public_transport_share)&igroup=='HIC')$public_transport_share,"beta",start=list(shape1=1,shape2=1))
# umic <- fitdistr(subset(p2data,!is.na(public_transport_share)&igroup=='UMIC')$public_transport_share,"beta",start=list(shape1=1,shape2=1))
# llmic <- fitdistr(subset(p2data,!is.na(public_transport_share)&igroup%in%c('LMIC','LIC'))$public_transport_share,"beta",start=list(shape1=1,shape2=1))
# 
# pt_distributions <- data.frame(parameter_name='pt',
#                                      igroup=c('LLMIC','UMIC','HIC'),
#                                      distribution='betainv',
#                                      `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape1']]),
#                                      `Parameter 2`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape2']]))

## contact matrices ######################
print('contact matrices')

cmcols <- which(grepl('CM',colnames(p2data)))
populations <- squire::population
countries <- unique(populations$country)

p2countries[p2countries=='Palestine'] <- 'State of Palestine'
p2countries[p2countries=='Laos'] <- 'Lao PDR'
p2countries[p2countries=='Taiwan'] <- 'China, Taiwan Province of China'
p2countries[p2countries=='Congo'] <- 'Republic of the Congo'
p2countries[p2countries=='Brunei'] <- 'Brunei Darussalam'
p2countries[p2countries=='Kyrgyzstan'] <- 'Kyrgyz Republic'
p2countries[p2countries=='Hong Kong'] <- 'Hong Kong SAR, China'
p2countries[p2countries=='Macao'] <- 'Macao SAR, China'
p2countries[p2countries=='Virgin Islands US'] <- 'United States Virgin Islands'
p2countries <- gsub('St ','St. ',p2countries)

p2countries[!p2countries%in%countries]
countries[!countries%in%p2countries]

rowindices <- p2countries%in%countries
countries_to_get <- p2countries[rowindices]

# squire:::process_contact_matrix(squire::get_mixing_matrix(cn), squire::get_population(cn,simple_SEIR = T)$n)
matrices <- lapply(countries_to_get,function(cn)
  c(squire:::process_contact_matrix(squire::get_mixing_matrix(cn), squire::get_population(cn,simple_SEIR = T)$n)))
names(matrices) <- countries_to_get



p2data[!rowindices,cmcols] <- NA
for(i in 1:length(matrices)){
  cn <- names(matrices)[i]
  rowindex <- which(p2countries==cn)
  p2data[rowindex,cmcols] <- matrices[[i]]
}
# View(p2data[,c(1,cmcols)])
# 
# write.csv(p2data,countrydatafile,row.names = F, na="")

sqmatrices <- lapply(countries_to_get,function(cn)
  squire:::process_contact_matrix(squire::get_mixing_matrix(cn), squire::get_population(cn,simple_SEIR = T)$n))

popcols <- which(grepl('Npop',colnames(p2data)))
pops <- lapply(p2data$country[rowindices],function(cn)
  subset(p2data,country==cn)[popcols])

## school a2 ############################
print('school a2')
# what fraction of contacts for 5-19 year olds come from school

out <- c()
illevels <- c()
for(i in 1:sum(rowindices)){
  cn <- p2data$country[which(rowindices)[i]]
  sqmat <- sqmatrices[[i]]
  pop <- pops[[i]]
  schoc <- subset(p2data,country==cn)$schoolA2
  if(nrow(pop)>0&!is.na(schoc)){
    ado <- rowSums(sqmat[2:4,2:4])
    contacts <- sum(ado*pop[2:4])/sum(pop[2:4])
    out <- rbind(out,c(contacts,schoc))
    illevels <- c(illevels, subset(p2data,country==cn)$igroup)
  }
}

frac_school_contacts <- out[,2]/out[,1]
ilevels <- illevels[frac_school_contacts<1]
frac_school_contacts <- frac_school_contacts[frac_school_contacts<1]
a2 <- data.frame(igroup=ilevels,frac_school_contacts=frac_school_contacts)

hic <- fitdistr(subset(a2,!is.na(frac_school_contacts)&igroup=='HIC')$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))
umic <- fitdistr(subset(a2,!is.na(frac_school_contacts)&igroup=='UMIC')$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))
llmic <- fitdistr(subset(a2,!is.na(frac_school_contacts)&igroup%in%c('LMIC','LIC'))$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))
dist2 <- fitdistr(subset(a2,!is.na(frac_school_contacts))$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))

schoolA2_distributions <- data.frame(parameter_name='schoolA2_frac',
                                     igroup=c('LLMIC','UMIC','HIC'),
                                     distribution='betainv',
                                     `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape1']]),
                                     `Parameter 2`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape2']]))


## school a1 #################################
print('school a1')
# what fraction of contacts for 0-4 year olds come from school

out <- c()
illevels <- c()
for(i in 1:sum(rowindices)){
  cn <- p2data$country[which(rowindices)[i]]
  contacts <- p2data$CMaa[which(rowindices)[i]]
  schoc <- subset(p2data,country==cn)$schoolA1
  if(!is.na(contacts)&!is.na(schoc)){
    out <- rbind(out,c(contacts,schoc))
    illevels <- c(illevels, subset(p2data,country==cn)$igroup)
  }
}

frac_school_contacts <- out[,2]/out[,1]
ilevels <- illevels[frac_school_contacts<1]
frac_school_contacts <- frac_school_contacts[frac_school_contacts<1]
a1 <- data.frame(igroup=ilevels,frac_school_contacts=frac_school_contacts)

hic <- fitdistr(subset(a1,!is.na(frac_school_contacts)&igroup=='HIC')$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))
umic <- fitdistr(subset(a1,!is.na(frac_school_contacts)&igroup=='UMIC')$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))
llmic <- fitdistr(subset(a1,!is.na(frac_school_contacts)&igroup%in%c('LMIC','LIC'))$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))
dist1 <- fitdistr(subset(a1,!is.na(frac_school_contacts))$frac_school_contacts,"beta",start=list(shape1=1,shape2=1))

schoolA1_distributions <- data.frame(parameter_name='schoolA1_frac',
                                     igroup=c('LLMIC','UMIC','HIC'),
                                     distribution='betainv',
                                     `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape1']]),
                                     `Parameter 2`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape2']]))


## international tourism ###################
print('international tourism')

data <- as.data.frame(readODS::read_ods(file.path(datapath,'tourism.ods'),sheet=1))
# colnames(data) <- data[1,]
# data <- data[-1,]
for(i in 2:ncol(data)) data[,i] <- as.numeric(data[,i])


data$`International tourism as a share of GDP` <- data$`Tourism as a share of GDP (%)`/100 * data$`International tourism as share of total tourism (%)`


ytd <- as.data.frame(readODS::read_ods(file.path(datapath,'tourism.ods'),sheet=2))

# colnames(ytd) <- ytd[1,]
# ytd <- ytd[-1,]
for(i in 2:ncol(ytd)) ytd[,i] <- as.numeric(ytd[,i])


meltytd <- melt(ytd,id.vars='Country')
meltdata <- melt(data[,c(1:5)],id.vars='Country')
tourismhist <- ggplot() + 
  geom_histogram(data=meltytd,aes(x=value),colour='navyblue',fill='navyblue',bins=30) + 
  geom_histogram(data=meltdata,aes(x=value),colour='white',fill='grey',bins=30) + 
  facet_wrap(~variable,scales='free') +
  theme_bw(base_size = 15) + labs(x='',y='')

dat <- ytd$`International tourist arrivals, YTD change (%)`/100 + 1
dat <- dat[!is.na(dat)]
fit_params <- fitdistr(dat,"lognormal")

# generate values given our fit parameters
x <- seq(0,max(dat),length=10000)
# hst <- hist(dat, breaks=x)
fit <- dlnorm(x, fit_params$estimate['meanlog'], fit_params$estimate['sdlog'])

tourism_distribution <- data.frame(parameter_name='remaining_international_tourism',
                                   igroup=c('all'),
                                   distribution='logninv',
                                   `Parameter 1`=fit_params$estimate[['meanlog']],
                                   `Parameter 2`=fit_params$estimate[['sdlog']])


ggplot(data) + geom_point(aes(x=`International tourism as share of total tourism (%)`,y=`Tourism as a share of GDP (%)`))
ggplot(data) + geom_point(aes(x=`International tourism as share of total tourism (%)`,y=`Food and accommodation services`))

input <- data$`Tourism as a share of GDP (%)`
output <- data$`International tourism as share of total tourism (%)`
nonas <- !is.na(input+output) & output!=0
input <- input[nonas]/100
output <- output[nonas]/100
fit_tourism_data <- function(par){
  # print(par)
  pointiness <- par[1]
  gradient <- par[2]
  intercept <- par[3]
  estoutmean <- pmax(1e-16,pmin(gradient*input + intercept, 1-1e-16))
  alpha <- estoutmean*pointiness
  beta <- pointiness - alpha
  pout <- dbeta(output,alpha,beta)
  -sum(log(pout))
}

optimfit <- optim(par=c(5,3,.1),fit_tourism_data,method='L-BFGS-B',lower=c(0,0,0),upper=c(Inf,10,max(output)))
pointiness <- optimfit$par[1]
betaparams <- optimfit$par[2:3]

tourism_pointiness <- data.frame(parameter_name='tourism_pointiness',
                                 igroup=c('all'),
                                 distribution=NA,
                                 `Parameter 1`=pointiness,
                                 `Parameter 2`=NA)

sec_to_international <- data.frame(parameter_name='sec_to_international',
                                   igroup=c('all'),
                                   distribution=NA,
                                   `Parameter 1`=betaparams[1],
                                   `Parameter 2`=betaparams[2])



## gdp to gni ppp #############################
print('gdp to gni ppp')

indicators <- wb_indicators()
indicators[grepl('GNI',indicators$indicator)&grepl('PPP',indicators$indicator),]
indicators[grepl('GDP',indicators$indicator)&grepl('current',indicators$indicator)&grepl('\\$',indicators$indicator),]

Sys.sleep(2)
gdpdata <- setDT(wb_data("NY.GDP.PCAP.CD",country = "countries_only", start_date = 2018, end_date = 2024))
Sys.sleep(2)
gnipppdata <- setDT(wb_data("NY.GNP.PCAP.PP.CD",country = "countries_only", start_date = 2018, end_date = 2024))
joineddata <- left_join(gdpdata[,.(iso3c,country,date,NY.GDP.PCAP.CD)],gnipppdata[,.(iso3c,country,date,NY.GNP.PCAP.PP.CD)],by=c('iso3c','country','date'))
joineddata <- subset(joineddata,!is.na(NY.GNP.PCAP.PP.CD)&!is.na(NY.GDP.PCAP.CD))
joineddata[,mostrecent:=max(date),by=country]
joineddata <- subset(joineddata,date==mostrecent)
joineddata[,gdp_to_gnippp:=NY.GNP.PCAP.PP.CD/NY.GDP.PCAP.CD]
joineddata <- left_join(joineddata,setDT(incomelevels)[,.(Country.Code,IncomeGroup)],by=c('iso3c'="Country.Code"))

ggplot(joineddata) + geom_point(aes(x=NY.GDP.PCAP.CD,y=gdp_to_gnippp,colour=IncomeGroup)) + theme_bw(base_size = 15) +
  labs(x='GDP pc',y='GDP pc PPP / GDP pc',colour='')

ggplot(joineddata) + geom_histogram(aes(x=gdp_to_gnippp,fill=IncomeGroup),position='dodge') 

hic <- fitdistr(subset(joineddata,!is.na(gdp_to_gnippp)&IncomeGroup=='High income')$gdp_to_gnippp,"gamma")
umic <- fitdistr(subset(joineddata,!is.na(gdp_to_gnippp)&IncomeGroup%in%c('Upper middle income'))$gdp_to_gnippp,"gamma")
llmic <- fitdistr(subset(joineddata,!is.na(gdp_to_gnippp)&IncomeGroup%in%c('Lower middle income','Low income'))$gdp_to_gnippp,"gamma")

gdp_to_gnippp_distributions <- data.frame(parameter_name='gdp_to_gnippp',
                                          igroup=c('LLMIC','UMIC','HIC'),
                                          distribution='gaminv',
                                          `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape']]),
                                          `Parameter 2`=sapply(list(llmic,umic,hic),function(x)1/x$estimate[['rate']]))


## pupil to teacher ratio ###################
print('pupil to teacher ratio')

# allyears_PTRpreprimary=  wb(indicator = "SE.PRE.ENRL.TC.ZS", startdate = 2000, enddate = 2020)
Sys.sleep(1)
allyears_PTRprimary=  setDT(wb_data("SE.PRM.ENRL.TC.ZS",country = "countries_only", 2000, 2024)) 
Sys.sleep(1)
allyears_PTRsecondary=  setDT(wb_data("SE.SEC.ENRL.TC.ZS",country = "countries_only", 2000, 2024))
# allyears_PTRtertiary=  wb(indicator = "SE.TER.ENRL.TC.ZS", startdate = 2000, enddate = 2020)

setnames(allyears_PTRprimary,'SE.PRM.ENRL.TC.ZS','value')
setnames(allyears_PTRsecondary,'SE.SEC.ENRL.TC.ZS','value')
allyears_PTRprimary[!is.na(value),maxyear:=max(date),by=iso3c]
allyears_PTRsecondary[!is.na(value),maxyear:=max(date),by=iso3c]
PTRprimary  <- subset(allyears_PTRprimary,date==maxyear)
PTRsecondary  <- subset(allyears_PTRsecondary,date==maxyear)

ptrdt <- setDT(do.call(rbind,list(PTRprimary,PTRsecondary)))
ptrdt[,meanptr:=mean(value,na.rm=T),by=iso3c]
meanptrs <- unique(ptrdt[,.(country,iso3c,meanptr)])
setorder(meanptrs,meanptr)
# meanptrs <- subset(meanptrs,!is.na(countrycode(country,origin='country.name',destination = 'iso3c')))
meanptrs$country <- factor(meanptrs$country,levels=meanptrs$country)
countries <- c('United Kingdom')
ggplot(meanptrs) + 
  geom_point(aes(x=1:length(country),y=meanptr),colour='navyblue') + 
  coord_flip() + 
  geom_point(data=with(meanptrs,data.frame(x=which(country%in%countries),y=meanptr[which(country%in%countries)])),
             aes(x=x,y=y),colour='orange',size=4) + 
  theme_bw(base_size = 15) +
  labs(x='',y='Pupil-teacher ratio')+
  scale_x_continuous(expand=c(0,0),
                     breaks = seq(1,nrow(meanptrs),by=2),
                     labels = meanptrs$country[seq(1,nrow(meanptrs),by=2)],
                     sec.axis = dup_axis(breaks=seq(2,nrow(meanptrs),by=2),labels = meanptrs$country[seq(2,nrow(meanptrs),by=2)])) +
  theme(axis.text=element_text(size=10)) -> ptrp
ggsave(ptrp,filename='store/pupil_teacher_ratio.png',width=8,height=9)

meanptrs <- left_join(meanptrs,incomelevels,by=c('iso3c'="Country.Code"))

hic <- fitdistr(subset(meanptrs,!is.na(meanptr)&IncomeGroup=='High income')$meanptr,"gamma")
umic <- fitdistr(subset(meanptrs,!is.na(meanptr)&IncomeGroup%in%c('Upper middle income'))$meanptr,"gamma")
llmic <- fitdistr(subset(meanptrs,!is.na(meanptr)&IncomeGroup%in%c('Lower middle income','Low income'))$meanptr,"gamma")

ptr_distributions <- data.frame(parameter_name='pupil_teacher_ratio',
                                igroup=c('LLMIC','UMIC','HIC'),
                                distribution='gaminv',
                                `Parameter 1`=sapply(list(llmic,umic,hic),function(x)x$estimate[['shape']]),
                                `Parameter 2`=sapply(list(llmic,umic,hic),function(x)1/x$estimate[['rate']]))

ggplot(meanptrs) + geom_histogram(aes(x=meanptr,fill=IncomeGroup),position='dodge') 


## contacts #########################
print('contacts')

source('workplace_related_contacts.R')

## compliance ###########################

# files <- list.files(yougovdata)
# files <- files[files!='tmp']
# csv_countries <- sapply(files[grepl('csv',files)],function(x) strsplit(x,'\\.')[[1]][1])
# ygdata <- list()
# for(i in files){
#   cn <- strsplit(i,'\\.')[[1]][1]
#   
#   if(grepl('csv',i)){
#     cndata <- read.csv(file.path(yougovdata,i))
#     cndata$country <- cn
#     ygdata[[cn]] <- cndata[,colnames(cndata)%in%c('country','qweek','i11_health','weight')]
#   }else if(!cn %in% csv_countries){
#     tmppath <- file.path(yougovdata,'tmp')
#     if(dir.exists(tmppath))
#       system(paste0('rm ',tmppath,'/*'))
#     else
#       system(paste0('mkdir ',tmppath))
#     unzip(file.path(yougovdata,i),exdir=tmppath)
#     print(list.files(tmppath,full.names=T))
#     cndata <- read.csv(list.files(tmppath,full.names=T), fileEncoding="latin1")
#     write.csv(cndata,file.path(yougovdata,paste0(cn,'.csv')))
#   }
# }
# 
# ygdf <- setDT(do.call(rbind,ygdata))
# ygdf <- subset(ygdf,!i11_health%in%c(' ','Not sure')) %>%
#   mutate(compliance=as.numeric(plyr::mapvalues(i11_health,from=c('Very willing',
#                                                     'Somewhat willing',
#                                                     'Neither willing nor unwilling',
#                                                     'Somewhat unwilling',
#                                                     'Very unwilling'),
#                                     to=c(1,.75,.5,.25,0))))
# ygdf[,week:=as.numeric(gsub('week ','',qweek)),by=qweek]
# ygdf[,sum(compliance*weight,na.rm=T)/sum(weight)]
# meanvals <- ygdf[,sum(compliance*weight,na.rm=T)/sum(weight),by=country]
# # uk: 87%
# hist(meanvals$V1)
# range(meanvals$V1)
# ygdf[,.N,by=country]
# ygdf[,mean(compliance,na.rm=T),by=week]
# sapply(ygdf,class)
# table(ygdf$i11_health)

## end ############################################

(parameter_distributions <- rbind(#internet_distributions,
                                  tourism_distribution,
                                  labsh_distributions,
                                  gdp_to_gnippp_distributions,
                                  Hmax_distributions,
                                  bmi_distributions,
                                  # pt_distributions,
                                  tourism_pointiness,
                                  sec_to_international,
                                  ptr_distributions,
                                  school1_distributions,
                                  school2_distributions,
                                  work_frac_distributions,
                                  hospitality_frac_distributions,
                                  hospitality_age_distributions,
                                  workforce_in_place
))

write.csv(parameter_distributions,file.path(datapath,'parameter_distributions.csv'),row.names = F)

#######################################

# vsl plot

Sys.sleep(2)
gdpdata2 <- setDT(wb_data("NY.GDP.PCAP.CD",country = "countries_only", start_date = 2018, end_date = 2024))
Sys.sleep(5)
gnipppdata2 <- setDT(wb_data("NY.GDP.PCAP.PP.CD",country = "countries_only", start_date = 2018, end_date = 2024))
joineddata <- left_join(gdpdata2[,.(iso3c,country,date,NY.GDP.PCAP.CD)],gnipppdata2[,.(iso3c,country,date,NY.GDP.PCAP.PP.CD)],by=c('iso3c','country','date'))
joineddata <- subset(joineddata,!is.na(NY.GDP.PCAP.PP.CD)&!is.na(NY.GDP.PCAP.CD))
joineddata[,mostrecent:=max(date),by=country]
joineddata <- subset(joineddata,date==mostrecent)
joineddata[,gdp_to_gnippp:=NY.GDP.PCAP.PP.CD/NY.GDP.PCAP.CD]
joineddata <- left_join(joineddata,setDT(incomelevels)[,.(Country.Code,IncomeGroup)],by=c('iso3c'="Country.Code"))


joineddata[,gnippp:=gdp_to_gnippp*NY.GDP.PCAP.CD]
usagdppc <- with(subset(joineddata,iso3c=='USA'),NY.GDP.PCAP.CD)
yint <- with(subset(joineddata,iso3c=='USA'),(10^6*10*(NY.GDP.PCAP.CD/usagdppc)^1.)/NY.GDP.PCAP.CD)
joineddata[,minvsl:=20*NY.GDP.PCAP.CD]
joineddata[,basicvsl:=max(minvsl,10^6*10*(NY.GDP.PCAP.CD/usagdppc)^1.),by=iso3c]
joineddata[,vsle:=max(minvsl,10^6*10*(NY.GDP.PCAP.CD/usagdppc)^1.5),by=iso3c]
joineddata[,vslpppe:=max(minvsl,10^6*10*(gnippp/usagdppc)^1.5),by=iso3c]
joineddata[,vslppp:=max(minvsl,10^6*10*(gnippp/usagdppc)^1.),by=iso3c]
subset(joineddata,iso3c=="COD")

p11 <- ggplot(joineddata,aes(x=NY.GDP.PCAP.CD,y=basicvsl/NY.GDP.PCAP.CD)) + 
  geom_hline(yintercept=yint,linewidth=2,colour='grey') +
  annotate('rect',ymin=10^(1.5),ymax=10^(2),xmin=-Inf,xmax=Inf,fill='goldenrod',alpha=.25) +
  geom_point(aes(colour=IncomeGroup),show.legend=F) + theme_bw(base_size = 13) +
  labs(x='GDP pc',y='VSL/GDPpc',title='MER, e=1',colour='') 
p12 <- ggplot(joineddata,aes(x=(NY.GDP.PCAP.CD),y=vslppp/NY.GDP.PCAP.CD)) + 
  geom_hline(yintercept=yint,linewidth=2,colour='grey') +
  annotate('rect',ymin=10^(1.5),ymax=10^(2),xmin=-Inf,xmax=Inf,fill='goldenrod',alpha=.25) +
  geom_point(aes(colour=IncomeGroup),show.legend=F) + theme_bw(base_size = 13) +
  labs(x='GDP pc',y='VSL/GDPpc',title='PPP, e=1',colour='') 
p21 <- ggplot(joineddata,aes(x=(NY.GDP.PCAP.CD),y=vsle/NY.GDP.PCAP.CD)) + 
  geom_hline(yintercept=yint,linewidth=2,colour='grey') +
  annotate('rect',ymin=10^(1.5),ymax=10^(2),xmin=-Inf,xmax=Inf,fill='goldenrod',alpha=.25) +
  geom_point(aes(colour=IncomeGroup),show.legend=F) + theme_bw(base_size = 13) +
  labs(x='GDP pc',y='VSL/GDPpc',title='MER, e=1.5',colour='') 
p22 <- ggplot(joineddata,aes(x=(NY.GDP.PCAP.CD),y=vslpppe/NY.GDP.PCAP.CD)) + 
  geom_hline(yintercept=yint,linewidth=2,colour='grey') +
  annotate('rect',ymin=10^(1.5),ymax=10^(2),xmin=-Inf,xmax=Inf,fill='goldenrod',alpha=.25) +
  geom_point(aes(colour=IncomeGroup),show.legend=F) + theme_bw(base_size = 13) +
  labs(x='GDP pc',y='VSL/GDPpc',title='PPP, e=1.5',colour='') 


plot_grid(p11,p21,p12,p22,ncol = 2, nrow = 2)







cod <- subset(gnipppdata,iso3c%in%c('COD','USA')&date==2023)
vsl_usa <- 10e6
gnipc_cod <- subset(cod,iso3c=='COD')$NY.GNP.PCAP.PP.CD
gnipc_usa <- subset(cod,iso3c=='USA')$NY.GNP.PCAP.PP.CD
vsl_cod <- max((gnipc_cod/gnipc_usa)^1.5 * vsl_usa, 20*gnipc_cod)






