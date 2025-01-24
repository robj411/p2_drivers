

## Packages
library(MASS)
library(splines)
library(data.table)
library(lubridate)
library(magrittr)
library(ggplot2)

# using data from https://github.com/jarvisc1/cmix_post_pandemic
# checkout a6c98ab
# many columns were omitted since
cmixpath <- '../../../cmix_post_pandemic/'
datapath <- '../../data/'

## functions ############################
betavals <- function(vals){
  meanv <- mean(vals)
  sdv <- sd(vals)
  # aoveraplusb <- meanv
  # meanv*(a+b) = a
  # meanv*b = a - meanv*a
  # b = a/meanv - a
  # sdv^2 = ab/((a+b)^2*(a+b+1))
  # sdv^2*(a/meanv+1) = (1/meanv - 1)/(1/meanv^2)
  a <- ((1/meanv - 1)/(1/meanv^2*sdv^2) - 1)*meanv
  b <- a/meanv - a
  # x <- rbeta(10000,a,b)
  c(a,b)
}

## data ################################

popdata <- subset(read.csv(file.path(datapath,'IHME_GBD_2019_POP_2019_Y2020M10D15.csv')),sex_name=='both')
agecodes <- list('<20 years',c('65 to 74',"75 plus"),"All Ages")
countrycodes <- c(UK='United Kingdom',NL='Netherlands',CH='Switzerland',BE='Belgium')
popsizes <- sapply(countrycodes,function(x){
  pop <- subset(popdata,location_name==x)
  pop1 <- subset(pop,age_group_name%in%agecodes[[1]])$val
  pop3 <- sum(subset(pop,age_group_name%in%agecodes[[2]])$val)
  pop2 <- subset(pop,age_group_name%in%agecodes[[3]])$val - pop1 - pop3
  c(`Under 20`=pop1,`Working age`=pop2,`Over 64`=pop3)
})

## maps #######################################⎄m

occmap <- read.csv(file.path(datapath,'isco.csv'))
occmap$major_label <- tolower(trimws(occmap$major_label))
occmap$description <- tolower(trimws(occmap$description))
occmap$sub_major_label <- tolower(trimws(occmap$sub_major_label))
occmap$minor_label <- tolower(trimws(occmap$minor_label))
subsm <- subset(occmap,ISCO_version=='ISCO-88')



# Load data ---------------------------------------------------------------
dt <- qs::qread(file.path(cmixpath,'data/wrapup_part_cnts.qs'))
# subset to this survey
dt <- dt[survey_round == 1000]

# store participant age info
partage <- dt[,.(part_wave_uid,part_age,part_age_est_min,part_age_est_max)]
# assign daedalus age groups
partage[part_age_est_max<20,agegrp:='Under 20']
partage[part_age_est_max<65&part_age_est_min>17,agegrp:='Working age']
partage[part_age<65&part_age_est_min>17,agegrp:='Working age']
partage[part_age>=65,agegrp:='Over 64']
partage[part_age_est_min>=65,agegrp:='Over 64']
subset(partage,is.na(agegrp))

dt[,dae_age:=cut(as.numeric(part_age_est_max),breaks=c(0,5,18,65,as.numeric(max(dt$part_age[!is.na(dt$part_age)]))+1),right=F)]
(dae_ages <- levels(setorder(dt,part_age_est_min)$dae_age))
# save for later
saveRDS(partage,'store/partage.Rds')


colnames(dt)[grepl('full',colnames(dt))]
inc <- copy(dt[,.(n_cnt_work,income_3,part_income)])
levels(as.factor(inc$income_3))
mappart<-c("£10,000 - £14,999"=3, "£100,000 or more"=9,  "£15,000 - £19,999"=4, "£20,000 - £24,999"=4, "£25,000 - £34,999"=5, "£35,000 - £44,999"=6, "£45,000 - £54,999"=7, "£5,000 - £9,999"=2,   "£55,000 - £99,999"=8, "under £5,000"=1)
map3 <- c("€ 0 - € 499"=1,       "€ 1 000 - € 1 499"=3, "€ 1 500 - € 2 499"=4, "€ 10 000 or more"=9,  "€ 2 500 - € 3 499"=5, "€ 3 500 - € 4 999"=6, "€ 5 000 - € 6 999"=7, "€ 500 - € 999"=2,     "€ 7 000 - € 9 999"=8)
inc[,income_3:=map3[income_3]]
inc[,part_income:=mappart[part_income]]
inc$income <- inc$part_income
inc$income[is.na(inc$income)] <- inc$income_3[is.na(inc$income)]
inc <- subset(inc,tolower(income)!='prefer not to answer')
with(subset(dt,part_attend_work_yesterday=='yes'),table(n_cnt_work,income_3))
with(subset(dt,part_attend_work_yesterday=='yes'),table(n_cnt_work,part_income))

## public transport ###############################
# relationship between mode share and # pt contacts (not used)

# others <- dt[,list(pt=mean(n_cnt_public_transport)),by=.(country)]
# pttab <- xlsx::read.xlsx(file.path(datapath,'mode_shares.xlsx'),sheetIndex = 2)
# pt4 <- subset(pttab,Country.name%in%c('Belgium','Netherlands','Switzerland','United Kingdom'))
# ptshare <- apply(pt4,1,function(x)sum(as.numeric(x[c(2,3)])))
# others$ptshare <- ptshare[c(1,3,2,4)]
# (pt <- ggplot(others,aes(x=ptshare,y=pt)) + geom_smooth(method=glm,formula=y~ns(x,df=1)) +
#     geom_label(aes(label=country),size=5) +
#     theme_bw(base_size=15) +
#     labs(x='Public transport mode share',y='Public transport contacts'))
# ggsave(pt,filename='store/ptplot.png',width=6,height=4)  
# summary(glm(pt~ptshare,data=others))

## summarise fractions spent in places ##########################
# not used

summ <- dt[,list(total=mean(n_cnt),home=mean(n_cnt_home),school=mean(n_cnt_school),work=mean(n_cnt_work),other=mean(n_cnt_other)),by=.(country,dae_age)]
summ2 <- dt[,list(total=mean(n_cnt),home=mean(n_cnt_home),school=mean(n_cnt_school),work=mean(n_cnt_work),other=mean(n_cnt_other)),by=.(country,survey_round)]
colnames(summ2)[2] <- 'dae_age'
summ2$dae_age <- as.character(summ2$dae_age)
summ2$dae_age <- 'everyone'
summ <- rbind(subset(summ,!is.na(dae_age)),summ2)

summ[,total:=home+work+school+other]
summ[,home:=home/total]
summ[,work:=work/total]
summ[,school:=school/total]
summ[,other:=other/total]

meltsum <- reshape2::melt(summ,id.vars=c('country','dae_age'))
(p1 <- ggplot(subset(meltsum,variable!='total')) + 
    theme_bw(base_size=15) +
    labs(x='',y='Percent of contacts',fill='') +
    geom_bar(aes(x=variable,y=value*100,fill=dae_age),position='dodge',stat='identity') +
    facet_wrap(~country))
ggsave(p1,filename='store/fractions.png',width=8,height=5)

## get school contacts as a fraction of total contacts #########################

allcs <- betavals(subset(meltsum,variable=='school'&dae_age=='[0,5)')$value)
schoolfracp <- ggplot() + geom_histogram(aes(x=rbeta(100000,allcs[1],allcs[2]),y=..density..),fill='navyblue') + 
  scale_x_continuous(limits=c(0,1)) +
  theme_bw(base_size = 15) +
  labs(x='Fraction contacts from school (0-4)',y='')
ggsave(schoolfracp,filename='store/school1frac.png',width=6,height=5)
(school1_distributions <- data.frame(parameter_name='school1_frac',
                                     igroup=c('all'),
                                     distribution='betainv',
                                     `Parameter 1`=allcs[1],
                                     `Parameter 2`=allcs[2]))

allcs <- betavals(subset(meltsum,variable=='school'&dae_age=='[5,18)')$value)
schoolfracp <- ggplot() + geom_histogram(aes(x=rbeta(100000,allcs[1],allcs[2]),y=..density..),fill='navyblue') + 
  scale_x_continuous(limits=c(0,1)) +
  theme_bw(base_size = 15) +
  labs(x='Fraction contacts from school (5-19)',y='')
ggsave(schoolfracp,filename='store/school2frac.png',width=6,height=5)

(school2_distributions <- data.frame(parameter_name='school2_frac',
                                     igroup=c('all'),
                                     distribution='betainv',
                                     `Parameter 1`=allcs[1],
                                     `Parameter 2`=allcs[2]))

## frac work contacts #######################
# use time at work to extend frac of contacts at work

tus <- xlsx::read.xlsx(file.path(datapath,'OECD_TIME_USE.xlsx'),sheetIndex = 1, startRow = 3)
tus <- tus[!is.na(tus$Other),-c(2,ncol(tus))]
for(i in 2:ncol(tus))
  tus[,i] <- tus[,i]/1440*100
setDT(tus)
setorder(tus,Paid.work.or.study)
subset(tus,Measure%in%c('Belgium','Netherlands','United Kingdom'))
subset(meltsum,variable=='work'&dae_age=='[18,65)')

countrycodes <- c(UK='United Kingdom',NL='Netherlands',CH='Switzerland',BE='Belgium')
df <- data.frame(do.call(rbind,lapply(1:length(countrycodes),function(x){
  nm <- countrycodes[x]
  cd <- names(countrycodes)[x]
  workcontacts <- subset(meltsum,variable=='work'&dae_age=='[18,65)'&country==cd)$value
  schoolA2contacts <- subset(meltsum,variable=='school'&dae_age=='[5,18)'&country==cd)$value
  schoolA1contacts <- subset(meltsum,variable=='school'&dae_age=='[0,5)'&country==cd)$value
  timeuse <- subset(tus,Measure%in%nm)$Paid.work.or.study
  vals <- c(cd,workcontacts,schoolA2contacts,schoolA1contacts,timeuse)
  if(length(vals)==5)
    return(vals)
  else
    return(NULL)
})))
colnames(df) <- c('country','work','schoolA1','schoolA2','timeuse')  
for(i in 2:ncol(df)) df[,i] <- as.numeric(df[,i])

vals <- predict(glm(work~timeuse,data=df),newdata=data.frame(timeuse=tus$Paid.work.or.study))
vals <- tus$Paid.work.or.study * with(df,mean(work)/mean(timeuse))
weights <- sapply(vals,function(x)sum(abs(vals-x)))^2

allcs <- betavals(sample(vals,1000,replace = T,prob=weights))
workfracp <- ggplot() + geom_histogram(aes(x=rbeta(100000,allcs[1],allcs[2]),y=..density..),fill='navyblue') + 
  scale_x_continuous(limits=c(0,1)) +
  theme_bw(base_size = 15) +
  labs(x='Fraction contacts from work',y='')
ggsave(workfracp,filename='store/workfrac.png',width=6,height=5)

(work_frac_distributions <- data.frame(parameter_name='work_frac',
                                       igroup=c('all'),
                                       distribution='betainv',
                                       `Parameter 1`=allcs[1],
                                       `Parameter 2`=allcs[2]))


## hospitality contacts as fraction of non-work/school contacts #################

con  <- dt[,list(hosp=mean(n_cnt_shop+n_cnt_bar_rest+n_cnt_salon),total=mean(n_cnt-n_cnt_work-n_cnt_school)),
           by=.(country,dae_age)]
con <- subset(con,!is.na(dae_age))

con[,hosp:=hosp/total]
meltcon <- reshape2::melt(con,id.vars=c('country','dae_age'))
(p2 <- ggplot(subset(meltcon,variable!='total')) + 
    theme_bw(base_size=15) +
    labs(x='',y='Hospitality as percent of non-work/school contacts',fill='') +
    geom_bar(aes(x=country,y=value*100,fill=dae_age),position='dodge',stat='identity'))# +
# facet_wrap(~country))
ggsave(p2,filename='store/confractions.png',width=8,height=6)

ages <- levels(subset(con,dae_age!='everyone')$dae_age)
betatab <- sapply(ages[ages!='everyone'],function(x){
  vals <- subset(con,dae_age==x)$hosp
  betavals(vals)
})
samples <- as.data.frame(do.call(rbind,lapply(colnames(betatab),function(x){
  xt <- betatab[,x]
  a <- max(1,xt[1]-1)
  b <- max(1,xt[2]-1)
  xv <- rbeta(10000,a,b)
  cbind(xv,x)
})))
samples$xv <- as.numeric(samples$xv)
hospfrac <- ggplot(samples) + geom_histogram(aes(x=as.numeric(xv),y=..density..,fill=x),show.legend = F) +
  facet_wrap(~x) + theme_bw(base_size = 15) +
  scale_x_continuous(limits=c(0,1)) +
  labs(x='Fraction non-work/school contacts from hospitality',y='')
ggsave(hospfrac,filename='store/hospfrac.png',width=6,height=5)


paste0(sapply(colnames(dt)[grepl('n_cnt_',colnames(dt))],function(x)gsub('n_cnt_','',x)),collapse=', ')

(hospitality_frac_distributions <- data.frame(parameter_name=paste0('hospitality',1:4,'_frac'),
                                              igroup=c('all'),
                                              distribution='betainv',
                                              `Parameter 1`=betatab[1,],
                                              `Parameter 2`=betatab[2,]))

## contact data ####################

cnts <- setDT(readRDS('store/contacts1000.Rds'))
cntocc <- setDT(readRDS('store/contactsoccupation.Rds'))

untab <- unique(cntocc[,.(part_wave_uid,occcode2,n_id,n_ind,n_tot,n_WW,n_CW,n_under18,n_workingage,n_65plus,country)])
untab[,occupation:=subsm$major_label[match(occcode2,subsm$major)],by=.(occcode2,country)]
untab[is.na(occupation),occupation:=subsm$sub_major_label[match(occcode2,subsm$sub_major)],by=.(occcode2,country)]

# mfrac <- melt(fracage,id.vars=c('occupation'))
fracage <- untab[,.(under18=mean(n_under18),workingage=mean(n_workingage),`65plus`=mean(n_65plus)),by=.(occupation,country)]
mfrac <- setDT(melt(fracage,id.vars=c('country','occupation')))
setorder(mfrac,value)
mfrac[,occupation:=factor(occupation,levels=unique(occupation))]
(page <- ggplot(mfrac) +
    geom_bar(aes(x=occupation,y=value,fill=variable,group=country),stat='identity') + #,position='dodge'
    labs(x='',y='Number of contacts',fill='') + facet_wrap(~country,scale='free_x') +
    coord_flip() + theme_bw(base_size=15) + theme(legend.position = 'top'))
ggsave(page,filename = paste0('store/age.png'),height=6,width=15)


## mixing proportions in hospitality ################################

partage <- readRDS('store/partage.Rds')
concnt <- partage[subset(cnts,cnt_bar_rest>0),,on='part_wave_uid']
concnt[,n_under18:=sum(cnt_age_est_max<18),by=part_wave_uid]
concnt[,n_workingage:=sum(cnt_age_est_max<65&cnt_age_est_min>17),by=part_wave_uid]
concnt[,n_65plus:=sum(cnt_age_est_min>=65),by=part_wave_uid]
concnt <- subset(concnt,!is.na(part_age_est_max))
concnt[,agegrp:=cut(as.numeric(part_age_est_max),breaks=c(0,5,18,65,as.numeric(max(dt$part_age[!is.na(dt$part_age)]))+1),right=F)]
mcon <- setDT(melt(concnt,id.vars=c('country','agegrp'),measure.vars = c('n_under18','n_workingage','n_65plus')))
mcon[,total:=sum(value),by=.(agegrp,country)]
mcon[,agetotal:=sum(value),by=.(agegrp,country,variable)]
mcon[,frac:=agetotal/total,by=.(agegrp,country,variable)]
mcon$variable <- factor(mcon$variable,levels=unique(mcon$variable),labels=c('Under 18','Working age','Over 64'))
mcon$agegrp <- factor(mcon$agegrp,levels=rev(levels(mcon$agegrp)))
snames <- unique(mcon$agegrp)
onames <- unique(mcon$variable)
confrac <- unique(mcon[,.(variable,agegrp,country,frac)])
(pconfrac <- ggplot(confrac) + 
    # geom_tile(data=mcon,aes(x=factor(variable,levels=onames),y=factor(agegrp,levels=rev(snames))),fill='white', show.legend = FALSE) +
    facet_wrap(~country) +
    scale_y_discrete(expand=c(0,0)) +
    geom_tile(aes(x=variable,y=agegrp,fill=round(frac*100)), show.legend = FALSE) +
    geom_text(aes(x=variable,y=agegrp,label=round(frac*100)),size=8) +
    # geom_text(data=mcon,aes(x=variable,y=factor(agegrp,levels=rev(snames[1:length(onames)])),label=contacts),size=8) +
    theme_bw(base_size = 15) +
    scale_x_discrete(expand=c(0,0),position='top') +
    scale_fill_gradient(low="white", high="goldenrod1") +
    labs(x='',y=''))# +
# annotate("segment", x = 0.5, xend = length(cnames)+.5, y = 1.5, yend = 1.5, colour = "black", size=1, alpha=1)
ggsave(pconfrac,filename = paste0('store/conagefrac.png'),height=7,width=7)

dird <- confrac[,mean(frac),by=.(agegrp,variable)]


(hospitality_age_distributions <- data.frame(parameter_name=paste0('hospitality_age',1:4),
                                             igroup=c('all'),
                                             distribution=NA,
                                             `Parameter 1`=sapply(dae_ages,function(x)subset(dird,agegrp==x&variable=='Working age')$V1),
                                             `Parameter 2`=sapply(dae_ages,function(x)subset(dird,agegrp==x&variable=='Over 64')$V1)))


## who do workers contact ###########################

# contacts from workers
hospworkers <- c("Housekeeping and restaurant services workers","Other personal services workers","Sales and services elementary occupations")
hospworkersids <- dt[(part_occupation%in%hospworkers|nl01occr%in%hospworkers|ch01occr%in%hospworkers)&part_attend_work_yesterday=='yes',part_wave_uid]

concnt <- cnts[part_wave_uid%in%hospworkersids&cnt_work>0,]
concnt[cnt_age_est_max<20,age_est:='Under 20']
concnt[cnt_age_est_max<65&cnt_age_est_min>17,age_est:='Working age']
concnt[cnt_age_est_min>=65,age_est:='Over 64']
workercontacts <- concnt[!is.na(age_est),sum(cnt_work),by=.(age_est,country)]
workercontacts[,total:=sum(V1),by=country]
workercontacts <- workercontacts[dt[country%in%workercontacts$country,sum(dae_age=='[18,65)',na.rm=T),by=country],,on='country']
workercontacts[,contactspp:=total/i.V1]
workercontacts[,agefrac:=V1/total]
workercontacts[,agecontactspp:=V1/i.V1]


## matching numbers in hospitality ########################
# do the contacts reported between under 20s and over 64s match as reported by each group?
# no

partage[,cntry:=substr(part_wave_uid,1,2)]
partage[,N:=.N,by=.(agegrp,cntry)]

concnt <- partage[cnts,,on='part_wave_uid']
concnt <- subset(concnt,!is.na(agegrp))
concnt[,age_est:='']
concnt[cnt_age_est_max<20,age_est:='Under 20']
concnt[cnt_age_est_max<65&cnt_age_est_min>17,age_est:='Working age']
concnt[cnt_age_est_min>=65,age_est:='Over 64']
concnt <- subset(concnt,age_est!='')

mcon <- unique(concnt[cnt_type!='work',sum(cnt_salon|cnt_bar_rest|cnt_shop|cnt_public_transport)/N,by=.(age_est,agegrp,country)])
mcon[,otherpop:=popsizes[age_est,as.character(country)],by=.(age_est,country)]
mcon[,mypop:=popsizes[agegrp,as.character(country)],by=.(country)]
mcon[,contactperother:=V1*mypop/otherpop]
mcon2 <- mcon[agegrp%in%c('Over 64','Under 20')&age_est%in%c('Over 64','Under 20')&(age_est!=agegrp),]
mcon2[,V1:=rev(V1),by=country]
(pagecon <- ggplot(mcon2) + 
    geom_abline(intercept=0,slope=1,colour='grey',linewidth=2) +
    geom_point(aes(y=contactperother,x=V1,colour=country,shape=age_est),size=4)+
    theme_bw(base_size = 15) +
    labs(y='As object',x='As subject',colour='Country',shape='Contacts per'))
ggsave(pagecon,filename='store/ageconmismatch.png',width=6,height=5)

## contacts to working age ##############################
# can we distinguish contacts with customers who are working age from contacts who are working and working age?
# no

wa_col_contacts <- unique(mcon[age_est=='Working age'&agegrp!='[0,5)',])
wa_col_contacts[,contactpp:=V1]
wa_col_contacts[agegrp=='Under 20',age_est:='Under 20']
wa_col_contacts[agegrp=='Over 64',age_est:='Over 64']

subset(wa_col_contacts,age_est=='Working age')

joined <- left_join(workercontacts[,.(age_est,country,agecontactspp)],
                    wa_col_contacts[,.(age_est,country,contactpp)],by=c('country','age_est'))

(phosp <- ggplot(joined) +
    geom_abline(intercept=0,slope=1,colour='grey',linewidth=2) +
    geom_point(aes(y=agecontactspp,x=contactpp,colour=age_est,shape=country),size=4)+
    theme_bw(base_size = 15) +
    labs(y='Reported by working working-age people',x='Hospitality contacts reported with working-age people',colour='Age'))
ggsave(phosp,filename='store/ageconworker.png',width=6,height=5)

joined[,.(mean(agecontactspp),mean(contactpp)),by=age_est]

(phosp <- ggplot(joined[,.(mean(agecontactspp),mean(contactpp)),by=age_est]) +
    geom_abline(intercept=0,slope=1,colour='grey',linewidth=2) +
    geom_point(aes(y=V1,x=V2,colour=age_est),size=4)+
    theme_bw(base_size = 15) +
    labs(y='Reported by working working-age people',x='Hospitality contacts reported with working-age people',colour='Age'))

