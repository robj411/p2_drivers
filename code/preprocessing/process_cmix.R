

## Packages
library(MASS)
library(splines)
library(data.table)
library(lubridate)
library(magrittr)
library(viridis)
library(ggplot2)
library(xlsx)

# using data from https://github.com/jarvisc1/cmix_post_pandemic
# checkout a6c98ab
# many columns were omitted since
cmixpath <- '../../../cmix_post_pandemic/'
datapath <- '../../data/'


# Load data ---------------------------------------------------------------
dt <- qs::qread(file.path(cmixpath,'data/wrapup_part_cnts.qs'))
partage <- dt[survey_round == 1000,.(part_wave_uid,part_age,part_age_est_min,part_age_est_max)]
partage[part_age_est_max<20,agegrp:='Under 20']
partage[part_age_est_max<65&part_age_est_min>17,agegrp:='Working age']
partage[part_age<65&part_age_est_min>17,agegrp:='Working age']
partage[part_age>=65,agegrp:='Over 64']
partage[part_age_est_min>=65,agegrp:='Over 64']
# subset(partage,is.na(agegrp))

saveRDS(partage,'store/partage.Rds')
dt <- dt[survey_round == 1000]

## contact data ####################

# Load data ---------------------------------------------------------------
cnts <- qs::qread(file.path(cmixpath,'data/wrapup_contacts.qs'))

cnts <- cnts[survey_round == 1000]
country_levs <- c("all", "uk", "be", "nl", "ch")
country_labs <- c("All", "UK", "BE", "NL", "CH")
cnts[, country := factor(country, levels = country_levs, labels = country_labs)]

saveRDS(cnts,'store/contacts1000.Rds')

## occupations #######################

occcols <- colnames(dt)[grepl('occ',colnames(dt))]
sapply(occcols,function(x)sum(!is.na(dt[[x]])))
nldt <- dt[!is.na(nl01occr)&part_attend_work_yesterday == 'yes',]
bedt <- dt[!is.na(be02occhi)&part_attend_work_yesterday == 'yes',]
chdt <- dt[!is.na(ch01occr)&part_attend_work_yesterday == 'yes',]
ukdt <- dt[!is.na(part_occupation)&part_attend_work_yesterday == 'yes',]


## get maps ##########################
occmap <- read.csv(file.path(datapath,'isco.csv'))
occmap$major_label <- tolower(trimws(occmap$major_label))
occmap$description <- tolower(trimws(occmap$description))
occmap$sub_major_label <- tolower(trimws(occmap$sub_major_label))
occmap$minor_label <- tolower(trimws(occmap$minor_label))
occuptions <- unique(c(occmap$major_label,occmap$description,occmap$sub_major_label,occmap$minor_label))

om <- subset(occmap,ISCO_version%in%'ISCO-88')

# soc2isco1 <- xlsx::read.xlsx('../data/soc2010toisco08.xls',sheetIndex = 3)
# soc2isco50 <- xlsx::read.xlsx('../data/soc2010toisco08.xls',sheetIndex = 4)
# soc2isco40 <- xlsx::read.xlsx('../data/soc2010toisco08.xls',sheetIndex = 5)
# soc2iscou <- xlsx::read.xlsx('../data/soc2010toisco08.xls',sheetIndex = 6)

# isco88 to isco08
# https://www.ilo.org/public/english/bureau/stat/isco/isco08/
# added map 1213->1239
isco8to0 <- setDT(xlsx::read.xlsx(file.path(datapath,'index08-draft.xlsx'),sheetIndex = 1))
isco8to0[,isco1:=substr(ISCO.88,1,1),by=ISCO.88]
isco8to0[,isco2:=substr(ISCO.88,1,2),by=ISCO.88]
isco8to0[,isco3:=substr(ISCO.88,1,3),by=ISCO.88]

isco2soc <- xlsx::read.xlsx(file.path(datapath,'lookup_soc2010_to_isco_v3.xlsx'),sheetIndex = 1)
colnames(isco2soc)[2:1] <- c('ISCO.08','SOC.2010')

# isco08 to soc (ONS)
# https://www.ons.gov.uk/methodology/classificationsandstandards/standardoccupationalclassificationsoc/soc2020/classifyingthestandardoccupationalclassification2020soc2020totheinternationalstandardclassificationofoccupationsisco08
isco2soc <- xlsx::read.xlsx(file.path(datapath,'onetsocconversion.xlsx'),sheetIndex = 4,startRow = 3)
newisco <- data.frame(ISCO.08=subset(occmap,ISCO_version=='ISCO-08'&!unit%in%isco2soc$ISCO.08)$unit,
                      SOC.2010=subset(occmap,ISCO_version=='ISCO-08'&!unit%in%isco2soc$ISCO.08)$unit,
                      Proportion=1)
tisco2soc <- setDT(rbind(isco2soc,newisco))
tisco2soc[,isco1:=substr(ISCO.08,1,1),by=ISCO.08]
tisco2soc[,isco2:=substr(ISCO.08,1,2),by=ISCO.08]
tisco2soc[,isco3:=substr(ISCO.08,1,3),by=ISCO.08]

# soc to sic
# https://www.ons.gov.uk/employmentandlabourmarket/peopleinwork/employmentandemployeetypes/adhocs/008102numbersofworkersinoccupationsbyindustryofemployment
soc2sic <- xlsx::read.xlsx(file.path(datapath,'soc4xsic12supfinal.xls'),sheetIndex = 2,startRow = 8, endRow = 378)
secnames <- colnames(soc2sic)[-c(1,ncol(soc2sic))]
sectornames <- unname(sapply(secnames,function(x)paste0(strsplit(x,'\\.\\.')[[1]][-1],collapse=' ')))
sectornames <- gsub('\\.',' ',sectornames)
colnames(soc2sic) <- c('soc2sic',sapply(secnames,function(x)strsplit(x,'\\.')[[1]][2]),'Total')
for(i in 2:ncol(soc2sic)){
  soc2sic[,i] <- as.numeric(soc2sic[,i])
  soc2sic[is.na(soc2sic[,i]),i] <- 0
}
setDT(soc2sic)[,soc:=substr(soc2sic,1,4)]

sec45names <- trimws(xlsx::read.xlsx(file.path(datapath,'6.economic_closures.xlsx'),sheetIndex = 1)[(1:45)+2,1])
sec45names <- unname(sapply(sec45names,function(x)strsplit(x,';')[[1]][1]))
soc2sic88 <- xlsx::read.xlsx(file.path(datapath,'soc4xsic12supfinal.xls'),sheetIndex = 3,startRow = 8, endRow = 378)
secnames88 <- colnames(soc2sic88)[-c(1,ncol(soc2sic88))]


{sectornames88 <- unname(sapply(secnames88,function(x)paste0(strsplit(x,'\\.\\.')[[1]][-1],collapse=' ')))
  sectornames88 <- gsub('\\.',' ',sectornames88)
  sectornames88 <- gsub('Manufacture of ','',sectornames88)
  joinstr <- origstr <- 1:88
  joinstr[2] <- 1
  joinstr[5] <- 4
  joinstr[7] <- 6
  joinstr[9:11] <- 9
  joinstr[12:14] <- 12
  joinstr[17] <- 16
  joinstr[30:32] <- 30
  joinstr[35:37] <- 34
  joinstr[39:40] <- 38
  joinstr[42:43] <- 41
  joinstr[50] <- 49
  joinstr[52:53] <- 51
  joinstr[56] <- 55
  joinstr[58:59] <- 57
  joinstr[62:67] <- 61
  joinstr[68:73] <- 68
  # http://siccodesupport.co.uk/sic-division.php?division=99
  joinstr[88] <- 74
  joinstr[77:78] <- 76
  joinstr[80:83] <- 79
  joinstr[84:85] <- 84
  joinstr[87] <- 86
  mapids <- joinstr
  dups <- duplicated(joinstr)
}
to_numeric <- function(x){
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}
for(i in unique(joinstr)){
  sectornames88[i] <- paste0(sectornames88[joinstr==i],collapse='+')
  cols <- which(joinstr==i) + 1
  vals <- 0
  for(cl in cols)
    vals <- vals + to_numeric(soc2sic88[,cl])
  soc2sic88[,cols[1]] <- vals
  joinstr[i] <- paste0(origstr[joinstr==i],collapse='+')
}
sectornames88 <- sectornames88[!dups]
soc2sic45 <- soc2sic88[,c(1,which(!dups)+1)]
joinstr <- joinstr[!dups]
# View(cbind(joinstr,sectornames88,sec45names))

colnames(soc2sic45) <- c('soc2sic',sec45names)
soc2sic45$Total <- apply(soc2sic45[,2:46],1,sum)
setDT(soc2sic45)[,soc:=substr(soc2sic,1,4)]
target_dist <- soc2sic45$Total[-c(nrow(soc2sic45))]/soc2sic45$Total[nrow(soc2sic45)]
soccodes <- soc2sic45$soc[-c(nrow(soc2sic45))]

## process codes ################################

get_occ_codes <- function(dt, colname){
  dt[tolower(get(colname))%in%om$major_label,
     occcode:=om$major[match(tolower(get(colname)),om$major_label)],by=.(get(colname))]
  dt[tolower(get(colname))%in%om$sub_major_label,
     occcode:=om$sub_major[match(tolower(get(colname)),om$sub_major_label)],by=.(get(colname))]
  dt[tolower(get(colname))%in%om$minor_label,
     occcode:=om$minor[match(tolower(get(colname)),om$minor_label)],by=.(get(colname))]
  dt[tolower(get(colname))%in%om$description,
     occcode:=om$unit[match(tolower(get(colname)),om$description)],by=.(get(colname))]
  
  dt[,occcode1:=substr(occcode,1,1),by=occcode]
  dt[,occcode2:=substr(occcode,1,2),by=occcode]
  dt <- subset(dt,!is.na(occcode))
  return(dt)
}

ukdt <- get_occ_codes(ukdt,'part_occupation')
chdt <- get_occ_codes(chdt,'ch01occr')
nldt <- get_occ_codes(nldt,'nl01occr')
sapply(list(nldt,bedt,chdt,ukdt),nrow)


## get weights for sectors #####################

dtlist <- list(nl=nldt,ch=chdt,uk=ukdt)
dtlist$all <- do.call(rbind,dtlist)
dt_for_freq <- dtlist$all[,.(part_wave_uid,occcode,occcode1,occcode2,n_cnt_work)]

## to whom (by age) ################################
# who do workers contact? proportions by age

cntocc <- cnts[dt_for_freq,,on='part_wave_uid']
cntocc <- subset(cntocc,cnt_main_type=='Work')
cntocc[,n_ind:=sum(cnt_mass=='individual'),by=part_wave_uid]
cntocc[,n_id:=sum(!is.na(contact)),by=part_wave_uid]
cntocc[,n_tot:=.N,by=part_wave_uid]
cntocc[,freq:=cnt_frequency%in%c('1-2 days', '3-7 days')]
cntocc[,n_WW:=sum(freq==1),by=part_wave_uid]
cntocc[,n_CW:=sum(freq==0&!is.na(cnt_frequency)),by=part_wave_uid]
cntocc[,n_under18:=sum(cnt_age_est_max<18),by=part_wave_uid]
cntocc[,n_workingage:=sum(cnt_age_est_max<65&cnt_age_est_min>17),by=part_wave_uid]
cntocc[,n_65plus:=sum(cnt_age_est_min>=65),by=part_wave_uid]
saveRDS(cntocc,'store/contactsoccupation.Rds')

fracage <- unique(cntocc[,.(part_wave_uid,occcode,n_under18,n_workingage,n_65plus)])[,.(
  under18=sum(n_under18)/sum(n_under18+n_workingage+n_65plus),
  workingage=sum(n_workingage)/sum(n_under18+n_workingage+n_65plus),
  `65plus`=sum(n_65plus)/sum(n_under18+n_workingage+n_65plus)),by=.(occcode)]

## back to weights ###############################


sectorids <- list()
sectorids45 <- list()
for(j in 1:length(dtlist)){
  dtworkers <- dtlist[[j]]
  
  ## get weights for soc
  isco88codes <- unique(dtworkers$occcode)
  socisco88_contengency <- matrix(0,nrow=length(isco88codes),ncol=length(soccodes))
  for(i in 1:length(unique(dtworkers$occcode))){
    isco8 <- unique(dtworkers$occcode)[i]
    if (nchar(isco8) == 1){
      isco0 <- unique(subset(isco8to0,isco8==isco1)$ISCO.08)
    }else if (nchar(isco8) == 2){
      isco0 <- unique(subset(isco8to0,isco2==isco8)$ISCO.08)
    }else if (nchar(isco8) == 3){
      isco0 <- unique(subset(isco8to0,isco3==isco8)$ISCO.08)
      if(isco8==110)
        break
    }else if (nchar(isco8) == 4){
      isco0 <- unique(subset(isco8to0,ISCO.88==isco8)$ISCO.08)
    }
    isco2socsub <- setDT(subset(tisco2soc,ISCO.08%in%isco0))
    socweights <- isco2socsub[SOC.2010%in%soccodes,sum(Proportion),by=SOC.2010]
    
    iscoid <- which(isco88codes==isco8)
    socids <- match(socweights$SOC.2010,soccodes)
    socisco88_contengency[iscoid,socids] <- 1#socweights$V1/sum(socweights$V1)
  }
  occcode_contingency <- sapply(unique(dtworkers$occcode),function(x)sum(dtworkers$occcode==x))
  sumsoc <- sapply(1:ncol(socisco88_contengency),function(x)sum(socisco88_contengency[,x]*occcode_contingency))
  socscalar <- target_dist/sumsoc
  
  ## get weights for sectors
  sectorid <- matrix(0,nrow=nrow(dtworkers),ncol=21)
  sectorid45 <- matrix(0,nrow=nrow(dtworkers),ncol=45)
  for(i in 1:nrow(dtworkers)){
    isco8 <- dtworkers$occcode[i]
    iscoid <- which(isco88codes==isco8)
    socmatches <- soccodes[which(socisco88_contengency[iscoid,]>0)]
    socweights <- socscalar[which(soccodes%in%socmatches)]
    ##!! this is the uk map
    sicsocrows <- subset(soc2sic,soc%in%socmatches)
    allsectors <- as.matrix(sicsocrows[,-c(1,(ncol(soc2sic)-1):ncol(soc2sic)),with=F])
    sectorrep <- t(allsectors) %*% socweights
    sectorid[i,] <- t(sectorrep)/sum(sectorrep)
    ##!! 45 sectors
    sicsocrows <- subset(soc2sic45,soc%in%socmatches)
    allsectors <- as.matrix(sicsocrows[,-c(1,(ncol(soc2sic45)-1):ncol(soc2sic45)),with=F])
    sectorrep <- t(allsectors) %*% socweights
    sectorid45[i,] <- sectorrep/sum(sectorrep)
  }
  sectorids[[j]] <- sectorid
  sectorids45[[j]] <- sectorid45
  # print(which(is.na(apply(sectorid,1,sum))))
  # print(apply(sectorid,2,sum))
}

## make occ and sec plots ###############################################

nsamples <- 1000
get_mean2 <- function(dt_, group_var_, cnt_var = "n_cnt"){
  top <- dt_[, .(num = .N), by = .(country, get(group_var_))]
  contact_ <- dt_[, .(n_cnt = mean(get(cnt_var), na.rm = TRUE),
                      n_cnt25 = quantile(get(cnt_var), .25, na.rm = TRUE),
                      n_cnt75 = quantile(get(cnt_var), .75, na.rm = TRUE)),
                  by = .(country, get(group_var_))]
  x1 <- merge(top, contact_)
  x1[!(get %in% c("Unknown", "Other", "remove")),
     text:= paste0(formatC(n_cnt, digits = 1, format = "f"),
                   " (",
                   formatC(n_cnt25, digits = 1, format = "f"),
                   " - ",
                   formatC(n_cnt75, digits = 1, format = "f"),
                   ")")]
  x1[get %in% c("Unknown", "Other"),
     num:= formatC(num, big.mark = ",")]
  x1 <- x1[get != "remove"]
  x1
}

# to get distribution to ages

countries <- names(dtlist)
dtlist$all$country <- 'all'
for(j in 1:length(dtlist)){
  dtworkers <- dtlist[[j]][,.(n_cnt_work,occcode,occcode1,occcode2,country)]
  cn <- countries[j]
  sectorid <- sectorids[[j]]
  sectorid45 <- sectorids45[[j]]
  
  resample <- list()
  for(i in 1:nsamples){
    cdt <- copy(dtworkers)
    cdt$sector <- apply(sectorid,1,function(p)sample(1:21,1,F,prob=p))
    cdt$sector45 <- apply(sectorid45,1,function(p)sample(1:45,1,F,prob=p))
    cdt <- fracage[cdt,,on='occcode']
    resample[[i]] <- cdt
  }
  newdt <- do.call(rbind,resample)
  seccontactdist <- newdt[,lapply(.SD,mean,na.rm=T),.SDcols=c('under18','workingage','65plus'),by=sector45]
  setorder(seccontactdist,'sector45')
  write.csv(seccontactdist,paste0(datapath,'sec_contact_dist_',toupper(cn),'.csv'),row.names = F)
  
  seccontactdist$sector45 <- sec45names
  msecdist <- melt(seccontactdist,id.vars = 'sector45')
  msecdist$variable <- factor(msecdist$variable, levels=c("under18","workingage","65plus"),labels=c('Under 18','Working age','Over 64'))
  msecdistp <- ggplot(msecdist) + 
    geom_bar(aes(x=sector45,y=value,fill=variable),stat='identity',position = position_fill(reverse = TRUE)) +
    coord_flip() + theme_bw(base_size = 15) +
    labs(x='',y='',fill='') +
    theme(legend.position = 'top') +
    scale_y_continuous(expand=c(0,0),labels=c(0,.25,.5,.75,1)) + 
    scale_fill_viridis(discrete=T, name="",option='inferno',direction=1) +
    scale_x_discrete(limits = rev(levels(msecdist$sector45)))
  ggsave(msecdistp,filename = paste0('store/',cn,'sec_dist_age.png'),height=10,width=10)
  
  ylab <- paste0('Contacts, ',toupper(cn))
  
  # sector 21
  per_sec1 <- get_mean2(newdt, "sector", cnt_var = "n_cnt_work")
  setorder(per_sec1,n_cnt)
  per_sec1[,sector:=sectornames[get],by=get]
  per_sec1[,sector:=factor(sector,levels=unique(sector))]
  p0 <- ggplot(per_sec1) +
    geom_errorbar(aes(x=sector,ymin=n_cnt25,ymax=n_cnt75),colour='azure4',width=.5) + 
    geom_point(aes(x=sector,y=n_cnt),colour='blue4',shape=18,size=4) +
    labs(x='',y=ylab) +
    coord_flip() + theme_bw(base_size=15)
  ggsave(p0,filename = paste0('store/',cn,'sector21.png'),height=4.5,width=9)
  
  # sector 45
  per_sec1 <- get_mean2(newdt, "sector45", cnt_var = "n_cnt_work")
  setorder(per_sec1,n_cnt)
  per_sec1[,sector:=sec45names[get],by=get]
  per_sec1[,sector:=factor(sector,levels=unique(sector))]
  p0 <- ggplot(per_sec1) +
    geom_errorbar(aes(x=sector,ymin=n_cnt25,ymax=n_cnt75),colour='azure4',width=.5) + 
    geom_point(aes(x=sector,y=n_cnt),colour='blue4',shape=18,size=4) +
    labs(x='',y=ylab) +
    coord_flip() + theme_bw(base_size=15)
  ggsave(p0,filename = paste0('store/',cn,'sector45.png'),height=9,width=15)
  
  write.csv(per_sec1[match(sec45names,per_sec1$sector),.(sector,n_cnt)],
            file.path(datapath,'sectorcontacts.csv'),row.names = F)
  
  # occupation
  per_occ1 <- get_mean2(dtworkers, "occcode1", cnt_var = "n_cnt_work")
  per_occ2 <- get_mean2(dtworkers, "occcode2", cnt_var = "n_cnt_work")
  subsm <- subset(occmap,ISCO_version=='ISCO-88')
  per_occ1[,occupation:=subsm$major_label[match(get,subsm$major)],by=get]
  per_occ2[,occupation:=subsm$sub_major_label[match(get,subsm$sub_major)],by=get]
  setorder(per_occ1,n_cnt)
  setorder(per_occ2,n_cnt)
  per_occ1[,occupation:=factor(occupation,levels=unique(occupation))]
  per_occ2[,occupation:=factor(occupation,levels=unique(occupation))]
  
  p1 <- ggplot(per_occ1) +
    geom_errorbar(aes(x=occupation,ymin=n_cnt25,ymax=n_cnt75),colour='azure4',width=.5) + 
    geom_point(aes(x=occupation,y=n_cnt),colour='blue4',shape=18,size=4) +
    labs(x='',y=ylab) +
    coord_flip() + theme_bw(base_size=15)
  ggsave(p1,filename = paste0('store/',cn,'occupation1.png'),height=4.5,width=9)
  
  
  p2 <- ggplot(subset(per_occ2,!is.na(occupation))) +
    geom_errorbar(aes(x=occupation,ymin=n_cnt25,ymax=n_cnt75),colour='azure4',width=.5) + 
    geom_point(aes(x=occupation,y=n_cnt),colour='blue4',shape=18,size=4) +
    labs(x='',y=ylab) +
    coord_flip() + theme_bw(base_size=15)
  ggsave(p2,filename = paste0('store/',cn,'occupation2.png'),height=4.5,width=12)
}


## belgium ##############################################

be_set <- get_mean2(bedt, "be02occhi", cnt_var = "n_cnt_work")
setorder(be_set,n_cnt)
be_set[,get:=factor(get,levels=unique(get))]
bep <- ggplot(be_set) +
  geom_errorbar(aes(x=get,ymin=n_cnt25,ymax=n_cnt75),colour='azure4',width=.5) + 
  geom_point(aes(x=get,y=n_cnt),colour='blue4',shape=18,size=4) +
  labs(x='',y='Contacts') +
  coord_flip() + theme_bw(base_size=15)
ggsave(bep,filename = 'store/be.png',height=4.5,width=12)

nrow(dt[(!is.na(ch01occr)&tolower(ch01occr)%in%occuptions|
           !is.na(part_occupation)&tolower(part_occupation)%in%occuptions|
           !is.na(nl01occr)&tolower(nl01occr)%in%occuptions|
           !is.na(be02occhi)&tolower(be02occhi)%in%be_set[-c(1:3),get])&survey_round == 1000,.(country)])


