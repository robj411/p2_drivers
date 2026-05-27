
## start / load #####################################

if(!file.exists('process_results.R')){
  setwd(getSrcDirectory(function(){})[1])
}

source('process_results.R')

topresults <- readRDS('results/topresults.Rds')
reslist <- readRDS('results/reslist.Rds')
alldiffs <- readRDS('results/alldiffs.Rds')
vsl_gdp <- readRDS('results/vsl_gdp.Rds')
vsls <- vsl_gdp$vsls
gdps <- vsl_gdp$gdps

scen_to_keep <- c(1:11)#[-2]

## global values ###############################


# rank correlations between countries based on chosen strategies, across countries, which are not correlated
corstr <- sapply(1:nScen,function(sl)
  sapply(income_levels,function(cn) 
    sapply(income_levels[!income_levels%in%cn] ,function(cn2) 
      cor(subset(topresults[[sl]],igroup==cn)$Costpc, subset(topresults[[sl]],igroup==cn2)$Costpc,method='spearman')
    )))

print(median(corstr))
rankcorrelation <- sqrt(median(corstr))
sample_size <- dim(topresults[[1]])[1]/length(income_levels)

popsizes <- c(718255072+3398187527, 2503136362, 1240629858)
average_vsl <- sum(sapply(vsls,mean)*popsizes)/sum(popsizes)*1e6
average_gdp <- sum(sapply(gdps,mean)*popsizes)/sum(popsizes)
popfrac <- popsizes/sum(popsizes)

sigma <- matrix(c(1, rankcorrelation,rankcorrelation, 1), nrow=2)

## begin comparisons #####################

bauresults <- lapply(bau_scens,function(x)reslist[[x]])

# vectors of numbers of deaths
deathvecs <- lapply(1:nbscens,function(bs) 
  lapply(income_levels,function(x){
    best <- subset(bauresults[[bs]],igroup==x&mincost==1)
    setorder(best,Costpc)
    best$Deaths
  }))

gdpvecs <- lapply(1:nbscens,function(bs) 
  lapply(income_levels,function(x){
    best <- subset(bauresults[[bs]],igroup==x&mincost==1)
    setorder(best,Costpc)
    best$gdplossusd
  }))

# vectors of costs (pc gdp)
# costvecs <- lapply(income_levels,function(x){
#   best <- subset(bauresults,igroup==x&mincost==1)
#   setorder(best,Costpc)
#   best$Costpc
# })

# vectors of costs (pc sl)
costslvecs <- lapply(1:nbscens,function(bs)
  lapply(income_levels,function(x){
    best <- subset(bauresults[[bs]],igroup==x&mincost==1)
    setorder(best,Costpc)
    best$Costsl
  }))

# vectors of costs (pc sl)
# orderedsamples <- lapply(income_levels,function(x){
#   best <- subset(bauresults,igroup==x&mincost==1)
#   setorder(best,Costpc)
#   best$samplei
# })

# sampling values
ncountries <- 200
samplefracs <- round(ncountries*popfrac)
ncountries <- sum(samplefracs)

# values of scenarios
slvaluelist <- ordered_outcomes <- list()
outcome_names <- c('slvalue','pcvalue','ylls','education','gdploss','scen_costsl')
for(j in 1:ncscens){
  slvaluelist[[j]] <- list()
  ordered_outcomes[[j]] <- list()
  for(refsl in 1:nbscens){
    # :(j-1)
    slvaluelist[[j]][[refsl]] <- list()
    for(i in 1:length(outcome_names)) slvaluelist[[j]][[refsl]][[i]] <- matrix(0,nrow=sample_size,ncol=ncountries)
    names(slvaluelist[[j]][[refsl]]) <- outcome_names
    ordered_outcomes[[j]][[refsl]] <- lapply(income_levels,function(x){
      best <- subset(alldiffs[[j]][[refsl]],igroup==x)
      setorder(best,Costorder)
      list(Costpc=best$Costpc,
           Costsl=best$Costsl,
           YLL=best$YLL/(best$vsl * best$gdp),
           School=best$School/(best$vsl * best$gdp),
           GDP_loss=best$GDP_loss/(best$vsl * best$gdp),
           scen_costsl=best$scenCostsl
      )
    })
  }
}

# resample for scenarios
deathsamples <- costsamples <- gdpsamples <- sampleorders <- list()
for(refsl in 1:nbscens)
  deathsamples[[refsl]] <- costsamples[[refsl]] <- gdpsamples[[refsl]] <- sampleorders[[refsl]] <- matrix(0,nrow=sample_size,ncol=ncountries)
for(i in 1:ncountries){
  whichig <- rep(1:3,times=samplefracs)[i]
  z <- mvrnorm(sample_size,mu=rep(0, 2),Sigma=sigma,empirical=T)
  z <- z[order(z[,1]),]
  rank1 <- rank(z[,1])
  rank2 <- rank(z[rank1,2])
  newsample <- order(z[rank1,2])
  for(refsl in 1:nbscens){
    gdpsamples[[refsl]][,i] <- gdpvecs[[refsl]][[whichig]][newsample]
    costsamples[[refsl]][,i] <- costslvecs[[refsl]][[whichig]][newsample]
    deathsamples[[refsl]][,i] <- deathvecs[[refsl]][[whichig]][newsample]
    sampleorders[[refsl]][,i] <- newsample
    for(j in 1:ncscens){
      # :(j-1)
      slvaluelist[[j]][[refsl]]$slvalue[,i] <- ordered_outcomes[[j]][[refsl]][[whichig]]$Costsl[newsample]
      slvaluelist[[j]][[refsl]]$ylls[,i] <- ordered_outcomes[[j]][[refsl]][[whichig]]$YLL[newsample]
      slvaluelist[[j]][[refsl]]$education[,i] <- ordered_outcomes[[j]][[refsl]][[whichig]]$School[newsample]
      slvaluelist[[j]][[refsl]]$gdploss[,i] <- ordered_outcomes[[j]][[refsl]][[whichig]]$GDP_loss[newsample]
      slvaluelist[[j]][[refsl]]$scen_costsl[,i] <- ordered_outcomes[[j]][[refsl]][[whichig]]$scen_costsl[newsample]
    }
  }
}


## get marani epidemic / pandemic dataset #########################

xlsxdata <- readxl::read_xlsx('../../mevd/Epidemics dataset 21 March 2021.xlsx',n_max = 541)[,1:8]
colnames(xlsxdata)[2:5] <- c('yearstart','yearend','deaths','pop')
# fix bug
xlsxdata$deaths[xlsxdata$deaths==1800] <- 1.8
setDT(xlsxdata)
# take same subset as marani
xlsxdata <- subset(xlsxdata,yearstart>1600)
# what's the fraction of deaths since 1900 that come from sars
fracs <- xlsxdata[deaths>=0&!Disease%in%'Pneumonia'&yearstart>1900,sum(deaths/pop),by=Disease%in%c('SARS')]$V1
sarsprob <- fracs[2]/sum(fracs)
n1900 <- nrow(xlsxdata[deaths>0&!Disease%in%'Pneumonia'&yearstart>1900,])


allyears = seq(min(xlsxdata$yearstart),max(xlsxdata$yearend))
for(i in 1:nrow(xlsxdata))
  allyears = allyears[!allyears %in% seq(xlsxdata$yearstart[i],xlsxdata$yearend[i])]
c(min(xlsxdata$yearstart),max(xlsxdata$yearend),nrow(xlsxdata))


# define functions to use gpd
library(evir)
use_evir = T
get_gpd_params <- function(smps, threshold){
  if(use_evir){
    # fp <- gpd(smps, threshold=threshold, method = "pwm")
    # fp <- gpd(smps, threshold=threshold, method = "ml")
    tryCatch(
      #try to do this
      {
        fp <- gpd(smps, threshold=threshold, method = "ml")
        sigmau <- fp$par.ests[['beta']] # 
        xi <-  fp$par.ests[['xi']] # 
      },
      #if an error occurs, tell me the error
      error=function(e) {
        message('An Error Occurred')
        print('using pwm')
        fp <- gpd(smps, threshold=threshold, method = "pwm")
        sigmau <- fp$par.ests[['beta']] # 
        xi <-  fp$par.ests[['xi']] # 
      },
      #if a warning occurs, tell me the warning
      warning=function(w) {
        # message('A Warning Occurred')
        # print(w)
        # return(NA)
      }
    )
  }else{
    fp <- fgpd(smps-threshold)
    sigmau <- fp$sigmau # 
    xi <-  fp$xi #
  }
  c(xi,sigmau)
}

get_gpd_probs = function(smps, threshold, gpdparams){
  if(use_evir){
    evir::pgpd(smps, mu=threshold, xi=gpdparams[1], beta=gpdparams[2])
  }else{
    evmix::pgpd(smps-threshold, xi=gpdparams[1], sigmau=gpdparams[2])
  }
}

# plot marani data and basic gpd fit
mu <- 0.001
prob_bulk <- sum(xlsxdata[,deaths/pop*1000<mu])/nrow(xlsxdata)
data <- subset(xlsxdata,deaths/pop*1000>mu)
severity <- data$deaths/data$pop*1000
gpdparams = get_gpd_params(severity, mu)
exceedance_probs_marani <- get_gpd_probs(severity, mu, gpdparams)
p <- ggplot() + geom_point(aes(x=severity,y=1-(prob_bulk+(1-prob_bulk)*exceedance_probs_marani)),size=5,colour='grey')  +
  scale_x_log10(labels=label_log(digits=1)) +
  scale_y_log10(labels=label_log(digits=1)) +
  theme_bw(base_size=16) +
  labs(x='Deaths per thousand population',y='Exceedance probability')

deathspermil <- lapply(deathsamples,function(x) sort(rowSums(x)/(50*1e6*ncountries)*1e3,decreasing=F))
exceedance_probs <- lapply(deathspermil,function(x) get_gpd_probs(x,mu,gpdparams))#get_gpd_params(x,mu)))
yvals <- lapply(exceedance_probs,function(x) 1-(prob_bulk+(1-prob_bulk)*x))
# deaths per thousand people
(expdeaths <- lapply(1:length(bau_scens),function(x) sum(-diff(yvals[[x]])*deathspermil[[x]][-1])+yvals[[x]][1]*deathspermil[[x]][1]))
# total deaths per year
lapply(expdeaths,function(x) x/1e3*(50*1e6*ncountries))

for(bau_scen in 1:nbscens){
  colind <- rank(yvals[[bau_scen]])%in%c(round(sample_size/4):round(3*sample_size/4))
  p1 <- p + geom_point(aes(x=deathspermil[[bau_scen]],y=yvals[[bau_scen]],colour=colind),alpha=1,size=1.5,show.legend = F) +
    scale_colour_manual(values=c(`FALSE`='chocolate3',`TRUE`='midnightblue'))
  ggsave(p1,filename=paste0('results/exceedance-',bau_names[bau_scen],'.png'))
  print(quantile(deathspermil[[bau_scen]],c(1,3)/4))
  print(quantile(yvals[[bau_scen]],c(1,3)/4))
}



## plot bootstrap marani ##################
nyears = max(xlsxdata$yearstart)-min(xlsxdata$yearstart)+1
exceedance_probs <- deathsx <- sam <- c()
nrep <- 50
mu_min = -1.35
mu_max = -0.75
set.seed(0)
for(i in 1:nrep){
  mu <- 10^runif(1,mu_min,mu_max)
  datares <- xlsxdata[sample(1:nrow(xlsxdata),nrow(xlsxdata),replace=T),]
  severity <- datares$deaths/datares$pop*1000
  sev_tail = severity[severity>mu]
  n_excess = length(sev_tail)
  prob_tail = n_excess/nyears
  prob_bulk <- 1 - prob_tail # here the "bulk" includes "no event" (because the denominator is nyears, not nevents)
  # print(i)
  fp <- get_gpd_params(sev_tail, mu) # fgpd(newsample-mu)
  deathspermilex <- 10^seq(log10(mu),2,by=.1)
  sam <- c(sam,rep(i,each=length(deathspermilex)))
  deathsx <- c(deathsx, deathspermilex)
  # exceedance_probs <- c(exceedance_probs,pgpd(deathspermilex-mu,xi=fp$xi,sigmau=fp$sigmau))
  exceedance_probs <- c(exceedance_probs,get_gpd_probs(deathspermilex, mu, fp) )
}

p <- ggplot(data.frame(y=exceedance_probs,x=deathsx,sam=sam)) + 
  geom_line(aes(x=x,y=1-(prob_bulk+(1-prob_bulk)*y),group=sam),linewidth=1,colour='grey',alpha=.5)  +
  scale_x_log10(labels=label_log(digits=1)) +
  scale_y_log10(labels=label_log(digits=1)) +
  theme_bw(base_size=16) +
  labs(x='Deaths per thousand population',y='Exceedance probability')
p
ggsave(p,filename='results/bootstrapmarani.png',width=6,height=5)


## probabilities of synthetic pandemics ##################################

## function to get expected values based on exceedance probabilities
get_ep <- function(vals,deathorder,probs){
  orderedvals <- vals[deathorder]
  expval <- sum(-diff(probs)*orderedvals[-1])+probs[1]*orderedvals[1]
  expval
}

## resample data and compute probabilities
deathspermillist <- lapply(deathsamples,function(x) rowSums(x)/(50*1e6*ncountries)*1e3)
abscostslist <- lapply(costsamples,function(x) rowSums(x)*average_vsl/(average_gdp*ncountries) * 100)
samples_for_pgpd <- with(subset(xlsxdata,Location=='World'),deaths/pop*1000)
boot <- 10000; uval <- 1.05; lval <- 0.95
allyvals <- abscosttab <- list()
expvalues <- data.frame()
for(bau_scen in 1:nbscens){
  abscosts <- abscostslist[[bau_scen]]
  deathspermil <- deathspermillist[[bau_scen]]
  deathorder <- order(deathspermil)
  allyvals[[bau_scen]] <- matrix(0,ncol=boot,nrow=length(deathspermil))
  deathests <- costests <- gpdparams <- c()
  # abscosttab <- expvals <- data.frame()
  cov_ex <- c()
  abscosttablist <- expvalslist <- scencostlist <- list()
  for(i in 1:boot){
    # print(i)
    set.seed(i)
    
    # new gpd function
    mu <- 10^runif(1,mu_min,mu_max)
    datares <- xlsxdata[sample(1:nrow(xlsxdata),nrow(xlsxdata),replace=T),]
    
    if('SARS'%in%subset(datares,deaths>mu&yearstart>1930)$Disease){
      setorder(datares,yearstart)
      testdata = copy(datares)
    }else{
      testdata = copy(xlsxdata)
    }
    mod <- nnet::multinom(Disease=='SARS' ~ yearstart, subset(testdata,deaths>mu&yearstart>1930),trace=F)
    fittab <- fitted(mod)
    # psars[i] = c(tail(fittab[,colnames(fittab)=='SARS'],1))
    sarsprobs = c(tail(fittab[,1],1))
    # psars[i] = sarsprobs
    
    severity <- datares$deaths/datares$pop*1000
    sev_tail = severity[severity>mu]
    n_excess = length(sev_tail)
    prob_tail = n_excess/nyears
    prob_bulk <- 1 - prob_tail # here the "bulk" includes "no event" (because the denominator is nyears, not nevents)
    # fit lnorm to bulk
    # fit_ln <- fitdist(severity[severity>0], "lnorm")
    # lnormdist = distr::Lnorm(fit_ln$estimate[1], fit_ln$estimate[2])
    # fit gpd to tail
    fp <- get_gpd_params(sev_tail, mu) # fgpd(newsample-mu)
    sigmau <- fp[2] # 
    xi <-  fp[1] # 
    gpdparams <- rbind(gpdparams,c(sigmau,xi,mu,prob_bulk,sarsprobs))
    cov_ex[i] <- (1-(prob_bulk+prob_tail*get_gpd_probs(samples_for_pgpd, mu, fp))) # apply to covid
    
    # apply to counterfactual (deaths per mil)
    exceedance_probs <- get_gpd_probs(deathspermil, mu, fp)*prob_tail + prob_bulk
    if(any(is.na(exceedance_probs)))
      exceedance_probs[is.na(exceedance_probs)] <- 1#distr::p(lnormdist)(deathspermil[is.na(exceedance_probs)])*prob_bulk
    yvals <- 1-exceedance_probs
    allyvals[[bau_scen]][,i] <- yvals
    
    # translate to return time
    for(return_time in seq(30,100,by=10)){
      newrow <- cbind(abscosts[yvals<1/return_time*uval&yvals>1/return_time*lval],sigmau,xi,mu,prob_bulk,return_time)
      if(length(newrow)>5) 
        abscosttablist[[length(abscosttablist)+1]] <- newrow
    }
    
    # order to get the expectation
    yvals <- yvals[deathorder]
    # sarsprobs <- rbeta(1,sarsprob*n1900,(1-sarsprob)*n1900)
    deathests[i] <- get_ep(vals=deathspermil,deathorder,probs=yvals * sarsprobs)
    costests[i] <- get_ep(abscosts,deathorder,yvals * sarsprobs)
    
    # get values and costs of scenarios under each return time
    for(j in 1:ncscens){
      vals <- rowSums(slvaluelist[[j]][[bau_scen]]$slvalue)
      expvalslist[[length(expvalslist)+1]] <- c(just_scen_names[j],bau_names[bau_scen], get_ep(vals,deathorder,yvals * sarsprobs) )
      # if(bau_scen==1){
      scencosts <- rowSums(slvaluelist[[j]][[bau_scen]]$scen_costsl)
      scencostlist[[length(scencostlist)+1]] <- c(just_scen_names[j],get_ep(scencosts,deathorder,yvals * sarsprobs) )
      # }
    }
  }
  sort(sapply(ls(),function(x)object.size(get(x))),decreasing = F)
  
  abscosttab[[bau_scen]] <- as.data.frame(do.call(rbind,abscosttablist))
  rm(abscosttablist)
  print('estimated deaths, millions')
  print(signif(summary(deathests/1000*8.2e9)/1e6,2))
  print('estimated deaths per thousand')
  print(signif(summary(deathests),2))
  print('covid return time')
  print(signif(summary(1/cov_ex),2))
  expvalues <- rbind(expvalues,
                     c(paste0(signif(quantile(costests,c(1,3)/4),2),collapse='--'), 
                       paste0(signif(quantile(deathests,c(1,3)/4),2),collapse='--')))
  colnames(gpdparams) <- c('sigmau','xi','mu','p_bulk','p_sars')
  
  # collate expected values (Delta LIR)
  expvals <- as.data.frame(do.call(rbind,expvalslist))
  rm(expvalslist)
  colnames(expvals) <- c('to','from','value')
  setDT(expvals)
  expvals[,LIR:=as.numeric(value)*average_vsl/(average_gdp*ncountries) * 100]
  expvals[,LQ:=signif(quantile(LIR,c(1)/4),2),by=.(from,to)]
  expvals[,UQ:=signif(quantile(LIR,c(3)/4),2),by=.(from,to)]
  expvaltab <- dcast(expvals[,paste0(signif(quantile(LIR,c(1,3)/4),2),collapse='--{}'),by=.(from,to)],formula=from~to,fill='')
  expvaltab <- expvals[,paste0(signif(quantile(LIR,c(1,3)/4),2),collapse='--\u200B'),by=.(from,to)]
  expvaltab$from <- NULL
  colnames(expvaltab) <- c('Scenario',bau_names[bau_scen])
  saveRDS(expvaltab,paste0('results/expvals_',bau_names[bau_scen],'.Rds'))
  write.csv(unique(expvals[,.(from,to,LQ,UQ)]),paste0('../results/Delta_LIR_IQR_pc_GDP_',bau_names[bau_scen],'.csv'),row.names = F, quote = F)
  
  # collate expected costs of each scenario
  # if(bau_scen==1){
  absscencosts  <- as.data.frame(do.call(rbind,scencostlist))
  rm(scencostlist)
  colnames(absscencosts) <- c('scenario','cost')
  setDT(absscencosts)
  absscencosts[,LIR:=as.numeric(cost)*average_vsl/(average_gdp*ncountries) * 100]
  absscencosts[,LQ:=signif(quantile(LIR,c(1)/4),2),by=.(scenario)]
  absscencosts[,UQ:=signif(quantile(LIR,c(3)/4),2),by=.(scenario)]
  scencosttab <- absscencosts[,paste0(signif(quantile(LIR,c(1,3)/4),2),collapse='--'),by=.(scenario)]
  saveRDS(scencosttab,paste0('results/expvals_scenarios.Rds'))
  write.csv(unique(absscencosts[,.(scenario,LQ,UQ)]),paste0('../results/LIR_IQR_pc_GDP',bau_names[bau_scen],'.csv'),row.names = F, quote = F)
  # }
  
  # plot LIR vs Delta LIR as % of counterfactual
  p <- ggplot( data.frame(loss=abscosts,value=rowSums(slvaluelist[[9]][[bau_scen]]$slvalue)*average_vsl/(average_gdp*ncountries) * 100 / abscosts * 100)) +
    geom_point(aes(x=loss,y=value),colour='midnightblue') +
    theme_bw(base_size = 16) +
    labs(x='LIR, % global GDP',y=expression(Delta*"LIR, % of counterfactual LIR"))
  ggsave(p,filename=paste0('../results/valueloss.png'))
  
  # dominance among BPSv scenarios
  sapply(2:3,function(x) sapply(1:(x-1),function(y)sum(rowSums(slvaluelist[[x]][[bau_scen]]$slvalue) - rowSums(slvaluelist[[y]][[bau_scen]]$slvalue) > 0))/length(rowSums(slvaluelist[[x]][[bau_scen]]$slvalue))*100)
  
  ## probability sars
  
  print(summary(as.data.frame(gpdparams)[['p_sars']]))
  ps = as.data.frame(gpdparams)[['p_sars']]
  timehor = 50 # years
  horizon_end = as.numeric(format(Sys.Date(),'%Y')) + timehor
  ps_th <- expNevents <- p_two_or_more <- c()
  for(i in 1:boot){
    peryear = sample(ps[i],timehor,replace=T)
    expNevents[i] = sum(peryear*1/30)
    ps_th[i] = 1-prod(1-1/30*peryear)
    p_two_or_more[i] = 1 - pbinom(1,size=timehor,prob = 1/30*peryear[1])
  }
  cat('probability of at least one event\n')
  print(summary(ps_th))
  cat('probability of at least two events\n')
  print(summary(p_two_or_more))
  cat('expected number of events\n')
  print(summary(expNevents))
  
  ## join to costs
  # prep LB preparedness and response costs
  expvalplot = copy(expvals)
  # expvals = expvalsave
  if(bau_scen==1){
    
    # group scenarios by colour
    colourlists <- lapply(list(bpsv = c(1), 
                               capres = c(2:3), 
                               ssv200 = c(4:6), 
                               ssv100 = c(7:9), 
                               eq = c(10)),function(x) just_scen_names[x])
    
    discountrate = runif(boot, 0.02, 0.06) # sample from 0.02 to 0.06 # 0.04 #  
    gdp2025 = 113.8e6 # million USD
    discountedth = sapply(discountrate, function(x) sum(1/(1+x)^(1:timehor-1)))
    
    # lb_cost_sheet = readxl::read_xlsx('../data/20251104 updated scenario delivery and costing.xlsx',sheet = "Cost Breakdown")
    costscens = scenario_names[scen_scens]
    alllq <- alluq <- segs <- c() 
    
    for(i in 1:ncscens){
      scenocosts = cost_diffs[[1]][,i]/1e3 # billions
      scenacosts = cost_diffs[[2]][,i]/1e3 # billions
      scenrespcosts = cost_diffs[[3]][,i] # millions
      # multiply LIR by expected number of events in timehor years; add response costs, as % GDP, multiplied by expected number of events in th years
      expvalplot$LIRplus[expvalplot$to==costscens[i]] = expvals$LIR[expvals$to==costscens[i]] - 
        scenrespcosts / gdp2025 * 100 * expNevents * discountedth / timehor
      print(summary(scenrespcosts / gdp2025 * 100 * expNevents * discountedth))
      # one-time cost plus annual costs, discounted over time horizon
      qs = quantile(scenocosts + scenacosts * discountedth, c(1,3)/4)
      alllq[i] = qs[1]
      alluq[i] = qs[2]
      # test
      x = expvalplot$LIRplus[expvalplot$to==costscens[i]] * timehor
      # PLUS the response cost, times the expected number of pandemics in the time horizon, averaging over the time it likely arrives, discounting later times
      y = scenocosts + scenacosts * discountedth #+ scenrespcosts * expNevents * discountedth / timehor
      segs = rbind(segs, cbind(get_pcs(data.frame(x=x,y=y)),scenario=i))
    }
    
    segs$to = just_scen_names[segs$scenario]
    setDT(segs)
    segs[,scennumber:=which(scenario_names==to),by=to]
    segs[,scencol:=names(colourlists)[sapply(colourlists,function(x)to%in%x)],by=scennumber]
    segs[,angle:=atan2((yend-y)/1000,xend-x) * 180 / pi]
    
    (segplot <- ggplot(segs) +
        geom_segment(
          aes(x = x, y = y/1000, xend = xend, yend = yend/1000, group = scenario, colour=scencol),
          linewidth = 1.5, show.legend=F
        ) +
        geom_textsegment(
          aes(x = x, y = y/1000, xend = xend, yend = yend/1000, colour = scencol,label = to),
          hjust = 0.5,        # centered along the segment
          vjust = 1.4,       # b