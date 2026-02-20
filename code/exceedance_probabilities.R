
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
deathsamples <- costsamples <- sampleorders <- list()
for(refsl in 1:nbscens)
  deathsamples[[refsl]] <- costsamples[[refsl]] <- sampleorders[[refsl]] <- matrix(0,nrow=sample_size,ncol=ncountries)
for(i in 1:ncountries){
  whichig <- rep(1:3,times=samplefracs)[i]
  z <- mvrnorm(sample_size,mu=rep(0, 2),Sigma=sigma,empirical=T)
  z <- z[order(z[,1]),]
  rank1 <- rank(z[,1])
  rank2 <- rank(z[rank1,2])
  newsample <- order(z[rank1,2])
  for(refsl in 1:nbscens){
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
  write.csv(unique(expvals[,.(from,to,LQ,UQ)]),paste0('../cepi_results/Delta_LIR_IQR_pc_GDP_',bau_names[bau_scen],'.csv'),row.names = F, quote = F)
  
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
  write.csv(unique(absscencosts[,.(scenario,LQ,UQ)]),paste0('../cepi_results/LIR_IQR_pc_GDP',bau_names[bau_scen],'.csv'),row.names = F, quote = F)
  # }
  
  # plot LIR vs Delta LIR as % of counterfactual
  p <- ggplot( data.frame(loss=abscosts,value=rowSums(slvaluelist[[9]][[bau_scen]]$slvalue)*average_vsl/(average_gdp*ncountries) * 100 / abscosts * 100)) +
    geom_point(aes(x=loss,y=value),colour='midnightblue') +
    theme_bw(base_size = 16) +
    labs(x='LIR, % global GDP',y=expression(Delta*"LIR, % of counterfactual LIR"))
  ggsave(p,filename=paste0('../cepi_results/valueloss.png'))
  
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
    
    lb_cost_sheet = readxl::read_xlsx(lb_file,sheet = "Cost Breakdown")
    costscens = colnames(lb_cost_sheet)[-c(1:3)]
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
          vjust = 1.4,       # “above” the segment (perpendicular offset)
          upright = TRUE,     # keeps text readable (flips if needed)
          text_smoothing = 0, # no curve smoothing (straight line)
          size=5, show.legend=F # 
        ) +
        labs(x = 'Decrease in LIR vs. BAU, % GDP', y = 'Increase in costs vs BAU, billion $') +
        theme_bw(base_size = 15))
    #ggsave(segplot,filename=paste0('../cepi_results/dominancediag_',bau_names[bau_scen],'.png'),width=7,height=7)
    
    expvalplot[,LQplus:=signif(quantile(LIRplus * timehor,c(1)/4),2),by=.(from,to)]
    expvalplot[,UQplus:=signif(quantile(LIRplus * timehor,c(3)/4),2),by=.(from,to)]
    dom = unique(expvalplot[,.(to,LQplus,UQplus)])
    setorder(dom,to)
    colourlists$bau = bau_names
    # dominance analysis
    dom[,scennumber:=which(scenario_names==to),by=to]
    dom[,scencol:=names(colourlists)[sapply(colourlists,function(x)to%in%x)],by=scennumber]
    dom$lq = alllq[match(dom$to,costscens)] - 0.2 # adjust to aid visibility
    dom$uq = alluq[match(dom$to,costscens)] + 0.2
    dom <- subset(dom,!is.na(lq))
    dom[,labelx:=UQplus + 0.17] # adjust to aid visibility
    dom[,labely:=lq]
    cscheme = "plasma"
    dom$scencol = factor(dom$scencol,
                          levels=c('bpsv','capres','ssv200','ssv100','eq'),
                         labels=c('BPSV', '365 days','200 days','100 days', 'Equality + Delivery'))
    
    (domplot = ggplot(dom) + 
        geom_hline(yintercept=0,linewidth=1.25,colour='grey') +
        geom_vline(xintercept=0,linewidth=1.25,colour='grey') +
        geom_rect(aes(xmin=LQplus,xmax=UQplus,
                      ymin=lq,#/gdp2025*100 - .0001,
                      ymax=uq,#/gdp2025*100 + .0001,
                      fill=scencol),colour='white',linewidth=0,alpha=0.5)  +
        # scale_colour_manual(values=c(`FALSE`='white',`TRUE`=NA)) + 
        # geom_textsegment(
        #   aes(x = LQplus, y = lq,#/gdp2025*100,
        #       xend = UQplus, yend = uq,#/gdp2025*100, 
        #       colour = scencol,label = to),
        #   hjust = 0.5,        # centered along the segment
        #   vjust = 0,       # “above” the segment (perpendicular offset)
        #   upright = TRUE,     # keeps text readable (flips if needed)
        #   text_smoothing = 0, # no curve smoothing (straight line)
        #   size=5, show.legend=F # 
        # ) +
        ggrepel::geom_label_repel(aes(x=labelx,y=labely,#/gdp2025*100,
                                      label=to,fill=scencol,
                                      color = after_scale(prismatic::best_contrast(fill))),
                                  force=.01,vjust = 0.5,hjust = 0.5,show.legend = F) +
        theme_bw(base_size=16) + 
        theme(legend.position = 'top') +
        # scale_fill_viridis_d(option = cscheme) +
        # scale_colour_viridis_d(option = cscheme) +
        labs(x="Decrease in pandemic loss, % GDP",
             y='Increase in preparedness cost, billion USD',fill='')) #+
        # ggtitle(paste0('Costs and impacts accumulated over ',timehor,' years, relative to ',bau_names[bau_scen])))
    ggsave(domplot,filename=paste0('../cepi_results/dominancehoriz_',bau_names[bau_scen],'.png'),width=7,height=7)
    
    
    comparisons <- lapply(list(bpsv = c(1:3), 
                               ssv = c(4:11), 
                               eq = c(12),
                               supp1 = c(3,5,8,11),
                               supp2 = c(2,4,7,10)),function(x) just_scen_names[x])
  }
}



## global saving, %, by scenario #################################
# slvaluelist is in terms of statistical lives to get back to gdp, assume one vsl is population-weighted average
bootsls <- ilevelvals <- list()
index <- c(0,cumsum(samplefracs))
for(j in 1:ncscens){
  # bootsls[[j]] <- ilevelvals[[j]] <- list()
  for(refsl in 1:nbscens){
    # :(j-1)
    # bootsls[[j]][[refsl]] <- ilevelvals[[j]][[refsl]] <- list()
    bootsls <- list()
    for(i in 1:length(slvaluelist[[j]][[refsl]])){
      thisval <- names(slvaluelist[[j]][[refsl]])[i]
      value <- rowSums(slvaluelist[[j]][[refsl]][[thisval]])*average_vsl/(average_gdp*ncountries) * 100
      bootsls[[thisval]] <- rep(value,boot)
      # bootsls[[j]][[refsl]][[thisval]] <- rep(value,boot)
    }
    saveRDS( bootsls, paste0('tmp/bootsls',j,'-',refsl,'.Rds'))
    for(i in 1:length(samplefracs)){
      ilvalue <- rowSums(slvaluelist[[j]][[refsl]]$slvalue[,(index[i]+1):index[i+1]])
      value <- rowSums(slvaluelist[[j]][[refsl]]$slvalue)
      saveRDS(rep(ilvalue/value*100,boot),paste0('tmp/ilevelvals',j,'-',refsl,'-',i,'.Rds'))
      #ilevelvals[[j]][[refsl]][[i]] <- rep(ilvalue/value*100,boot)
    }
  }
}

tens = seq(30,100,by=10)
plotlist <- list()
valuetablist <- list()
absvalues <- s11lir <- data.frame()
for(ex in 1:length(tens)){
  allvaluesil <- allvaluescost <- valuetabcounter <- valuetabgdp <- data.frame()
  for(refsl in 1:nbscens){
    exlist <- lapply(tens,function(x) allyvals[[refsl]]<1/x*uval&allyvals[[refsl]]>1/x*lval)
    oneinXindex <- exlist[[ex]]
    
    absvaluelist <- lapply(list(abscostslist[[refsl]],deathspermillist[[refsl]]),function(x)rep(x,boot)[oneinXindex])
    absvaluesummary <- sapply(absvaluelist,function(x)paste0(signif(quantile(x,c(1,3)/4),2),collapse='--'))
    absvalues <- rbind(absvalues,absvaluesummary)
    for(j in 1:ncscens){
      # j = scen_to_keep[ji]
      scenname = just_scen_names[j]
      # :(j-1)
      # refsl = scen_to_keep[refsli]
      ilevelvals1 = readRDS(paste0('tmp/ilevelvals',j,'-',refsl,'-1.Rds'))
      ilevelvals2 = readRDS(paste0('tmp/ilevelvals',j,'-',refsl,'-2.Rds'))
      ilevelvals3 = readRDS(paste0('tmp/ilevelvals',j,'-',refsl,'-3.Rds'))
      # in sl
      value1 <- ilevelvals1#[[j]][[refsl]][[1]] 
      printval1 <- quantile(value1[oneinXindex],c(1:3)/4,na.rm=T)
      value2 <- ilevelvals2#[[j]][[refsl]][[2]] 
      printval2 <- quantile(value2[oneinXindex],c(1:3)/4,na.rm=T)
      value3 <- ilevelvals3#[[j]][[refsl]][[3]] 
      printval3 <- quantile(value3[oneinXindex],c(1:3)/4,na.rm=T)
      allvaluesil <- rbind(allvaluesil,
                           c(income_levels[1],bau_names[refsl],scenname,printval1),
                           c(income_levels[2],bau_names[refsl],scenname,printval2),
                           c(income_levels[3],bau_names[refsl],scenname,printval3))
      
      bootsls = readRDS(paste0('tmp/bootsls',j,'-',refsl,'.Rds'))
      value1 <- bootsls$ylls/bootsls$slvalue*100 
      printval1 <- quantile(value1[oneinXindex],c(1:3)/4,na.rm=T)
      value2 <- bootsls$gdploss/bootsls$slvalue*100 
      printval2 <- quantile(value2[oneinXindex],c(1:3)/4,na.rm=T)
      value3 <- bootsls$education/bootsls$slvalue*100 
      printval3 <- quantile(value3[oneinXindex],c(1:3)/4,na.rm=T)
      allvaluescost <- rbind(allvaluescost,
                             c(c('YLL','GDP','Education')[1],bau_names[refsl],scenname,printval1),
                             c(c('YLL','GDP','Education')[2],bau_names[refsl],scenname,printval2),
                             c(c('YLL','GDP','Education')[3],bau_names[refsl],scenname,printval3))
      value4 <- bootsls$slvalue[oneinXindex]
      ## index 1 of absvaluelist is cost. (2 is deaths)
      printval <- paste0(signif(quantile(value4/absvaluelist[[1]]*100,c(1,3)/4),2),collapse='--')
      printvalabs <- paste0(signif(quantile(value4,c(1,3)/4),2),collapse='--')
      valuetabcounter <- rbind(valuetabcounter,c(bau_names[refsl],scenname,printval))
      valuetabgdp <- rbind(valuetabgdp,c(bau_names[refsl],scenname,printvalabs))
      newcolname <- paste0('val',refsl,'to',j)
      if(length(value4)!=sum(abscosttab[[refsl]]$return==tens[ex])) break
      abscosttab[[refsl]][[newcolname]][abscosttab[[refsl]]$return==tens[ex]] <- value4
      # in % gdp
      # print(summary(rowSums(valuelist[[j]][[refsl]])/ncountries))
      if(scenname=='S09'){
        printval <- paste0(signif(quantile(absvaluelist[[1]]-value4,c(1,3)/4),2),collapse='--')
        s11lir <- rbind(s11lir,c(bau_names[refsl],paste0('Once in ',tens[ex],' years'),printval))
      }
    }
  }
  colnames(allvaluescost) <- c('Cost','From scenario','To','lower','Value','upper')
  colnames(allvaluesil) <- c('Level','From scenario','To','lower','Value','upper')
  for(i in 4:6){
    allvaluesil[,i] <- as.numeric(allvaluesil[,i])
    allvaluescost[,i] <- as.numeric(allvaluescost[,i])
  }
  allvaluesil$Level <- factor(allvaluesil$Level,levels=income_levels)
  p2 <- ggplot(subset(allvaluesil,`From scenario`%in%bau_names[1]&To%in%scenario_names[scen_to_keep])) +
    geom_bar(aes(x=To,fill=Level,y=as.numeric(Value)),stat='identity',position = 'dodge',show.legend = ex==1) +
    # geom_boxplot(aes(x=To,fill=Level,lower=lower,middle=upper,ymin=lower,ymax=upper,upper=upper),
    #            stat='identity',position = 'dodge',show.legend = ex==1) +
    scale_colour_viridis(discrete=T, name="",option='plasma',end=.9) +
    geom_hline(data=data.frame(y=100*popfrac[1],Level=factor('LLMIC',levels=income_levels)),aes(yintercept=y,colour=Level),show.legend = F,linewidth=1.2) +
    geom_hline(yintercept=0,colour='black',linewidth=1) +
    scale_fill_viridis(discrete=T, name="",option='plasma',end=.9) +
    # facet_grid(`From scenario`~.,scales='free') +
    theme_bw(base_size = 15) + 
    labs(x='',y='')  +
    theme(legend.position = 'top') + guides(colour = 'none')
  
  p1 <- ggplot(subset(allvaluescost,`From scenario`%in%bau_names[1]&To%in%scenario_names[scen_to_keep])) +
    geom_hline(yintercept=0,colour='black',linewidth=1) +
    geom_bar(aes(x=To,fill=Cost,y=as.numeric(Value)),stat='identity',position = 'dodge',show.legend = ex==1) +
    # geom_boxplot(aes(x=To,fill=Cost,lower=lower,middle=upper,ymin=lower,ymax=upper,upper=upper),
    #              stat='identity',position = 'dodge',show.legend = ex==1) +
    scale_fill_viridis(discrete=T, name="",option='viridis',end=.9) +
    # facet_grid(`From scenario`~.,scales='free') +
    theme_bw(base_size = 15) + 
    labs(x='',y='Median % contribution') +
    theme(legend.position = 'top')
  
  plotlist[[ex]] <- p1 + p2
  
  colnames(valuetabcounter) <- c('From scenario','To','Value')
  colnames(valuetabgdp) <- c('From scenario','To','Value')
  valuetabcounter$counter = 'counter'
  valuetabgdp$counter = 'gdp'
  valuetablist[[ex]] <- do.call(rbind,list(valuetabcounter,valuetabgdp))
  
}
colnames(s11lir) <- c('Counterfactual','Probability','LIR')
colnames(absvalues) <- colnames(expvalues) <- c('LIR, % global GDP','Deaths per thousand')
absvalues <- rbind(absvalues,expvalues)
decades = c('thirty','forty','fifty','sixty','seventy','eighty','ninety','one hundred')
rownames(absvalues) <- paste0(c(rep(paste0('Once in ',decades,' years'),each=nbscens),rep('Expectation',nbscens)),', ',bau_names)
names(valuetablist) <- names(exlist)
plotlist[[2]]
saveRDS(plotlist[[1]] +  p1 + p2 + plot_layout(ncol = 2),'results/exXbars.Rds')
saveRDS(list(valuetablist,absvalues),'results/exXvalues.Rds')
write.csv(absvalues,'../cepi_results/counterfactual.csv')

ggsave(plotlist[[1]] +  p1 + p2 + plot_layout(ncol = 2), filename='../cepi_results/ex30bar.png',height=6.5,width=9)

## save values as % gdp and % counterfactuals ###################################################

redo <- do.call(rbind,lapply(1:length(valuetablist),function(x)
  cbind(subset(valuetablist[[x]],counter=='gdp'),x)))
redopcgdp <- subset(redo,Value!='')
redopcgdp$counter = NULL
redopcgdp$lower <- sapply(redopcgdp$Value,function(y) as.numeric(strsplit(y,'--')[[1]][1]))
redopcgdp$upper <- sapply(redopcgdp$Value,function(y) as.numeric(strsplit(y,'--')[[1]][2]))

panlabs <- seq(100,30,-10) # rev(paste0('Once in ',decades,' years')) # c('Once in thirty years','Once in one hundred years'))
(plotvalues <-
    ggplot(subset(redopcgdp,`From scenario`%in%'BAU'&To%in%just_scen_names[scen_to_keep])) +
    geom_boxplot(aes(x=To, #,
                     colour=factor(x,levels=length(valuetablist):1,labels=panlabs),
                     fill=factor(x,levels=length(valuetablist):1,labels=panlabs),
                     lower=lower,middle=upper,ymin=lower,ymax=upper,upper=upper),
                 stat='identity',position = 'dodge',show.legend = T,key_glyph=draw_key_rect) +
    scale_colour_viridis(discrete=T, name="",option='plasma',end=.9,direction=-1) +
    scale_fill_viridis(discrete=T, name="Return time",option='plasma',end=.9,direction=-1) +
    # facet_grid(~factor(`From scenario`),scales='free') +
    theme_bw(base_size = 15) + 
    labs(x='',y=expression(Delta*"LIR, % of GDP"), parse=TRUE)  +
    theme(panel.grid.major.y = element_blank() ,
          panel.grid.minor.y = element_blank() ,
          axis.ticks.y = element_blank() ,
          strip.background = element_blank(),
          legend.position = 'bottom',
          legend.spacing.x = unit(-10,'cm'),  
          legend.margin = margin(l = -.15, unit = "npc")) +
    guides(fill=guide_legend(reverse = TRUE,nrow=1,label.position='bottom'#,#keywidth=2,byrow=T
    ), colour='none' ) +
    scale_y_continuous(sec.axis=dup_axis(name='',labels='',breaks=NULL)) +
    coord_flip() )
ggsave(plotvalues,filename='../cepi_results/Delta_LIR_IQR_pc_GDP_return.png',width=4,height=8)



redopcgdp$Value <- NULL
redopcgdp$x <- tens[redopcgdp$x]
colnames(redopcgdp) <- c('From','To','Return time','LQ','UQ')
write.csv(redopcgdp,'../cepi_results/Delta_LIR_IQR_pc_GDP_return.csv',row.names = F, quote = F)



redo <- do.call(rbind,lapply(1:length(valuetablist),function(x)
  cbind(subset(valuetablist[[x]],counter=='counter'),x)))
redopcc <- subset(redo,Value!='')
redopcc$counter <- NULL
redopcc$lower <- sapply(redopcc$Value,function(y) as.numeric(strsplit(y,'--')[[1]][1]))
redopcc$upper <- sapply(redopcc$Value,function(y) as.numeric(strsplit(y,'--')[[1]][2]))
redopcc$Value <- NULL
redopcc$x <- tens[redopcc$x]
colnames(redopcc) <- c('From','To','Return time','LQ','UQ')
write.csv(redopcc,'../cepi_results/Delta_LIR_IQR_pc_counterfactual_return.csv',row.names = F, quote = F)


## reread and write to one file ###################################

sheet1 = read.csv('../cepi_results/counterfactual.csv',check.names = F)
sheet2 = read.csv('../cepi_results/Delta_LIR_IQR_pc_GDP_BAU.csv',check.names = F)
# sheet2.5 = read.csv('../cepi_results/Delta_LIR_IQR_pc_GDP_BAU2.csv',check.names = F)
sheet3 = read.csv('../cepi_results/Delta_LIR_IQR_pc_GDP_return.csv',check.names = F)
sheet4 = read.csv('../cepi_results/Delta_LIR_IQR_pc_counterfactual_return.csv',check.names = F)
sheet5 = read.csv('../cepi_results/LIR_IQR_pc_GDPBAU.csv',check.names = F)
# sheet5.5 = read.csv('../cepi_results/LIR_IQR_pc_GDPBAU2.csv',check.names = F)
# print(sheet1)
# print(sheet2)

# xlsx::write.xlsx(sheet1,file = '../cepi_results/cepi_results.xlsx',sheetName='Counterfactual (Table S5)', append=F,row.names = F)
# xlsx::write.xlsx(rbind(sheet2,sheet2.5),file = '../cepi_results/cepi_results.xlsx',sheetName='Delta LIR, % GDP (Table S6)', append=T,row.names = F)
# xlsx::write.xlsx(sheet3,file = '../cepi_results/cepi_results.xlsx',sheetName='given return, SARS-X (Table S7)', append=T,row.names = F)
# xlsx::write.xlsx(sheet4,file = '../cepi_results/cepi_results.xlsx',sheetName='as % counterfactual (Table S8)', append=T,row.names = F)
# xlsx::write.xlsx(sheet5,file = '../cepi_results/cepi_results.xlsx',sheetName='LIR, % GDP (BAU1)', append=F,row.names = F)
# xlsx::write.xlsx(sheet5.5,file = '../cepi_results/cepi_results.xlsx',sheetName='LIR, % GDP (BAU2)', append=T,row.names = F)





# frac yll of costs and values
orderedyll <- lapply(income_levels,function(x){
  best <- subset(bauresults[[1]],igroup==x&mincost==1)
  setorder(best,Costpc)
  best$YLL/(best$vsl * best$gdp)
})
orderedchoices <- lapply(income_levels,function(x){
  best <- subset(bauresults[[1]],igroup==x&mincost==1)
  setorder(best,Costpc)
  best$policy
})

# yll frac as cost increases
j <- 9; refsl <- 1
scatter <- data.frame(policy=apply(sapply(1:ncol(sampleorders[[1]]),function(x){
  ig <- rep(1:3,times=samplefracs)[x]
  orderedchoices[[ig]][sampleorders[[1]][,x]]
}),1,function(y)sum(y=='No Closures'))/sum(samplefracs)*100,
deaths=rowSums(deathsamples[[1]])/(50*1e6*ncountries)*1e3,
pcyll = rowSums(slvaluelist[[j]][[refsl]]$ylls)/rowSums(slvaluelist[[j]][[refsl]]$slvalue)*100,
weightedyllpc = apply(
  sapply(1:ncol(sampleorders[[1]]),function(x){
    ig <- rep(1:3,times=samplefracs)[x]
    orderedyll[[ig]][sampleorders[[1]][,x]]
  }),1,sum)/rowSums(costsamples[[1]])*100
)
cyl_labels <- c("weightedyllpc" = "LIR", "pcyll" = "Delta*'LIR'")
(pcyllplot <- ggplot(reshape2::melt(scatter,id.vars=c('deaths','policy'))) +
  geom_vline(xintercept = quantile((scatter$deaths),c(1,3)/4),linewidth=1.5,colour='grey') +
  geom_point(aes(x=(deaths),y=value,colour=policy)) +
    facet_grid(~ variable, labeller = as_labeller(cyl_labels, default = label_parsed)) +
  # facet_grid(~factor(variable,levels=c('weightedyllpc','pcyll'),labels=c('LIR',deparse(bquote("Delta"))))) +
  scale_x_continuous(transform='log',breaks=c(.1,1,10)) +
  theme_bw(base_size=18) +
  labs(x='Deaths per thousand',y='YLL percent',colour='% No Closures'))
# ggsave(pcyllplot,filename='results/pcYLL.png',width=10,height=6)



