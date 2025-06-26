
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

scen_to_keep <- 1:14

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


## get epidemic / pandemic dataset #########################

mu <- 0.004
xlsxdata <- readxl::read_xlsx('../../mevd/Epidemics dataset 21 March 2021.xlsx',n_max = 541)[,1:8]
colnames(xlsxdata)[2:5] <- c('yearstart','yearend','deaths','pop')
xlsxdata$deaths[xlsxdata$deaths==1800] <- 1.8
setDT(xlsxdata)
# xlsxdata <- subset(xlsxdata,Disease%in%c("Influenza","Whooping Cough","Tuberculosis","MERS","Pneumonia","SARS"))
xlsxdata <- subset(xlsxdata,yearstart>1600)
fracs <- xlsxdata[deaths>=0&!Disease%in%'Pneumonia'&yearstart>1900,sum(deaths/pop),by=Disease%in%c('SARS')]$V1
sarsprob <- fracs[2]/sum(fracs)
n1900 <- nrow(xlsxdata[deaths>0&!Disease%in%'Pneumonia'&yearstart>1900,])
print(c(sarsprob*100,n1900))
pb <- sum(xlsxdata[,deaths/pop*1000<mu])/nrow(xlsxdata)
data <- subset(xlsxdata,deaths/pop*1000>mu)


# csvdata <- read.csv('../../mevd/Epidemic16002020March2021.csv',header=F)
# colnames(csvdata)[1:4] <- c('yearstart','yearend','deaths','pop')
# csvdata$deaths[csvdata$deaths==1800] <- 1.8
# setDT(csvdata)
# csvdata[,duration:=yearend-yearstart+1]
# csvdata <- subset(csvdata,yearstart>1600&duration<10)

# pb <- sum(csvdata[,deaths/pop*1000<mu])/nrow(csvdata)
# data <- subset(csvdata,deaths/pop*1000>mu&yearstart>1600)

intensity <- data$deaths/data$pop*1000
fp <- fgpd(intensity-mu)
sigmau <- fp$sigmau # 0.0113 # 
xi <-  fp$xi # 1.41 # 
exceedance_probs_marani <- pgpd(intensity-mu,xi=xi,sigmau=sigmau)

p <- ggplot() + geom_point(aes(x=intensity,y=1-(pb+(1-pb)*exceedance_probs_marani)),size=5,colour='grey')  +
  scale_x_log10(labels=label_log(digits=1)) +
  scale_y_log10(labels=label_log(digits=1)) +
  theme_bw(base_size=16) +
  labs(x='Deaths per thousand population',y='Exceedance probability')

deathspermil <- lapply(deathsamples,function(x) sort(rowSums(x)/(50*1e6*ncountries)*1e3,decreasing=F))
exceedance_probs <- lapply(deathspermil,function(x) pgpd(x-mu,xi=xi,sigmau=sigmau))
yvals <- lapply(exceedance_probs,function(x) 1-(pb+(1-pb)*x))
# deaths per thousand people
(expdeaths <- lapply(1:length(bau_scens),function(x) sum(-diff(yvals[[x]])*deathspermil[[x]][-1])+yvals[[x]][1]*deathspermil[[x]][1]))
# total deaths per year
lapply(expdeaths,function(x) x/1e3*(50*1e6*ncountries))

for(bau_scen in 1:nbscens){
  colind <- rank(yvals[[bau_scen]])%in%c(round(sample_size/4):round(3*sample_size/4))
  p <- p + geom_point(aes(x=deathspermil[[bau_scen]],y=yvals[[bau_scen]],colour=colind),alpha=1,size=1.5,show.legend = F) +
    scale_colour_manual(values=c(`FALSE`='chocolate3',`TRUE`='midnightblue'))
  ggsave(p,filename=paste0('results/exceedance-',bau_names[bau_scen],'.png'))
  print(quantile(deathspermil[[bau_scen]],c(1,3)/4))
  print(quantile(yvals[[bau_scen]],c(1,3)/4))

}

1/(1-(pb+(1-pb)*pgpd(8.6-mu,xi=xi,sigmau=sigmau)))


## plot bootstrap marani ##################

exceedance_probs <- deathsx <- sam <- c()
nrep <- 50
for(i in 1:nrep){
  mu <- 10^runif(1,-3,-2)
  datares <- xlsxdata[sample(1:nrow(xlsxdata),nrow(xlsxdata),replace=T),]
  pb <- sum(datares[,deaths/pop*1000<mu])/nrow(datares)
  data <- subset(datares,deaths/pop*1000>mu)
  intensity <- data$deaths/data$pop*1000
  newsample <- sample(intensity,nrow(data),replace=T)
  # newsample[newsample==intensity[length(intensity)]] <- samplecovid(sum(newsample==intensity[length(intensity)]))
  
  fp <- fgpd(newsample-mu)
  deathspermilex <- 10^seq(log10(mu),2,by=.1)
  sam <- c(sam,rep(i,each=length(deathspermilex)))
  deathsx <- c(deathsx, deathspermilex)
  exceedance_probs <- c(exceedance_probs,pgpd(deathspermilex-mu,xi=fp$xi,sigmau=fp$sigmau))
}

p <- ggplot(data.frame(y=exceedance_probs,x=deathsx,sam=sam)) + 
  geom_line(aes(x=x,y=1-(pb+(1-pb)*y),group=sam),linewidth=1,colour='grey',alpha=.5)  +
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
boot <- 1000; uval <- 1.05; lval <- 0.95
allyvals <- abscosttab <- list()
expvalues <- data.frame()
for(bau_scen in 1:nbscens){
  abscosts <- abscostslist[[bau_scen]]
  deathspermil <- deathspermillist[[bau_scen]]
  deathorder <- order(deathspermil)
  allyvals[[bau_scen]] <- matrix(0,ncol=boot,nrow=length(deathspermil))
  deathests <- costests <- gpdparams <- c()
  # abscosttab <- expvals <- data.frame()
  set.seed(0)
  cov_ex <- c()
  abscosttablist <- expvalslist <- scencostlist <- list()
  for(i in 1:boot){
    print(i)
    
    # new gpd function
    mu <- 10^runif(1,-3,-2)
    datares <- xlsxdata[sample(1:nrow(xlsxdata),nrow(xlsxdata),replace=T),]
    pb <- sum(datares[,deaths/pop*1000<mu])/nrow(datares)
    data <- subset(datares,deaths/pop*1000>mu)
    intensity <- data$deaths/data$pop*1000
    newsample <- sample(intensity,nrow(data),replace=T)
    # newsample[newsample==intensity[length(intensity)]] <- samplecovid(sum(newsample==intensity[length(intensity)]))
    fp <- fgpd(newsample-mu)
    sigmau <- fp$sigmau # 0.0113 # 
    xi <-  fp$xi # 1.41 # 
    gpdparams <- rbind(gpdparams,c(sigmau,xi,mu,pb))
    cov_ex[i] <- (1-(pb+(1-pb)*pgpd(samples_for_pgpd-mu,xi=xi,sigmau=sigmau))) # apply to covid
    
    # apply to counterfactual (deaths per mil)
    exceedance_probs <- pgpd(deathspermil-mu,xi=xi,sigmau=sigmau)
    yvals <- (1-(pb+(1-pb)*exceedance_probs))
    allyvals[[bau_scen]][,i] <- yvals
    
    # translate to return time
    for(return_time in seq(30,100,by=10)){
      newrow <- cbind(abscosts[yvals<1/return_time*uval&yvals>1/return_time*lval],sigmau,xi,mu,pb,return_time)
      if(length(newrow)>5) #newrow <- cbind(NA,sigmau,xi,mu,pb,return_time){
        abscosttablist[[length(abscosttablist)+1]] <- newrow
        # abscosttab <- rbind(abscosttab,newrow)
    }
    
    # order to get the expectation
    yvals <- yvals[deathorder]
    sarsprobs <- rbeta(1,sarsprob*n1900,(1-sarsprob)*n1900)
    deathests[i] <- get_ep(deathspermil,deathorder,yvals * sarsprobs)
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
  # estimated deaths, millions
  summary(deathests/1000*8.1e9)/1e6
  # estimated deaths per thousand
  summary(deathspermil)
  # covid return time
  summary(1/cov_ex)
  expvalues <- rbind(expvalues,
                     c(paste0(signif(quantile(costests,c(1,3)/4),2),collapse='--'), 
                 paste0(signif(quantile(deathests,c(1,3)/4),2),collapse='--')))
  colnames(gpdparams) <- c('sigmau','xi','mu','pb')
  
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
  p <- ggplot( data.frame(loss=abscosts,value=rowSums(slvaluelist[[5]][[bau_scen]]$slvalue)*average_vsl/(average_gdp*ncountries) * 100 / abscosts * 100)) +
             geom_point(aes(x=loss,y=value),colour='midnightblue') +
             theme_bw(base_size = 16) +
             labs(x='LIR, % global GDP',y=expression(Delta*"LIR, % of LIR"))
  ggsave(p,filename=paste0('results/valueloss__',bau_names[bau_scen],'.png'))
  
  # prep LB preparedness and response costs
  lb_cost_sheet = readxl::read_xlsx(lb_file,sheet = 1+bau_scen)
  costscens = colnames(lb_cost_sheet)[-c(1:3)]
  firstdiff = c(which(lb_cost_sheet[,1]=="Difference to BAU"))
  costindex25 = c(which(lb_cost_sheet[firstdiff:nrow(lb_cost_sheet),2]=="0.25")) + firstdiff - 1
  costindex75 = c(which(lb_cost_sheet[firstdiff:nrow(lb_cost_sheet),2]=="0.75")) + firstdiff - 1
  lqrow = lb_cost_sheet[costindex25,]
  uqrow = lb_cost_sheet[costindex75,]
  alllq = as.numeric(lqrow[1,-c(1:3)])
  alluq = as.numeric(uqrow[1,-c(1:3)])
  dom = subset(unique(expvals[,.(from,to,LQ,UQ)]),from==bau_names[bau_scen])
  setorder(dom,to)
  # present value, up to the year 2040
  
  # group scenarios by colour
  colourlists <- lapply(list(bpsv = c(1:3), 
                             capres = c(4:5), 
                             ssv200 = c(6:8), 
                             ssv100 = c(9:11), 
                             eq = c(12)),function(x) just_scen_names[x])
  colourlists$bau = bau_names
  # dominance analysis
  dom[,scennumber:=which(scenario_names==to),by=to]
  dom[,scencol:=names(colourlists)[sapply(colourlists,function(x)to%in%x)],by=scennumber]
  dom$lq = alllq[match(dom$to,costscens)]
  dom$uq = alluq[match(dom$to,costscens)]
  dom <- subset(dom,!is.na(lq))
  dom[,labelx:=UQ]
  dom[,labely:=lq]
  # dom[scencol%in%c('bpsv','rnd','eq'),labely:=uq-2500]
  # dom[scencol%in%c('bpsv','rnd','eq'),labelx:=labelx-.0025]
  # dom[to%in%c('S18','S13','S01','S08','S07','S10'),labelx:=LQ]
  cscheme = "plasma"
  (domplot = ggplot(dom) + 
      geom_hline(yintercept=0,linewidth=1.25,colour='grey') +
      geom_vline(xintercept=0,linewidth=1.25,colour='grey') +
      geom_rect(aes(xmin=LQ*15,xmax=UQ*15,ymin=lq/1000,ymax=uq/1000,fill=scencol),colour='white',
                linewidth=1,alpha=0.5,show.legend = F) +
      ggrepel::geom_label_repel(aes(x=labelx*15,y=labely/1000,label=to,fill=scencol,
                     color = after_scale(prismatic::best_contrast(fill))),
                     force=.001,vjust = 0,hjust = 0,show.legend = F) +
      theme_bw(base_size=16) + 
      scale_fill_viridis_d(option = cscheme) +
      scale_colour_viridis_d(option = cscheme) +
      labs(x=paste0("Decrease in LIR relative to ",bau_names[bau_scen],", % of GDP"),
                    y=paste0('Total cost relative to ',bau_names[bau_scen],', billion $'),colour='') +
      ggtitle('Present value, up to the year 2040'))
  ggsave(domplot,filename=paste0('../cepi_results/dominance_',bau_names[bau_scen],'.png'),width=7,height=7)
  
  # Policy Area Relevant Scenarios Notes
  # 1. BPSV S1, S2, S3, S28, S29 S28 and S29 highlight poor cost-effectiveness; S1/S2 are near-identical for inclusion, S3 for comparison with S29
  # 2. Expedited R&D Timelines S4 (200-day), S5 (100-day), S16–S19 S16–S19 include equity and delivery variations; S5 has more benefit but at higher cost
  # 3. Capacity Reservations S8, S9 Both assume proportional distribution; S9 shows higher cost savings
  # 4. Combined R&D + Manufacturing (“Moonshot”) S4 (200-day), S9 (2B), S23 (combo), S20–S27 S27 is the “unicorn” scenario; demonstrates synergy between R&D speed and capacity savings
  # 5. Equity in Vaccine Distribution S12, S16, S17 S12 shows benefit without added cost; S16–S17 add expedited R&D (100-/200-day)
  
  comparisons <- lapply(list(bpsv = c(1:3), 
                             ssv = c(4:11), 
                             eq = c(12)),function(x) just_scen_names[x])
  domplotlist <- list()
  for(i in 1:length(comparisons)) {
    (domplotlist[[i]] = ggplot(dom) + 
       geom_hline(yintercept=0,linewidth=1.25,colour='grey') +
       geom_vline(xintercept=0,linewidth=1.25,colour='grey') +
       geom_rect(aes(xmin=LQ*15,xmax=UQ*15,ymin=lq/1000,ymax=uq/1000,fill=scencol),colour='white',linewidth=1,alpha=0.2,show.legend = F) +
       geom_rect(data=subset(dom,to%in%comparisons[[i]]),colour='white',linewidth=1,
                 aes(xmin=LQ*15,xmax=UQ*15,ymin=lq/1000,ymax=uq/1000,fill=scencol),alpha=0.7,show.legend = F) +
       # geom_errorbar(data=subset(dom,scennumber%in%comparisons[[i]]),width=0,
                     # aes(xmin=LQ*15,xmax=UQ*15,y=cost/1000,colour=scencol),linewidth=3,alpha=0.7,show.legend = F) +
       # geom_errorbar(aes(xmin=LQ*15,xmax=UQ*15,y=cost/1000,colour=scencol),linewidth=3,alpha=0.2,show.legend = F) +
       ggrepel::geom_label_repel(data=subset(dom,to%in%comparisons[[i]]),
                  force=.001,aes(x=labelx*15,y=labely/1000,label=to),vjust = 0,hjust = 0,size=4.5) +
       # geom_label(data=subset(dom,scennumber%in%comparisons[[i]]),
       #            aes(x=labelx*15,y=labely/1000,label=to,fill=scencol,
       #                color = after_scale(prismatic::best_contrast(fill))),
       #            vjust = 0,hjust = 0,show.legend = F) +
       scale_fill_viridis_d(option = cscheme) +
       scale_colour_viridis_d(option = cscheme) +
       theme_bw(base_size=16) +
       labs(x=paste0("Decrease in LIR relative to ",bau_names[bau_scen],", % of GDP"),
            y=paste0('Total cost relative to ',bau_names[bau_scen],', billion $'),colour='') +
       ggtitle('Present value, up to the year 2040'))
    print(domplotlist[[i]])
    ggsave(domplotlist[[i]],filename=paste0('../cepi_results/',names(comparisons)[i],'_',bau_names[bau_scen],'.png'),width=7,height=7)
  }
  domplot = ggarrange(domplotlist[[1]],domplotlist[[2]],domplotlist[[3]],nrow=1,labels=c('A','B','C'))
  print(domplot)
  ggsave(domplot,filename=paste0('../cepi_results/dom3_',bau_names[bau_scen],'.png'),width=16,height=5)
}

# # "final cost" voi
# colnames(abscosttab) <- c('cost',colnames(gpdparams),'return')
# abscosttab <- subset(abscosttab,!is.na(cost))
# for(i in seq(30,100,by=10)){
#   subabscosttab <- subset(abscosttab,return==i)
#   y <- log(subabscosttab$cost)
#   vary <- var(y) 
#   sourcesj <- subabscosttab[,2:5]
#   fittedvalues <- evppifit(y,sourcesj,pars=colnames(gpdparams),method='earth')
#   mi <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
#   (voi <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100)
#   print(voi)
# }
# 
# y <- deathests
# y <- costests
# vary <- var(y) 
# fittedvalues <- evppifit(y,gpdparams,pars=colnames(gpdparams),method='earth')
# mi <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
# (voi <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100)


ggplot(cbind(costests,gpdparams)) +
  geom_point(aes(x=log(xi),y=log(sigmau),colour=costests))

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
absvalues <- data.frame()
for(ex in 1:length(tens)){
  allvaluesil <- allvaluescost <- allvalues4 <- allvalues5 <- data.frame()
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
      allvalues4 <- rbind(allvalues4,c(bau_names[refsl],scenname,printval))
      allvalues5 <- rbind(allvalues5,c(bau_names[refsl],scenname,printvalabs))
      newcolname <- paste0('val',refsl,'to',j)
      if(length(value4)!=sum(abscosttab[[refsl]]$return==tens[ex])) break
      abscosttab[[refsl]][[newcolname]][abscosttab[[refsl]]$return==tens[ex]] <- value4
      # in % gdp
      # print(summary(rowSums(valuelist[[j]][[refsl]])/ncountries))
    }
  }
  colnames(allvaluescost) <- c('Cost','From scenario','To','lower','Value','upper')
  colnames(allvaluesil) <- c('Level','From scenario','To','lower','Value','upper')
  for(i in 4:6){
    allvaluesil[,i] <- as.numeric(allvaluesil[,i])
    allvaluescost[,i] <- as.numeric(allvaluescost[,i])
  }
  allvaluesil$Level <- factor(allvaluesil$Level,levels=income_levels)
  p2 <- ggplot(subset(allvaluesil,`From scenario`%in%bau_names&To%in%scenario_names[scen_to_keep])) +
    geom_bar(aes(x=To,fill=Level,y=as.numeric(Value)),stat='identity',position = 'dodge',show.legend = ex==1) +
    # geom_boxplot(aes(x=To,fill=Level,lower=lower,middle=upper,ymin=lower,ymax=upper,upper=upper),
    #            stat='identity',position = 'dodge',show.legend = ex==1) +
    scale_colour_viridis(discrete=T, name="",option='plasma',end=.9) +
    geom_hline(data=data.frame(y=100*popfrac[1],Level=factor('LLMIC',levels=income_levels)),aes(yintercept=y,colour=Level),show.legend = F,linewidth=1.2) +
    geom_hline(yintercept=0,colour='black',linewidth=1) +
    scale_fill_viridis(discrete=T, name="",option='plasma',end=.9) +
    facet_grid(`From scenario`~.,scales='free') +
    theme_bw(base_size = 15) + 
    labs(x='',y='')  +
    theme(legend.position = 'top') + guides(colour = 'none')
  
  p1 <- ggplot(subset(allvaluescost,`From scenario`%in%bau_names&To%in%scenario_names[scen_to_keep])) +
    geom_hline(yintercept=0,colour='black',linewidth=1) +
    geom_bar(aes(x=To,fill=Cost,y=as.numeric(Value)),stat='identity',position = 'dodge',show.legend = ex==1) +
    # geom_boxplot(aes(x=To,fill=Cost,lower=lower,middle=upper,ymin=lower,ymax=upper,upper=upper),
    #              stat='identity',position = 'dodge',show.legend = ex==1) +
    scale_fill_viridis(discrete=T, name="",option='viridis',end=.9) +
    facet_grid(`From scenario`~.,scales='free') +
    theme_bw(base_size = 15) + 
    labs(x='',y='Median % contribution') +
    theme(legend.position = 'top')
  
  plotlist[[ex]] <- p1 + p2
  
  colnames(allvalues4) <- c('From scenario','To','Value')
  colnames(allvalues5) <- c('From scenario','To','Value')
  valuetab <- allvalues4 # reshape2::dcast(allvalues4,formula=`From scenario`~To,value.var = 'Value',fill=''))
  valuetabgdp <- allvalues5 # reshape2::dcast(allvalues5,formula=`From scenario`~To,value.var = 'Value',fill=''))
  valuetab$counter = 'counter'
  valuetabgdp$counter = 'gdp'
  valuetablist[[ex]] <- do.call(rbind,list(valuetab,valuetabgdp))
  
}
colnames(absvalues) <- colnames(expvalues) <- c('LIR, % global GDP','Deaths per thousand')
absvalues <- rbind(absvalues,expvalues)
decades = c('thirty','forty','fifty','sixty','seventy','eighty','ninety','one hundred')
rownames(absvalues) <- paste0(c(rep(paste0('Once in ',decades,' years'),each=nbscens),rep('Expectation',nbscens)),', ',bau_names)
names(valuetablist) <- names(exlist)
plotlist[[2]]
saveRDS(plotlist[[1]] +  p1 + p2 + plot_layout(ncol = 2),'results/exXbars.Rds')
saveRDS(list(valuetablist,absvalues),'results/exXvalues.Rds')
write.csv(absvalues,'../cepi_results/counterfactual.csv')

## save values as % gdp and % counterfactuals ###################################################

redo <- do.call(rbind,lapply(1:length(valuetablist),function(x)
  cbind(subset(valuetablist[[x]],counter=='gdp'),x)))
redopcgdp <- subset(redo,Value!='')
redopcgdp$counter = NULL
redopcgdp$lower <- sapply(redopcgdp$Value,function(y) as.numeric(strsplit(y,'--')[[1]][1]))
redopcgdp$upper <- sapply(redopcgdp$Value,function(y) as.numeric(strsplit(y,'--')[[1]][2]))

panlabs <- seq(100,30,-10) # rev(paste0('Once in ',decades,' years')) # c('Once in thirty years','Once in one hundred years'))
(plotvalues <-
  ggplot(subset(redopcgdp,`From scenario`%in%bau_names&To%in%just_scen_names[scen_to_keep])) +
  geom_boxplot(aes(x=To, #,
                   colour=factor(x,levels=length(valuetablist):1,labels=panlabs),
                   fill=factor(x,levels=length(valuetablist):1,labels=panlabs),
                               lower=lower,middle=upper,ymin=lower,ymax=upper,upper=upper,
                   alpha=I(.25+x%in%3:6)),
             stat='identity',position = 'dodge',show.legend = T,key_glyph=draw_key_rect) +
  scale_colour_viridis(discrete=T, name="",option='plasma',end=.9,direction=-1) +
  scale_fill_viridis(discrete=T, name="Return time",option='plasma',end=.9,direction=-1) +
  facet_grid(~factor(`From scenario`),scales='free') +
  theme_bw(base_size = 15) + 
  labs(x='To:',y=expression(Delta*"LIR, % of GDP"), parse=TRUE)  +
  theme(panel.grid.major.y = element_blank() ,
        panel.grid.minor.y = element_blank() ,
        axis.ticks.y = element_blank() ,
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.spacing.x = unit(-10,'cm'),  
        legend.margin = margin(l = -.15, unit = "npc")) +
  guides(fill=guide_legend(reverse = TRUE,nrow=1,label.position='bottom'#,#keywidth=2,byrow=T
                           ), colour='none' ) +
  scale_y_continuous(sec.axis=dup_axis(name='From:',labels='',breaks=NULL)) +
  coord_flip() )
ggsave(plotvalues,filename='results/Delta_LIR_IQR_pc_GDP_return.png',width=4,height=6)



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
sheet2 = read.csv('../cepi_results/Delta_LIR_IQR_pc_GDP_BAU1.csv',check.names = F)
sheet2.5 = read.csv('../cepi_results/Delta_LIR_IQR_pc_GDP_BAU2.csv',check.names = F)
sheet3 = read.csv('../cepi_results/Delta_LIR_IQR_pc_GDP_return.csv',check.names = F)
sheet4 = read.csv('../cepi_results/Delta_LIR_IQR_pc_counterfactual_return.csv',check.names = F)
sheet5 = read.csv('../cepi_results/LIR_IQR_pc_GDPBAU1.csv',check.names = F)
sheet5.5 = read.csv('../cepi_results/LIR_IQR_pc_GDPBAU2.csv',check.names = F)

# xlsx::write.xlsx(sheet1,file = '../cepi_results/cepi_results.xlsx',sheetName='Counterfactual (Table S5)', append=F,row.names = F)
# xlsx::write.xlsx(rbind(sheet2,sheet2.5),file = '../cepi_results/cepi_results.xlsx',sheetName='Delta LIR, % GDP (Table S6)', append=T,row.names = F)
# xlsx::write.xlsx(sheet3,file = '../cepi_results/cepi_results.xlsx',sheetName='given return, SARS-X (Table S7)', append=T,row.names = F)
# xlsx::write.xlsx(sheet4,file = '../cepi_results/cepi_results.xlsx',sheetName='as % counterfactual (Table S8)', append=T,row.names = F)
xlsx::write.xlsx(sheet5,file = '../cepi_results/cepi_results.xlsx',sheetName='LIR, % GDP (BAU1)', append=F,row.names = F)
xlsx::write.xlsx(sheet5.5,file = '../cepi_results/cepi_results.xlsx',sheetName='LIR, % GDP (BAU2)', append=T,row.names = F)




# voi values
# voitab <- data.frame()
# for(ex in 1:length(exlist)){
#   voi <- c()
#   mi <- c()
#   subabscosttab <- subset(abscosttab,return==tens[ex])
#   for(j in 2:length(slvaluelist)){
#     for(refsl in 1:(j-1)){
#       newcolname <- paste0('val',refsl,'to',j)
#       sourcesj <- log(subabscosttab$cost)
#       y <- subabscosttab[[newcolname]]
#       vary <- var(y) 
#       sourcesj <- subabscosttab[,1,drop=F]
#       fittedvalues <- evppifit(y,sourcesj,pars='cost',method='earth')
#       mi[length(mi)+1] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
#       (voi[length(voi)+1] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100)
#       voitab <- rbind(voitab,c(rownames(absvalues)[ex],
#                                scenario_names[j],
#                                scenario_names[refsl],round(voi[length(voi)])))
#     }
#   }
#   print(round(voi))
# }
# colnames(voitab) <- c('Frequency','To','From scenario','VOI')
# saveRDS(reshape2::dcast(voitab,formula=`Frequency`+`From scenario`~To,value.var = 'VOI',fill=''),'results/valuevoi.Rds')


# resample_deathspermil <- sample(x=deathspermil,size=nsamples,replace = T,prob = yvals)
# ggplot() + 
#     geom_violin(aes(x=1,y=resample_deathspermil),position="dodge", scale='width')


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
# j <- nScen; refsl <- 1
# scatter <- data.frame(policy=apply(sapply(1:ncol(sampleorders[[1]]),function(x){
#   ig <- rep(1:3,times=samplefracs)[x]
#   orderedchoices[[ig]][sampleorders[[1]][,x]]
# }),1,function(y)sum(y=='No Closures'))/sum(samplefracs)*100,
# deaths=rowSums(deathsamples[[1]])/(50*1e6*ncountries)*1e3,
# pcyll = rowSums(slvaluelist[[j]][[refsl]]$ylls)/rowSums(slvaluelist[[j]][[refsl]]$slvalue)*100,
# weightedyllpc = apply(
#   sapply(1:ncol(sampleorders[[1]]),function(x){
#     ig <- rep(1:3,times=samplefracs)[x]
#     orderedyll[[ig]][sampleorders[[1]][,x]]
#   }),1,sum)/rowSums(costsamples[[1]])*100
# )
# cyl_labels <- c("weightedyllpc" = "LIR", "pcyll" = "Delta*'LIR'")
# (pcyllplot <- ggplot(reshape2::melt(scatter,id.vars=c('deaths','policy'))) + 
#   geom_vline(xintercept = quantile((scatter$deaths),c(1,3)/4),linewidth=1.5,colour='grey') +
#   geom_point(aes(x=(deaths),y=value,colour=policy)) +
#     facet_grid(~ variable, labeller = as_labeller(cyl_labels, default = label_parsed)) + 
#   # facet_grid(~factor(variable,levels=c('weightedyllpc','pcyll'),labels=c('LIR',deparse(bquote("Delta"))))) +
#   scale_x_continuous(transform='log',breaks=c(.1,1,10)) +
#   theme_bw(base_size=18) + 
#   labs(x='Deaths per thousand',y='YLL percent',colour='% No Closures'))
# ggsave(pcyllplot,filename='results/pcYLL.png',width=10,height=6)



