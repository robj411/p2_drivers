## start ################################

library(earth)
library(ggplot2)
library(infotheo)
library(data.table)
library(viridis)
library(GGally)
library(splines)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

if(Sys.info()[['sysname']]=='Linux'){
  setwd('~/projects/DAEDALUS/Daedalus-P2-Dashboard/code')
}else{
  setwd('C:/Users/rj411/OneDrive - Imperial College London/p2_drivers/code')
}

## functions and variables #####################################

strategies <- c('No Closures','School Closures','Economic Closures','Elimination')
income_levels <- c('LLMIC','UMIC','HIC')
vaccination_levels <- c(365,100)
bpsv_levels <- c(0,1)
result_cols <- c('Cost','YLL','School','GDP_loss')

#

evppivar_npreg <- function (outputs, inputs, pars, method, verbose, ...) {
  fitted <- voi:::fitted_npreg_fn(method)(y = outputs, inputs = inputs, 
                                    pars = pars, verbose = verbose, ...)
  list(data.frame(evppi = var(fitted)),fitted)
}

evppivar <- function (outputs, inputs, pars = NULL, method = NULL, nsim = NULL, 
                      verbose = TRUE, ...) {
  inputs <- voi:::check_inputs(inputs, iname = deparse(substitute(inputs)))
  voi:::check_outputs_vector(outputs, inputs)
  if (is.list(pars)) {
    return(evppivar_list(outputs = outputs, inputs = inputs, 
                         pars = pars, method = method, nsim = nsim, verbose = verbose, 
                         ...))
  }
  pars <- voi:::check_pars(pars, inputs)
  opts <- list(...)
  if (is.null(method)) 
    method <- voi:::default_evppi_method(pars)
  if (is.null(nsim)) 
    nsim <- nrow(inputs)
  outputs <- outputs[1:nsim]
  inputs <- inputs[1:nsim, , drop = FALSE]
  if (method %in% voi:::npreg_methods) {
    rese <- evppivar_npreg(outputs = outputs, inputs = inputs, 
                           pars = pars, method = method, verbose = verbose, 
                           ...)
  }
  else stop("Other methods not implemented yet")
  res <- cbind(pars = paste(pars, collapse = ","), rese[[1]])
  res
}

evppifit <- function (outputs, inputs, pars = NULL, method = NULL, nsim = NULL, 
                      verbose = TRUE, ...) {
  inputs <- voi:::check_inputs(inputs, iname = deparse(substitute(inputs)))
  voi:::check_outputs_vector(outputs, inputs)
  if (is.list(pars)) {
    return(evppivar_list(outputs = outputs, inputs = inputs, 
                         pars = pars, method = method, nsim = nsim, verbose = verbose, 
                         ...))
  }
  pars <- voi:::check_pars(pars, inputs)
  opts <- list(...)
  if (is.null(method)) 
    method <- voi:::default_evppi_method(pars)
  if (is.null(nsim)) 
    nsim <- nrow(inputs)
  outputs <- outputs[1:nsim]
  inputs <- inputs[1:nsim, , drop = FALSE]
  if (method %in% voi:::npreg_methods) {
    rese <- evppivar_npreg(outputs = outputs, inputs = inputs, 
                           pars = pars, method = method, verbose = verbose, 
                           ...)
  }
  else stop("Other methods not implemented yet")
  res <- rese[[2]]
  res
}




multisource <- list('Importation time'="t_import",
                    'Probability symptomatic'="ps",
                    'Asymptomatic infectiousness'="red",
                    "Frac presymptomatic"="frac_presymptomatic",
                    'Latent period'="Tlat",
                    'Asymptomatic period'="Tay",
                    'Infectious period'="Tsr",
                    'Time to hospitalisation'="Tsh",                              
                    'Time to discharge'="Threc",
                    'Time to death'="Thd",
                    "R0"="R0",                               
                    "R0 + seed size"=c("R0",'seedsize'),                               
                    "Candidate infectees"="CI",                               
                    "beta"="beta",       
                    'Doubling time'="Td",
                    'Generation time'="generation_time",
                    'Hosp + Duration + R0'=c('R0','Threc','Hmax'),
                    'mean IHR'=c('IHR_mean'),
                    'mean IHR + R0'=c('R0','IHR_mean'),
                    'mean IFR'=c('IFR_mean'),
                    'mean IFR + R0'=c('R0','IFR_mean'),
                    'mean IHR + mean IFR + R0'=c('R0','IFR_mean','IHR_mean'),
                    'mean HFR'=c('HFR_mean'),
                    'mean HFR + R0'=c('R0','HFR_mean'),
                    
                    'Hospital response'="Hres",
                    'Response time'="Tres",
                    "Response time quantile"="response_time_quantile",
                    'Testing rate'="trate",
                    "Self isolation compliance"="self_isolation_compliance",
                    'Fraction infectiousness averted'=c("frac_sym_infectiousness_averted","frac_presym_infectiousness_averted","frac_asym_infectiousness_averted"),
                    'Hospital threshold'="hosp_release_trigger",
                    'Hospital capacity'="Hmax",                              
                    'Social distancing function'=c("sd_mandate_coef","sd_death_coef","sd_baseline"),
                    "SD baseline"="sd_baseline",
                    "SD death coefficient"="sd_death_coef" ,                    
                    "SD mandate coefficient"="sd_mandate_coef",
                    
                    'GDP'="gdp" ,                              
                    'International tourism'="frac_tourism_international",
                    'Tourism'=c('frac_tourism_international','obj32'), #Food_sector
                    "Labour share"="labsh",
                    "Employment rate"='employmentrate',
                    "Remote teaching effectiveness"="remote_teaching_effectiveness",     
                    'VSY'="vsy",
                    'VLY'="vly",
                    'Valuations'=c('vly','vsy'),
                    'Education factors'=c('NNs47','vsy'),
                    'Age groups'=c('NNs47','NNs46','NNs49'),
                    "Remote working"="remote_quantile",
                    'Work frac'='work_frac',
                    'Nursery frac'='school1_frac',
                    'School frac'='school2_frac',
                    'Pupil-teacher ratio'='pupil_teacher_ratio',
                    'School/work contacts'=c('pupil_teacher_ratio','work_frac','school1_frac','school2_frac'),
                    "Hospitality frac"=c("hospitality_frac1","hospitality_frac2","hospitality_frac3","hospitality_frac4"),
                    'Timing'=c('t_import','Tres','Hres','response_time_quantile'),
                    # 'Testing + R0'=c('R0','trate',"frac_sym_infectiousness_averted","frac_presym_infectiousness_averted","frac_asym_infectiousness_averted"),
                    'Testing'=c('trate',"frac_sym_infectiousness_averted","frac_presym_infectiousness_averted","frac_asym_infectiousness_averted"),
                    'Seed size'='seedsize'
                    )


#
## saving 1 ###############################

topresults <- list()
choicestab <- c()
# income_levels <- 'HIC'
for(bl in 1:length(bpsv_levels)){
  topresults[[bl]] <- list()
  bpsv <- bpsv_levels[bl]
  for(v in 1:length(vaccination_levels)){
    vaccination_level <- vaccination_levels[v]
    allresults <- data.frame()
    for (i in 1:length(income_levels)){
      income_level <- income_levels[i];
      inputtab <- read.csv(paste0('results/inputs_',income_level,'.csv'),header=T);
      for (k in 1:length(strategies)){
        inp3 <- strategies[k];
        results <- read.csv(paste0('results/outputs_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
        results$policy <- inp3
        results$igroup <- income_level
        results$samplei <- 1:nrow(results)
        results <- cbind(results,inputtab)
        allresults <- rbind(allresults,results)
      }
    }
    setDT(allresults)
    allresults[,Costrnd:=Cost*(1+runif(nrow(allresults))*1e-8)]
    allresults[,mincost:=Costrnd==min(Costrnd),by=.(igroup,samplei)]
    choices <- allresults[,.(N=sum(mincost)),by=.(igroup,policy)]
    choices$vaccination_level <- vaccination_level
    choices$bpsv <- bpsv
    choicestab <- rbind(choicestab,choices)
    allresults[mincost==T,mean(Cost/gdp),by=.(igroup)]
    
    evalues <- allresults[,mean(Cost/gdp),by=.(igroup,policy)]
    topchoices <- evalues[,list(policy[which.min(V1)],round(V1[which.min(V1)]*100)),by=igroup]
    print(topchoices)
    # allresults[,keeprow:=policy==topchoices$V1[which(topchoices$igroup==igroup)],by=igroup]
    
    ##!! decision under uncertainty, or decision under certainty?
    # topresults[[bl]][[v]] <- subset(allresults,keeprow)
    topresults[[bl]][[v]] <- subset(allresults,mincost)
    print(dim(topresults[[bl]][[v]]))
    setorder(topresults[[bl]][[v]],igroup,samplei)
    
    if(v>1|bl>1){
      
      difftab <- copy(topresults[[1]][[1]])
      # difftab[,policy:=NULL]
      difftab[,Cost:=Cost-topresults[[bl]][[v]]$Cost]
      difftab[,YLL:=YLL-topresults[[bl]][[v]]$YLL]
      difftab[,School:=School-topresults[[bl]][[v]]$School]
      difftab[,GDP_loss:=GDP_loss-topresults[[bl]][[v]]$GDP_loss]
      difftab[,scen_Exit_wave:=topresults[[bl]][[v]]$Exit_wave]
      difftab[,scen_Mitigated_deaths:=topresults[[bl]][[v]]$Mitigated_deaths]
      
      saveRDS(difftab,paste0('results/difftab_',bpsv,'_',vaccination_level,'.Rds'))
      
      print(difftab[,round(mean(Cost/gdp*100)),by=.(igroup)])
      
      pctab <- copy(topresults[[1]][[1]])
      pctab[,policy:=NULL]
      pctab[,Cost:=(Cost-topresults[[bl]][[v]]$Cost)/Cost]
      pctab[,YLL:=(YLL-topresults[[bl]][[v]]$YLL)/YLL]
      pctab[,School:=(School-topresults[[bl]][[v]]$School)/School]
      pctab[,GDP_loss:=(GDP_loss-topresults[[bl]][[v]]$GDP_loss)/GDP_loss]
      
      saveRDS(pctab,paste0('results/pctab_',bpsv,'_',vaccination_level,'.Rds'))
    }
  }
}
dim(topresults[[1]][[1]])[1]/3
plotdur <- list()
for(i in 1:length(bpsv_levels)){
  plotdur[[i]] <- copy(topresults[[i]][[1]])
  plotdur[[i]]$bpsv <- c('No BPSV','BPSV')[i]
}
plotdur <- do.call(rbind,plotdur)
mplotdur <- melt(plotdur[,.(bpsv,igroup,End_mitigation,End_simulation)],measure.vars=c('End_mitigation','End_simulation'))
ggplot(mplotdur) + geom_density(aes(x=value/365,y=..density..),fill='darkorange',colour='midnightblue',alpha=.75,linewidth=2) +
  facet_grid(~factor(variable,labels=c('Mitigation end','Simulation end')),scales='free') +
  theme_bw(base_size=15) +
  scale_x_continuous(limits=c(0,NA)) +
  scale_y_continuous(expand=c(0,NA)) +
  labs(x='Day',y='Density')
unique(subset(topresults[[1]][[1]],End_simulation>2000)$samplei)

params <- unlist(multisource)    
params[!params%in%colnames(allresults)]

endhosp <- list()
for(j in 1:length(vaccination_levels)){
  endhosp[[j]] <- list()
  for(i in 1:length(bpsv_levels)){
    endhosp[[j]][[i]] <- copy(topresults[[i]][[j]])
    endhosp[[j]][[i]]$bpsv <- c('No BPSV','BPSV')[i]
  }
  endhosp[[j]] <- do.call(rbind,endhosp[[j]])
  endhosp[[j]]$sarsx <- c(365,100)[j]
}
endhosp <- do.call(rbind,endhosp)
dcast(endhosp,formula=igroup+samplei~bpsv+sarsx,value.var='End_hosp')

unique(subset(endhosp,End_hosp>1000)[,.(igroup,policy,samplei)])

## negatives ###########################

dispcols <- c('Cost','GDP_loss','YLL','School','igroup','gdp',
                                      'policy','samplei','Exit_wave','scen_Exit_wave')
nneg <- 0
for(vaccination_level in vaccination_levels){
  for(bpsv in bpsv_levels){
    if(vaccination_level==100|bpsv==1){
      print(paste0('results/difftab_',bpsv,'_',vaccination_level,'.Rds'))
      difftab <- readRDS(paste0('results/difftab_',bpsv,'_',vaccination_level,'.Rds'))
      for(income_level in income_levels){
        subtab <- subset(difftab,igroup==income_level)
        subtab$sample <- 1:nrow(subtab)
        print(income_level)
        print(subset(subtab,Cost/gdp< -.05 & scen_Mitigated_deaths > Mitigated_deaths )[,..dispcols])
        nneg <- nneg + nrow(subset(subtab,Cost< 0 ))
      }
    }
  }
} 
      
## saving 2 ###############################

topresults <- list()
choicestab <- c()
ilistvoi <- ilistmi <- list()
for(bl in 1:length(bpsv_levels)){
  topresults[[bl]] <- list()
  bpsv <- bpsv_levels[bl]
  for(v in 1:length(vaccination_levels)){
    vaccination_level <- vaccination_levels[v]
    allresults <- data.frame()
    for (i in 1:length(income_levels)){
      income_level <- income_levels[i];
      inputtab <- read.csv(paste0('results/inputs_',income_level,'.csv'),header=T);
      for (k in 1:length(strategies)){
        inp3 <- strategies[k];
        results <- read.csv(paste0('results/outputs_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
        results$policy <- inp3
        results$igroup <- income_level
        results$samplei <- 1:nrow(results)
        results <- cbind(results,inputtab)
        allresults <- rbind(allresults,results)
      }
    }
    setDT(allresults)
    allresults[,Costrnd:=Cost*(1+runif(nrow(allresults))*1e-8)]
    allresults[,mincost:=Costrnd==min(Costrnd),by=.(igroup,samplei)]
    choices <- allresults[,.(N=sum(mincost)),by=.(igroup,policy)]
    choices$vaccination_level <- vaccination_level
    choices$bpsv <- bpsv
    choicestab <- rbind(choicestab,choices)
    allresults[mincost==T,mean(Cost/gdp),by=.(igroup)]
    
    evalues <- allresults[,mean(Cost/gdp),by=.(igroup,policy)]
    topchoices <- evalues[,policy[which.min(V1)],by=igroup]
    print(topchoices)
    allresults[,keeprow:=policy==topchoices$V1[which(topchoices$igroup==igroup)],by=igroup]
    
    ##!! decision under uncertainty, or decision under certainty?
    # topresults[[bl]][[v]] <- subset(allresults,keeprow)
    topresults[[bl]][[v]] <- subset(allresults,mincost)
    setorder(topresults[[bl]][[v]],igroup,samplei)
    
    if(v>1|bl>1){
      
      difftab <- copy(topresults[[bl]][[v]])
      difftab[,keeprow:=NULL]
      difftab[,policy:=NULL]
      difftab[,Cost:=-Cost+topresults[[1]][[1]]$Cost]
      
      
      voilist <- list()
      milist <- list()
      for (il in 1:length(income_levels)){
        income_level <- income_levels[il]
        results <- as.data.frame(subset(difftab,igroup==income_level&R0>1))
        results$igroup <- NULL
        sourcelist <- list()
        for(src in 1:length(multisource)) {
          sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]],drop=F]
          if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
        }
        
        sourcemat <- results
        y <- results$Cost/sourcemat$gdp
        
          voi <- c()
          mi <- c()
          vary <- var(y) 
          for(j in 1:length(sourcelist)){
            sourcesj <- sourcelist[[j]]
            fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj),method='gam')
            mi[j] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
            # voi[j+ncol(sourcemat)] <- voi::evppivar(y,sourcesj,pars=colnames(sourcesj))[2]/vary*100
            voi[j] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
          }
          milist[[il]] <- mi
          voilist[[il]] <- voi

        
      }
      
      voitab <- do.call(rbind,voilist)
      colnames(voitab) <- c(names(multisource))
      rownames(voitab) <- paste0(vaccination_level,', ',c('No ','')[bpsv+1],'BPSV: ',income_levels)
      mitab <- do.call(rbind,milist)
      colnames(mitab) <- c(names(multisource))
      rownames(mitab) <- paste0(vaccination_level,', ',c('No ','')[bpsv+1],'BPSV: ',income_levels)
      
      ilistvoi[[length(ilistvoi)+1]] <- voitab
      ilistmi[[length(ilistmi)+1]] <- mitab
      
      
    }
    
    {
    # trdf <- topresults[[bl]][[v]]
    # sourcelist <- list()
    # for(src in 1:length(multisource)) {
    #   sourcelist[[src]] <- trdf[,colnames(trdf)%in%multisource[[src]],with=F]
    #   if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
    # }
    # 
    # mi <- c()
    # for(j in 1:length(sourcelist)){
    #   print(paste0(i,' ',names(multisource)[j]))
    #   sourcesj <- sourcelist[[j]]
    #   tryCatch({
    #     mi[j] <- infotheo::mutinformation(infotheo::discretize(sourcesj),(trdf$policy))
    #   },error = function(cond) {
    #     mi[j] <- 0
    #     # NA
    #   })
    # }
    # names(mi) <- names(multisource)
    # sort(mi)
  }
  }
}
voiall <- do.call(rbind,ilistvoi)
# roworder <- unlist(lapply(c("Cost","YLL","School","GDP loss"),
#                           function(x)which(grepl(paste0(x,':'),rownames(voiall)))))
# voiall <- voiall[roworder,]
miall <- do.call(rbind,ilistmi)
# miall <- miall[roworder,]

saveRDS(voiall,paste0('results/voidiff.Rds'))
saveRDS(miall,paste0('results/midiff.Rds'))

saveRDS(choicestab,paste0('results/choicestab.Rds'))

#
## income levels separately ####################################

get_voi_mi <- function(inp3,income_level,vaccination_level,bpsv){
  print(paste0('results/outputs_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'))
  pth <- 'C:/Users/rj411/OneDrive - Imperial College London/p2_drivers/code/results/'
  inputtab <- read.csv(paste0(pth,'inputs_',income_level,'.csv'),header=T);
  results <- read.csv(paste0(pth,'outputs_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
  results <- cbind(results,inputtab)
  sourcelist <- list()
  for(src in 1:length(multisource)) {
    sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]],drop=F]
    if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
  }
  
  firstreultcol <- which(colnames(results)==result_cols[1])
  
  colnames(results) <- gsub('_',' ',colnames(results))
  colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
  colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
  colnames(results)[colnames(results)=='Labour share'] <- 'Labour share of GVA'
  
  sourcemat <- results#[,1:(firstreultcol-1)]
  
  outcomes <- results[,firstreultcol:(firstreultcol+length(result_cols)-1)]
  for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/sourcemat$gdp
  
  voilist <- list()
  milist <- list()
  
  for(i in 1:ncol(outcomes)){
    voi <- c()
    mi <- c()
    y <- log(outcomes[,i])
    vary <- var(y) 
    for(j in 1:length(sourcelist)){
      print(paste0(i,' ',names(multisource)[j]))
      sourcesj <- sourcelist[[j]]
      tryCatch({
        fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj),method='gam')
        mi[j] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
        # voi[j+ncol(sourcemat)] <- voi::evppivar(y,sourcesj,pars=colnames(sourcesj))[2]/vary*100
        voi[j] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
      },error = function(cond) {
        mi[j] <- voi[j] <- 0
        # NA
      })
    }
    milist[[i]] <- mi
    voilist[[i]] <- voi
  }
  voitab <- do.call(rbind,voilist)
  colnames(voitab) <- c(names(multisource))
  rownames(voitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3,', ',vaccination_level,', ',bpsv)
  mitab <- do.call(rbind,milist)
  colnames(mitab) <- c(names(multisource))
  rownames(mitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3,', ',vaccination_level,', ',bpsv)
  list(voitab,mitab)
}

voimi <- list()
for(vl in 1:length(vaccination_levels)){
  vaccination_level <- vaccination_levels[vl]
  voimi[[vl]] <- list()
  for(bl in 1:length(bpsv_levels)){
    bpsv <- bpsv_levels[bl]
    voimi[[vl]][[bl]] <- foreach (ks = 1:length(strategies))%dopar%{
      inp3 <- strategies[ks]
      vmlist <- list()
      for (il in 1:length(income_levels)){
        income_level <- income_levels[il]
        vmlist[[il]] <- get_voi_mi(inp3,income_level,vaccination_level,bpsv)
      }
      vmlist
    }
  }
}

voilist <- milist <- list()
for(vl in 1:length(vaccination_levels)){
  vaccination_level <- vaccination_levels[vl]
  voilist[[vl]] <- milist[[vl]] <- list()
  for(bl in 1:length(bpsv_levels)){
    bpsv <- bpsv_levels[bl]
    voilist[[vl]][[bl]] <- milist[[vl]][[bl]] <- list()
    for (ks in 1:length(strategies)){
      voilist[[vl]][[bl]][[ks]] <- milist[[vl]][[bl]][[ks]] <- list()
      for (il in 1:length(income_levels)){
        voilist[[vl]][[bl]][[ks]][[il]] <- voimi[[vl]][[bl]][[ks]][[il]][[1]]
        milist[[vl]][[bl]][[ks]][[il]] <- voimi[[vl]][[bl]][[ks]][[il]][[2]]
      }
      voilist[[vl]][[bl]][[ks]] <- do.call(rbind,voilist[[vl]][[bl]][[ks]])
      milist[[vl]][[bl]][[ks]] <- do.call(rbind,milist[[vl]][[bl]][[ks]])
    }
    voilist[[vl]][[bl]] <- do.call(rbind,voilist[[vl]][[bl]])
    milist[[vl]][[bl]] <- do.call(rbind,milist[[vl]][[bl]])
    
    roworder <- unlist(lapply(c("Cost","YLL","School","GDP loss"),
                              function(x)which(grepl(paste0(x,':'),rownames(voilist[[vl]][[bl]])))))
    voilist[[vl]][[bl]] <- voilist[[vl]][[bl]][roworder,]
    milist[[vl]][[bl]] <- milist[[vl]][[bl]][roworder,]
    
    saveRDS(voilist[[vl]][[bl]],paste0('results/voi_',vaccination_level,'_',bpsv,'.Rds'))
    saveRDS(milist[[vl]][[bl]],paste0('results/mi_',vaccination_level,'_',bpsv,'.Rds'))
  }
}


#
## decision #####################################

if(F){
for(vl in 1:length(vaccination_levels)){
  vaccination_level <- vaccination_levels[vl]
  for(bl in 1:length(bpsv_levels)){
    bpsv <- bpsv_levels[bl]
    listout <- foreach (il = 1:length(income_levels))%dopar%{
      resultslist <- list()
      for (ks in 1:length(strategies)){
        inp3 <- strategies[ks];
        income_level <- income_levels[il];
        results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
        results <- subset(results,R0>1)
        firstreultcol <- which(colnames(results)=='Cost')
        resultslist[[ks]] <- results[,firstreultcol:(ncol(results))]
        for(i in 1:ncol(resultslist[[ks]])) resultslist[[ks]][,i] <- -resultslist[[ks]][,i]/results$gdp
        print(paste0('results/VOI_',inp3,'_',income_level,'.csv'))
      }
      
      sourcelist <- list()
      for(src in 1:length(multisource)) {
        sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
        if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
      }
      
      colnames(results) <- gsub('_',' ',colnames(results))
      colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
      colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
      colnames(results)[colnames(results)=='Labour share'] <- 'Labour share of GVA'
      
      sourcemat <- results[,1:(firstreultcol-1)]
      
      
      voilist <- list()
      
      for(i in 1:ncol(resultslist[[1]])){
        voi <- rep(0,ncol(sourcemat)+length(sourcelist))
        y <- do.call(cbind,lapply(resultslist,'[',i))
        colnames(y) <- paste0(colnames(y),1:ncol(y))
        keeprows <- rowSums(y < -.05)>0
        for(j in 1:ncol(sourcemat)){
          # model outcome as a function of input(s)
          sourcesj <- sourcemat[,j,drop=F]
          if(diff(range(sourcesj))>0)
            voi[j] <- voi::evppi(y[keeprows,],sourcesj[keeprows,,drop=F],pars=colnames(sourcesj))$evppi
        }
        for(j in 1:length(sourcelist)){
          sourcesj <- sourcelist[[j]]
          voi[j+ncol(sourcemat)] <- voi::evppi(y[keeprows,],sourcesj[keeprows,],pars=colnames(sourcesj),method='gam')$evppi
        }
        voilist[[i]] <- voi
      }
      
      voitab <- do.call(rbind,voilist)
      colnames(voitab) <- c(colnames(sourcemat),names(multisource))
      rownames(voitab) <- paste0(colnames(results)[firstreultcol:(ncol(results))],': ',income_level)
      
      # klistvoi[[il]] <- voitab
      voitab
    }
    
    voiall <- do.call(rbind,listout)
    roworder <- unlist(lapply(c("Cost","YLL","School","GDP loss"),
                              function(x)which(grepl(paste0(x,':'),rownames(voiall)))))
    voiall <- voiall[roworder,]
    saveRDS(voiall,paste0('results/decisionvoi',vaccination_level,'_',bpsv,'.Rds'))
  }
}
}



## log cost share plot ###################################

if(F){
library(scales)
ivoioutcomelist <- list()
for (il in 1:length(income_levels)){
  voioutcomelist <- list()
  for (ks in 1:length(strategies)){
    inp3 <- strategies[ks];
    income_level <- income_levels[il];
    results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'_BAU.csv'),header=T);
    
    firstreultcol <- which(colnames(results)=='Cost')
    outcomes <- results[,firstreultcol:(ncol(results))]
    
    # meltedres <- melt(outcomes[,2:4]/outcomes[,1]*100,value.name = 'value')
    for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/results$gdp*100
    for(i in 2:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/outcomes[,1]
    
    
    setDT(outcomes)
    outcomes[,log10cost:=log10(Cost)]
    x <- hist(outcomes$log10cost,plot = F, breaks=seq(-20,20,by=.15))
    breaks <- x$breaks
    outcomes[,costcat:=cut(log10cost,breaks,labels=1:(length(breaks)-1))]
    outcomes[,density:=x$density[costcat]/max(x$density),by=costcat]
    barplt <- outcomes[,.(School=mean(School)*mean(density),
                          GDP_loss=mean(GDP_loss)*mean(density),
                          YLL=mean(YLL)*mean(density)),by=costcat]
    barplt[,midpoint:=x$mids[costcat],by=costcat]
    
    
    meltbarplt <- reshape2::melt((barplt[,2:5]),value.name = 'value',measure.vars=c('School','YLL','GDP_loss'))
    
    # meltedres <- reshape2::melt((outcomes[,2:4]),value.name = 'value',measure.vars=c('School','YLL','GDP_loss'))
    meltbarplt$Income_group <- income_level
    meltbarplt$Policy <- inp3
    voioutcomelist[[ks]] <- meltbarplt
  }
  ivoioutcomelist[[il]] <- do.call(rbind,voioutcomelist)
}
voioutcomes <- do.call(rbind,ivoioutcomelist)
voioutcomes$Income_group <- factor(voioutcomes$Income_group,levels=c('LLMIC','MIC','HIC'))
voioutcomes$Policy <- factor(voioutcomes$Policy,levels=c('No Closures','School Closures','Economic Closures','Elimination'))
voioutcomes$variable <- factor(voioutcomes$variable,levels=c('School','YLL','GDP_loss'),labels=c('Education','YLL','GDP'))

ggplot(voioutcomes) + 
  annotate(geom = "rect",xmin = 10,xmax = 10^1.5,ymin = -Inf,ymax = +Inf,alpha = 0.2) +
  annotate(geom = "rect",xmin = 10^1.5,xmax = 100,ymin = -Inf,ymax = +Inf,alpha = 0.4) +
  annotate(geom = "rect",xmin = 100,xmax = 10^2.5,ymin = -Inf,ymax = +Inf,alpha = 0.6) +
  annotate(geom = "rect",xmin = 10^2.5,xmax = 1000,ymin = -Inf,ymax = +Inf,alpha = 0.8) +
  # geom_vline(aes(xintercept=100),colour='grey',size=1) +
  # geom_vline(aes(xintercept=10),colour='grey',size=0.5) +
  # geom_vline(aes(xintercept=1000),colour='grey',size=1.5) +
  geom_bar(aes(x=10^midpoint,y=value,fill=variable),stat="identity") +
  # geom_histogram(aes(x=value,colour=variable,fill=variable,y=..density..),alpha=.7)+#,position='identity') +
  facet_grid(Policy~Income_group,scales='free_y') +
  scale_fill_viridis(discrete=T, name="",option='inferno') +
  scale_colour_viridis(discrete=T, name="",option='inferno') +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits=c(10e-1,10e2),expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size = 18) +
  labs(x='Total cost, % of GDP',y='Relative frequency') +
  theme(legend.position = 'bottom',axis.text.x = element_text(angle = 60,  hjust=1,vjust=1)) -> p
ggsave(p,filename = paste0('logcostshare.pdf'),width=10,height=10)
}
    
