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
  setwd('~/overflow_dropbox/DAEDALUS/Daedalus-P2-Dashboard/code')
}else{
  setwd('C:/Users/rj411/OneDrive - Imperial College London/p2_drivers/code')
}

## functions and variables #####################################

strategies <- c('No Closures','School Closures','Economic Closures','Elimination')
income_levels <- c('LLMIC','UMIC','HIC')
vaccination_levels <- c(365,100)
bpsv_levels <- c(0,1)
result_cols <- c('Cost','dYLLs','School','GDP_loss')

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




multisource <- list('Sectors'=c('Agriculture','Food_sector'),
                    'Tourism'=c('International_tourism','Food_sector'),
                    'Valuations'=c('VLY','VSY'),
                    'Education factors'=c('School_age','VSY'),
                    'Age groups'=c('School_age','Working_age','Elders'),
                    'Working'=c('Working_age','Work_contacts'),#,'Remote_quantile'),
                    'Contacts'=c('Work_contacts','Hospitality_contacts','School_contacts'),
                    'Timing'=c('Response_time','Importation_time'),
                    'Timing + R0'=c('R0','Response_time','Importation_time'),
                    'Testing'=c('Test_rate','Fraction_infectiousness_averted'),
                    'Testing + R0'=c('R0','Test_rate','Fraction_infectiousness_averted'),
                    'Social distancing'=c('Social_distancing_baseline','Social_distancing_death','Social_distancing_mandate'),
                    #'BMI RR' = c('BMI','BMI_hospitalisation','BMI_infection','BMI_death'),
                    'IHR + R0'=c('Mean_IHR', 'R0'),
                    'HFR + R0'=c('Mean_HFR', 'R0'),
                    'IFR + R0'=c('Mean_IFR', 'R0'),
                    'IHR + IFR + R0'=c('Mean_IHR','Mean_IFR', 'R0'),
                    'Duration + IHR + R0'=c('R0','Mean_IHR','Time_to_discharge'),
                    'Hosp + Duration + R0'=c('R0','Hospital_capacity','Time_to_discharge'),
                    'Hosp + IHR + R0'=c('R0','Mean_IHR','Hospital_capacity')#,
                    # 'SD rate + IHR + R0'=c('Social_distancing_rate','R0','Mean_IHR'),
                    # 'SD rate + Hosp + IHR + R0'=c('Social_distancing_rate','R0','Mean_IHR','Hospital_capacity')
                    #'SD rate + Hosp + IHR + IFR + R0'=c('Social_distancing_rate','R0','Mean_IHR','Mean_IFR','Hospital_capacity')
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
      for (k in 1:length(strategies)){
        inp3 <- strategies[k];
        income_level <- income_levels[i];
        results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
        results$strategy <- inp3
        results$igroup <- income_level
        results$samplei <- 1:nrow(results)
        allresults <- rbind(allresults,results)
      }
    }
    setDT(allresults)
    allresults[,Costrnd:=Cost*(1+runif(nrow(allresults))*1e-8)]
    allresults[,mincost:=Costrnd==min(Costrnd),by=.(igroup,samplei)]
    choices <- allresults[,.(N=sum(mincost)),by=.(igroup,strategy)]
    choices$vaccination_level <- vaccination_level
    choices$bpsv <- bpsv
    choicestab <- rbind(choicestab,choices)
    allresults[mincost==T,mean(Cost/GDP),by=.(igroup)]
    
    evalues <- allresults[,mean(Cost/GDP),by=.(igroup,strategy)]
    topchoices <- evalues[,strategy[which.min(V1)],by=igroup]
    print(topchoices)
    # allresults[,keeprow:=strategy==topchoices$V1[which(topchoices$igroup==igroup)],by=igroup]
    
    ##!! decision under uncertainty, or decision under certainty?
    # topresults[[bl]][[v]] <- subset(allresults,keeprow)
    topresults[[bl]][[v]] <- subset(allresults,mincost)
    setorder(topresults[[bl]][[v]],igroup,samplei)
    
    if(v>1|bl>1){
      
      difftab <- copy(topresults[[1]][[1]])
      difftab[,keeprow:=NULL]
      # difftab[,strategy:=NULL]
      difftab[,Cost:=Cost-topresults[[bl]][[v]]$Cost]
      difftab[,dYLLs:=dYLLs-topresults[[bl]][[v]]$dYLLs]
      difftab[,School:=School-topresults[[bl]][[v]]$School]
      difftab[,GDP_loss:=GDP_loss-topresults[[bl]][[v]]$GDP_loss]
      difftab[,scen_Exit_wave:=topresults[[bl]][[v]]$Exit_wave]
      
      saveRDS(difftab,paste0('results/difftab_',bpsv,'_',vaccination_level,'.Rds'))
      
      pctab <- copy(topresults[[1]][[1]])
      pctab[,keeprow:=NULL]
      pctab[,strategy:=NULL]
      pctab[,Cost:=(Cost-topresults[[bl]][[v]]$Cost)/Cost]
      pctab[,dYLLs:=(dYLLs-topresults[[bl]][[v]]$dYLLs)/dYLLs]
      pctab[,School:=(School-topresults[[bl]][[v]]$School)/School]
      pctab[,GDP_loss:=(GDP_loss-topresults[[bl]][[v]]$GDP_loss)/GDP_loss]
      
      saveRDS(pctab,paste0('results/pctab_',bpsv,'_',vaccination_level,'.Rds'))
    }
  }
}
 

params <- unlist(multisource)    
params[!params%in%colnames(allresults)]

## negatives ###########################

dispcols <- colnames(allresults)%in%c('Cost','GDP_loss','Deaths','School','igroup','strategy','samplei','Exit_wave','Scen_Exit_wave')
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
        print(subset(subtab,Cost< 0 & !(scen_Exit_wave>(1-1e-2)) )[,..dispcols])
        nneg <- nneg + nrow(subset(subtab,Cost< 0 ))
      }
    }
  }
} 
      
## saving 2 ###############################

topresults <- list()
choicestab <- c()
for(bl in 1:length(bpsv_levels)){
  topresults[[bl]] <- list()
  bpsv <- bpsv_levels[bl]
  for(v in 1:length(vaccination_levels)){
    vaccination_level <- vaccination_levels[v]
    allresults <- data.frame()
    for (i in 1:length(income_levels)){
      for (k in 1:length(strategies)){
        inp3 <- strategies[k];
        income_level <- income_levels[i];
        results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
        results$strategy <- inp3
        results$igroup <- income_level
        results$samplei <- 1:nrow(results)
        allresults <- rbind(allresults,results)
      }
    }
    setDT(allresults)
    allresults[,Costrnd:=Cost*(1+runif(nrow(allresults))*1e-8)]
    allresults[,mincost:=Costrnd==min(Costrnd),by=.(igroup,samplei)]
    choices <- allresults[,.(N=sum(mincost)),by=.(igroup,strategy)]
    choices$vaccination_level <- vaccination_level
    choices$bpsv <- bpsv
    choicestab <- rbind(choicestab,choices)
    allresults[mincost==T,mean(Cost/GDP),by=.(igroup)]
    
    evalues <- allresults[,mean(Cost/GDP),by=.(igroup,strategy)]
    topchoices <- evalues[,strategy[which.min(V1)],by=igroup]
    print(topchoices)
    allresults[,keeprow:=strategy==topchoices$V1[which(topchoices$igroup==igroup)],by=igroup]
    
    ##!! decision under uncertainty, or decision under certainty?
    # topresults[[bl]][[v]] <- subset(allresults,keeprow)
    topresults[[bl]][[v]] <- subset(allresults,mincost)
    setorder(topresults[[bl]][[v]],igroup,samplei)
    
    if(v>1|bl>1){
      
      difftab <- copy(topresults[[bl]][[v]])
      difftab[,keeprow:=NULL]
      difftab[,strategy:=NULL]
      difftab[,Cost:=-Cost+topresults[[1]][[1]]$Cost]
      difftab[,dYLLs:=-dYLLs+topresults[[1]][[1]]$dYLLs]
      difftab[,School:=-School+topresults[[1]][[1]]$School]
      difftab[,GDP_loss:=-GDP_loss+topresults[[1]][[1]]$GDP_loss]
      
      saveRDS(difftab,paste0('results/difftab_',bpsv,'_',vaccination_level,'.Rds'))
      
      pctab <- copy(topresults[[bl]][[v]])
      pctab[,keeprow:=NULL]
      pctab[,strategy:=NULL]
      pctab[,Cost:=(-Cost+topresults[[1]][[1]]$Cost)/Cost]
      pctab[,dYLLs:=(-dYLLs+topresults[[1]][[1]]$dYLLs)/dYLLs]
      pctab[,School:=(-School+topresults[[1]][[1]]$School)/School]
      pctab[,GDP_loss:=(-GDP_loss+topresults[[1]][[1]]$GDP_loss)/GDP_loss]
      
      saveRDS(pctab,paste0('results/pctab_',bpsv,'_',vaccination_level,'.Rds'))
      
      ilistvoi <- ilistmi <- list()
      for (il in 1:length(income_levels)){
        income_level <- income_levels[il]
        results <- as.data.frame(subset(difftab,igroup==income_level&R0>1))
        results$igroup <- NULL
        sourcelist <- list()
        for(src in 1:length(multisource)) {
          sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
          if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
        }
        
        firstreultcol <- which(colnames(results)==result_cols[1])
        
        colnames(results) <- gsub('_',' ',colnames(results))
        colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
        colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
        colnames(results)[colnames(results)=='Labour share'] <- 'Labour share of GVA'
        
        sourcemat <- results[,1:(firstreultcol-1)]
        
        outcomes <- results[,firstreultcol:(firstreultcol+length(result_cols)-1)]
        for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/sourcemat$GDP
        
        voilist <- list()
        milist <- list()
        
        for(i in 1:ncol(outcomes)){
          voi <- c()
          mi <- c()
          y <- outcomes[,i]
          vary <- var(y) 
          for(j in 1:ncol(sourcemat)){
            # model outcome as a function of input(s)
            sourcesj <- sourcemat[,j]
            max_degree <- ifelse(is.vector(sourcesj),1,ncol(sourcesj))
            model <- earth::earth(y ~ sourcesj, degree=min(4,max_degree))
            # compute evppi as percentage
            voi[j] <- (vary - mean((y - model$fitted) ^ 2)) / vary * 100
            mi[j] <- infotheo::mutinformation(infotheo::discretize(sourcesj),infotheo::discretize(y))
          }
          for(j in 1:length(sourcelist)){
            sourcesj <- sourcelist[[j]]
            fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj),method='gam')
            mi[j+ncol(sourcemat)] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
            # voi[j+ncol(sourcemat)] <- voi::evppivar(y,sourcesj,pars=colnames(sourcesj))[2]/vary*100
            voi[j+ncol(sourcemat)] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
          }
          milist[[i]] <- mi
          voilist[[i]] <- voi
        }
        
        voitab <- do.call(rbind,voilist)
        colnames(voitab) <- c(colnames(sourcemat),names(multisource))
        rownames(voitab) <- paste0(colnames(outcomes),': ',income_level)
        mitab <- do.call(rbind,milist)
        colnames(mitab) <- c(colnames(sourcemat),names(multisource))
        rownames(mitab) <- paste0(colnames(outcomes),': ',income_level)
        
        ilistvoi[[il]] <- voitab
        ilistmi[[il]] <- mitab
      }
      
      
      voiall <- do.call(rbind,ilistvoi)
      roworder <- unlist(lapply(c("Cost","dYLLs","School","GDP loss"),
                                function(x)which(grepl(paste0(x,':'),rownames(voiall)))))
      voiall <- voiall[roworder,]
      miall <- do.call(rbind,ilistmi)
      miall <- miall[roworder,]
      
      saveRDS(voiall,paste0('results/voidiff_',bpsv,'_',vaccination_level,'.Rds'))
      saveRDS(miall,paste0('results/midiff_',bpsv,'_',vaccination_level,'.Rds'))
    }
  }
}

saveRDS(choicestab,paste0('results/choicestab.Rds'))

#
## income levels separately ####################################

get_voi_mi <- function(inp3,income_level,vaccination_level,bpsv){
  results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'),header=T);
  results <- subset(results,R0>1)
  print(paste0('results/VOI_',inp3,'_',income_level,'_',vaccination_level,'_',bpsv,'.csv'))
  sourcelist <- list()
  for(src in 1:length(multisource)) {
    sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
    if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
  }
  
  firstreultcol <- which(colnames(results)==result_cols[1])
  
  colnames(results) <- gsub('_',' ',colnames(results))
  colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
  colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
  colnames(results)[colnames(results)=='Labour share'] <- 'Labour share of GVA'
  
  sourcemat <- results[,1:(firstreultcol-1)]
  
  outcomes <- results[,firstreultcol:(firstreultcol+length(result_cols)-1)]
  for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/sourcemat$GDP
  
  voilist <- list()
  milist <- list()
  
  for(i in 1:ncol(outcomes)){
    voi <- c()
    mi <- c()
    y <- log(outcomes[,i])
    vary <- var(y) 
    for(j in 1:ncol(sourcemat)){
      # model outcome as a function of input(s)
      sourcesj <- sourcemat[,j]
      max_degree <- ifelse(is.vector(sourcesj),1,ncol(sourcesj))
      model <- earth::earth(y ~ sourcesj, degree=min(4,max_degree))
      # compute evppi as percentage
      voi[j] <- (vary - mean((y - model$fitted) ^ 2)) / vary * 100
      mi[j] <- infotheo::mutinformation(infotheo::discretize(sourcesj),infotheo::discretize(y))
    }
    for(j in 1:length(sourcelist)){
      sourcesj <- sourcelist[[j]]
      fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj),method='gam')
      mi[j+ncol(sourcemat)] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
      # voi[j+ncol(sourcemat)] <- voi::evppivar(y,sourcesj,pars=colnames(sourcesj))[2]/vary*100
      voi[j+ncol(sourcemat)] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
    }
    milist[[i]] <- mi
    voilist[[i]] <- voi
  }
  voitab <- do.call(rbind,voilist)
  colnames(voitab) <- c(colnames(sourcemat),names(multisource))
  rownames(voitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3,', ',vaccination_level,', ',bpsv)
  mitab <- do.call(rbind,milist)
  colnames(mitab) <- c(colnames(sourcemat),names(multisource))
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
    
    roworder <- unlist(lapply(c("Cost","dYLLs","School","GDP loss"),
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
        for(i in 1:ncol(resultslist[[ks]])) resultslist[[ks]][,i] <- -resultslist[[ks]][,i]/results$GDP
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
    roworder <- unlist(lapply(c("Cost","dYLLs","School","GDP loss"),
                              function(x)which(grepl(paste0(x,':'),rownames(voiall)))))
    voiall <- voiall[roworder,]
    saveRDS(voiall,paste0('results/decisionvoi',vaccination_level,'_',bpsv,'.Rds'))
  }
}
}



## log cost share plot ###################################

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
    for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/results$GDP*100
    for(i in 2:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/outcomes[,1]
    
    
    setDT(outcomes)
    outcomes[,log10cost:=log10(Cost)]
    x <- hist(outcomes$log10cost,plot = F, breaks=seq(-20,20,by=.15))
    breaks <- x$breaks
    outcomes[,costcat:=cut(log10cost,breaks,labels=1:(length(breaks)-1))]
    outcomes[,density:=x$density[costcat]/max(x$density),by=costcat]
    barplt <- outcomes[,.(School=mean(School)*mean(density),
                          GDP_loss=mean(GDP_loss)*mean(density),
                          dYLLs=mean(dYLLs)*mean(density)),by=costcat]
    barplt[,midpoint:=x$mids[costcat],by=costcat]
    
    
    meltbarplt <- reshape2::melt((barplt[,2:5]),value.name = 'value',measure.vars=c('School','dYLLs','GDP_loss'))
    
    # meltedres <- reshape2::melt((outcomes[,2:4]),value.name = 'value',measure.vars=c('School','dYLLs','GDP_loss'))
    meltbarplt$Income_group <- income_level
    meltbarplt$Strategy <- inp3
    voioutcomelist[[ks]] <- meltbarplt
  }
  ivoioutcomelist[[il]] <- do.call(rbind,voioutcomelist)
}
voioutcomes <- do.call(rbind,ivoioutcomelist)
voioutcomes$Income_group <- factor(voioutcomes$Income_group,levels=c('LLMIC','MIC','HIC'))
voioutcomes$Strategy <- factor(voioutcomes$Strategy,levels=c('No Closures','School Closures','Economic Closures','Elimination'))
voioutcomes$variable <- factor(voioutcomes$variable,levels=c('School','dYLLs','GDP_loss'),labels=c('Education','YLLs','GDP'))

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
  facet_grid(Strategy~Income_group,scales='free_y') +
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

    
