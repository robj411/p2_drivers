library(earth)
library(ggplot2)
library(infotheo)
library(data.table)
library(viridis)

if(Sys.info()[['sysname']]=='Linux'){
  setwd('~/overflow_dropbox/DAEDALUS/Daedalus-P2-Dashboard/code')
}else{
  setwd('C:/Users/rj411/OneDrive - Imperial College London/p2_drivers/code')
}

## voi #####################################

strategies <- c('No Closures','School Closures','Economic Closures','Elimination')
income_levels <- c('LLMIC','MIC','HIC')

#
## functions ##########

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

## compute ########################

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
registerDoParallel(cores=4)


ilistvoi <- list()
ilistmi <- list()
multisource <- list('Sectors'=c('Agriculture','Food_sector'),
                    'Tourism'=c('International_tourism','Food_sector'),
                    'Education factors'=c('School_age','Internet','Labour_share'),
                    'Age groups'=c('School_age','Working_age','Elders'),
                    'Working'=c('Working_age','Work_contacts'),
                    'Contacts'=c('Work_contacts','Hospitality_contacts','Community_contacts','School_contacts'),
                    'Testing'=c('Test_rate','Test_start','Self_isolation_compliance'),
                    'Social distancing'=c('Social_distancing_rate','Social_distancing_max'),
                    'BMI RR' = c('BMI','BMI_hospitalisation','BMI_infection','BMI_death'),
                    'IHR + R0'=c('Mean_IHR', 'R0'),
                    'HFR + R0'=c('Mean_HFR', 'R0'),
                    'IFR + R0'=c('Mean_IFR', 'R0'),
                    'IHR + IFR + R0'=c('Mean_IHR','Mean_IFR', 'R0'),
                    'Duration + IHR + R0'=c('R0','Mean_IHR','Time_to_discharge'),
                    'Hosp + Duration + R0'=c('R0','Hospital_capacity','Time_to_discharge'),
                    'Hosp + IHR + R0'=c('R0','Mean_IHR','Hospital_capacity'),
                    'SD rate + IHR + R0'=c('Social_distancing_rate','R0','Mean_IHR'),
                    'SD rate + Hosp + IHR + R0'=c('Social_distancing_rate','R0','Mean_IHR','Hospital_capacity')
                    # 'SD rate + Hosp + IHR + IFR + R0'=c('Social_distancing_rate','Work_contacts','R0','Mean_IHR','Mean_IFR','Hospital_capacity'),
                    )

## income levels separately ####################################

klistvoi <- list()
klistmi <- list()
listout <- foreach (ks = 1:length(strategies))%dopar%{
  for (il in 1:length(income_levels)){
    inp3 <- strategies[ks];
    income_level <- income_levels[il];
    results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'.csv'),header=T);
    print(paste0('results/VOI_',inp3,'_',income_level,'.csv'))
    sourcelist <- list()
    for(src in 1:length(multisource)) {
      sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
      if(ncol(sourcelist[[src]])<length(multisource[[src]])) print(src)
    }
    
    firstreultcol <- which(colnames(results)=='Cost')
    
    colnames(results) <- gsub('_',' ',colnames(results))
    colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
    colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
    colnames(results)[colnames(results)=='Labour share'] <- 'Labour share of GVA'
    
    sourcemat <- results[,1:(firstreultcol-1)]
    
    outcomes <- results[,firstreultcol:(ncol(results))]
    for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/sourcemat$GDP
    # for(i in which(colnames(outcomes)%in%c('Cost','Deaths')))outcomes[,i] <- log(outcomes[,i])
    
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
    rownames(voitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3)
    mitab <- do.call(rbind,milist)
    colnames(mitab) <- c(colnames(sourcemat),names(multisource))
    rownames(mitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3)
    
    klistvoi[[il]] <- voitab
    klistmi[[il]] <- mitab
  }
  # ilistvoi[[ks]] <- do.call(rbind,klistvoi)
  # ilistmi[[ks]] <- do.call(rbind,klistmi)
  list(do.call(rbind,klistvoi), do.call(rbind,klistmi))
}

ilistvoi <- list(listout[[1]][[1]], listout[[2]][[1]], listout[[3]][[1]], listout[[4]][[1]])
ilistmi <- list(listout[[1]][[2]], listout[[2]][[2]], listout[[3]][[2]], listout[[4]][[2]])

voiall <- do.call(rbind,ilistvoi)
roworder <- unlist(lapply(c("Cost","Deaths","School","GDP loss"),
                          function(x)which(grepl(paste0(x,':'),rownames(voiall)))))
voiall <- voiall[roworder,]
miall <- do.call(rbind,ilistmi)
miall <- miall[roworder,]

saveRDS(voiall,paste0('results/voi.Rds'))
saveRDS(miall,paste0('results/mi.Rds'))

library(GGally)
inp3 <- strategies[1];
income_level <- income_levels[1];
results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'.csv'),header=T);
colnames(results)
results$outcome <- results$Cost/results$GDP
results$high <- cut(results$outcome,breaks=c(0,1,5,30))
x11(); ggpairs(results,columns=c(1:17,41),aes(colour=high))
x11(); ggpairs(results,columns=c(18:36,41),aes(colour=high))

## income levels together ####################################

multisource <- list('mean IHR, beta'=c('Mean_IHR', 'beta'),
                    'mean IFR, beta'=c('Mean_IFR', 'beta'),
                    'IHR, IFR, beta'=c('Mean_IHR','Mean_IFR', 'beta'),
                    'Duration + IHR + IFR + beta'=c('beta','Mean_IHR','Mean_IFR','Time_to_discharge'))
listout <- foreach (ks = 1:length(strategies))%dopar%{
  klistvoi <- list()
  klistmi <- list()
  resultslist <- list()
  for (il in 1:length(income_levels)){
    inp3 <- strategies[ks];
    income_level <- income_levels[il];
    resultslist[[il]] <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'.csv'),header=T);
    print(paste0('results/VOI_',inp3,'_',income_level,'.csv'))
  }
  results <- do.call(rbind,resultslist)
  sourcelist <- list()
  for(src in 1:length(multisource)) sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
  
  firstreultcol <- which(colnames(results)=='Cost')
  
  colnames(results) <- gsub('_',' ',colnames(results))
  colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
  colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
  colnames(results)[colnames(results)=='Labour share'] <- 'Labour share of GVA'
  
  sourcemat <- results[,1:(firstreultcol-1)]
  
  outcomes <- results[,firstreultcol:(ncol(results))]
  for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/sourcemat$GDP
  # for(i in which(colnames(outcomes)%in%c('Cost','Deaths')))outcomes[,i] <- log(outcomes[,i])
  
  voilist <- list()
  milist <- list()
  bestfits <- list()
  for(i in 1:ncol(outcomes)){
    voi <- c()
    mi <- c()
    y <- outcomes[,i]
    vary <- var(y) 
    save_fits <- list()
    for(j in 1:ncol(sourcemat)){
      # model outcome as a function of input(s)
      sourcesj <- sourcemat[,j]
      max_degree <- ifelse(is.vector(sourcesj),1,ncol(sourcesj))
      model <- earth::earth(y ~ sourcesj, degree=min(4,max_degree))
      # compute evppi as percentage
      voi[j] <- (vary - mean((y - model$fitted) ^ 2)) / vary * 100
      mi[j] <- infotheo::mutinformation(infotheo::discretize(sourcesj),infotheo::discretize(y))
      save_fits[[j]] <- model$fitted
    }
    for(j in 1:length(sourcelist)){
      sourcesj <- sourcelist[[j]]
      fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj),method='gam')
      mi[j+ncol(sourcemat)] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
      # voi[j+ncol(sourcemat)] <- voi::evppivar(y,sourcesj,pars=colnames(sourcesj))[2]/vary*100
      voi[j+ncol(sourcemat)] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
      save_fits[[j+ncol(sourcemat)]] <- fittedvalues
    }
    milist[[i]] <- mi
    voilist[[i]] <- voi
    whichfit <- which.max(voi)
    bestfits[[i]] <- save_fits[[whichfit]]
  }
  
  bestfits <- do.call(cbind,bestfits)
  voitab <- do.call(rbind,voilist)
  colnames(voitab) <- c(colnames(sourcemat),names(multisource))
  rownames(voitab) <- paste0(colnames(outcomes),', ',inp3)
  mitab <- do.call(rbind,milist)
  colnames(mitab) <- c(colnames(sourcemat),names(multisource))
  rownames(mitab) <- paste0(colnames(outcomes),', ',inp3)
  
  colnames(bestfits) <- paste0(colnames(outcomes),'_pred')
  
  returnresults <- cbind(sourcemat,outcomes,bestfits)
  
  list(voitab, mitab, returnresults)
}

## rerun ###############################################

multisource <- list('Sectors'=c('Agriculture','Food sector'),
                    'Tourism'=c('International tourism','Food sector'),
                    'Education factors'=c('School age','Internet','Labour share of GVA','Unemployment rate'),
                    'Contacts'=c('Work contacts','School contacts','Community contacts','Hospitality contacts'),
                    'Age groups'=c('School age','Working age','Elders'),
                    'Working'=c('Working age','Work contacts'),
                    'Testing'=c('Test rate','Test start'),
                    'Health system'=c('Test rate','Test start','Hospital capacity'),
                    'Social distancing'=c('Social distancing rate','Social distancing max'))
listout2 <- foreach (ks = 1:length(strategies))%dopar%{
  klistvoi <- list()
  klistmi <- list()
  inp3 <- strategies[ks];
  
  results <- listout[[ks]][[3]]
  sourcelist <- list()
  for(src in 1:length(multisource)) sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
  
  firstreultcol <- which(colnames(results)=='Cost')
  
  sourcemat <- results[,1:(firstreultcol-1)]
  
  outcomes_raw <- results[,firstreultcol + 0:3]
  outcome_preds <- results[,(firstreultcol+4):(ncol(results))]
  outcomes <- outcomes_raw - outcome_preds
  
  voilist <- list()
  milist <- list()
  for(i in 1:ncol(outcomes)){
    voi <- c()
    mi <- c()
    y <- outcomes[,i]
    best_pred <- outcome_preds[,i]
    vary <- var(y) 
    for(j in 1:ncol(sourcemat)){
      # model outcome as a function of input(s)
      sourcesj <- cbind(sourcemat[,j],best_pred)
      colnames(sourcesj) <- c('v1','v2')
      fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj),method='gam')
      # compute evppi as percentage
      voi[j] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
      mi[j] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
    }
    for(j in 1:length(sourcelist)){
      sourcesj <- cbind(sourcelist[[j]],best_pred)
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
  rownames(voitab) <- paste0(colnames(outcomes),', ',inp3)
  mitab <- do.call(rbind,milist)
  colnames(mitab) <- c(colnames(sourcemat),names(multisource))
  rownames(mitab) <- paste0(colnames(outcomes),', ',inp3)
  
  list(voitab, mitab)
}

ilistvoi <- list(listout2[[1]][[1]], listout2[[2]][[1]], listout2[[3]][[1]], listout2[[4]][[1]])
ilistmi <- list(listout2[[1]][[2]], listout2[[2]][[2]], listout2[[3]][[2]], listout2[[4]][[2]])

voiall <- do.call(rbind,ilistvoi)
roworder <- unlist(lapply(c("Cost","Deaths","School","GDP loss"),
                          function(x)which(grepl(paste0(x,','),rownames(voiall)))))
voiall <- voiall[roworder,]
miall <- do.call(rbind,ilistmi)
miall <- miall[roworder,]



# saveRDS(voiall,paste0('results/voi.Rds'))
# saveRDS(miall,paste0('results/mi.Rds'))

library(GGally)
inp3 <- strategies[1];
income_level <- income_levels[1];
results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'.csv'),header=T);
colnames(results)
results$outcome <- results$Cost/results$GDP
results$high <- cut(results$outcome,breaks=c(0,1,5,30))
x11(); ggpairs(results,columns=c(1:14,38),aes(colour=high))
x11(); ggpairs(results,columns=c(15:29,38),aes(colour=high))

## something else ###################################

library(scales)
ivoioutcomelist <- list()
for (il in 1:length(income_levels)){
  voioutcomelist <- list()
  for (ks in 1:length(strategies)){
    inp3 <- strategies[ks];
    income_level <- income_levels[il];
    results <- read.csv(paste0('results/VOI_',inp3,'_',income_level,'.csv'),header=T);
    
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
                          Deaths=mean(Deaths)*mean(density)),by=costcat]
    barplt[,midpoint:=x$mids[costcat],by=costcat]
    
    
    meltbarplt <- reshape2::melt((barplt[,2:5]),value.name = 'value',measure.vars=c('School','Deaths','GDP_loss'))
    
    # meltedres <- reshape2::melt((outcomes[,2:4]),value.name = 'value',measure.vars=c('School','Deaths','GDP_loss'))
    meltbarplt$Income_group <- income_level
    meltbarplt$Strategy <- inp3
    voioutcomelist[[ks]] <- meltbarplt
  }
  ivoioutcomelist[[il]] <- do.call(rbind,voioutcomelist)
}
voioutcomes <- do.call(rbind,ivoioutcomelist)
voioutcomes$Income_group <- factor(voioutcomes$Income_group,levels=c('LLMIC','MIC','HIC'))
voioutcomes$Strategy <- factor(voioutcomes$Strategy,levels=c('No Closures','School Closures','Economic Closures','Elimination'))
voioutcomes$variable <- factor(voioutcomes$variable,levels=c('School','Deaths','GDP_loss'),labels=c('Education','YLLs','GDP'))

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

    

## voi vs mi ####################

set.seed(0)

x1 <- rnorm(1000)
y1 <- rnorm(1000)
model <- earth(y1 ~ x1)
vary <- var(y1) 
(vary - mean((y1 - model$fitted) ^ 2)) / vary * 100
# 0.1
mutinformation(discretize(x1),discretize(y1))
# 0.03428454

x2 <- c(x1,10)
y2 <- c(y1,10)
model <- earth(y2 ~ x2)
vary <- var(y2) 
(vary - mean((y2 - model$fitted) ^ 2)) / vary * 100
# 9.194201
mutinformation(discretize(x2),discretize(y2))
# 0.03926193



x1 <- rnorm(1000)
y <- list()
y[[1]] <- x1
y[[2]] <- rnorm(1000,x1,.1)
y[[3]] <- rnorm(1000,x1,.5)
y[[4]] <- rnorm(1000,x1,1)
y[[5]] <- rnorm(1000,x1,5)

infos <- sapply(y,function(yi){
  model <- earth::earth(yi ~ x1, degree=min(4,max_degree))
  # compute evppi as percentage
  voi <- (var(yi) - mean((yi - model$fitted) ^ 2)) / var(yi) * 100
  mi <- infotheo::mutinformation(infotheo::discretize(x1),infotheo::discretize(yi))
  c(voi,mi,cor(x1,yi))
})
scatter <- as.data.frame(do.call(rbind,lapply(y,function(yi)cbind(x1,yi))))
scatter$example <- rep(1:5,each=1000)
infos <- as.data.frame(t(infos))
infos$example <- 1:5
infos$text <- paste0('VOI = ',round(infos$V1),'\n MI = ',round((infos$V2),2),'\n corr = ',round(infos$V3,2))
ggplot() +
  geom_point(data=scatter,aes(x=x1,y=yi),colour='gainsboro',alpha=.5) +
  geom_text(data=infos,aes(x=0,y=0,label=text),colour='navyblue',size=7.5) +
  facet_grid(~example) +
  theme_bw(base_size = 15)+ 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(x='Input',y='Output') -> p
png("../figures/voimicorr.png", height = 200, width = 900)
print(p)
dev.off()



x1[1] <- 100
for(i in 1:length(y)) y[[i]][1] <- 100
infos <- sapply(y,function(yi){
  model <- earth::earth(yi ~ x1, degree=min(4,max_degree))
  # compute evppi as percentage
  voi <- (var(yi) - mean((yi - model$fitted) ^ 2)) / var(yi) * 100
  mi <- infotheo::mutinformation(infotheo::discretize(x1),infotheo::discretize(yi))
  c(voi,mi,cor(x1,yi))
})
scatter <- as.data.frame(do.call(rbind,lapply(y,function(yi)cbind(x1,yi))))
scatter$example <- rep(1:5,each=1000)
infos <- as.data.frame(t(infos))
infos$example <- 1:5
infos$text <- paste0('VOI = ',round(infos$V1),'\n MI = ',round((infos$V2),2),'\n corr = ',round(infos$V3,2))
ggplot() +
  geom_point(data=scatter,aes(x=x1,y=yi),colour='gainsboro',alpha=1) +
  geom_text(data=infos,aes(x=50,y=50,label=text),colour='navyblue',size=7.5) +
  facet_grid(~example) +
  theme_bw(base_size = 15)+ 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(x='Input',y='Output') -> p
png("../figures/voimicorr100.png", height = 200, width = 900)
print(p)
dev.off()

x1[1] <- 100
for(i in 1:length(y)) y[[i]][1] <- 0
infos <- sapply(y,function(yi){
  model <- earth::earth(yi ~ x1, degree=min(4,max_degree))
  # compute evppi as percentage
  voi <- (var(yi) - mean((yi - model$fitted) ^ 2)) / var(yi) * 100
  mi <- infotheo::mutinformation(infotheo::discretize(x1),infotheo::discretize(yi))
  c(voi,mi,cor(x1,yi))
})
scatter <- as.data.frame(do.call(rbind,lapply(y,function(yi)cbind(x1,yi))))
scatter$example <- rep(1:5,each=1000)
infos <- as.data.frame(t(infos))
infos$example <- 1:5
infos$text <- paste0('VOI = ',round(infos$V1),'\n MI = ',round((infos$V2),2),'\n corr = ',round(infos$V3,2))
ggplot() +
  geom_point(data=scatter,aes(x=x1,y=yi),colour='gainsboro',alpha=1) +
  geom_text(data=infos,aes(x=50,y=0,label=text),colour='navyblue',size=7.5) +
  facet_grid(~example) +
  theme_bw(base_size = 15)+ 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(x='Input',y='Output') -> p
png("../figures/voimicorr0.png", height = 200, width = 900)
print(p)
dev.off()
