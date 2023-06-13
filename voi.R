library(earth)
library(ggplot2)
library(infotheo)
library(data.table)
library(viridis)

# setwd('~/overflow_dropbox/DAEDALUS/Daedalus-P2-Dashboard/')
setwd('C:/Users/rj411/OneDrive - Imperial College London/p2_drivers')

## voi #####################################

diseases   <- c('Influenza 2009','Influenza 1957','Influenza 1918',
                'Covid Omicron','Covid Wildtype','Covid Delta',
                'SARS')
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

jd <- 5

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
registerDoParallel(cores=4)


foreach (jd = 1:length(diseases)) %do%{
  ilistvoi <- list()
  ilistmi <- list()
  multisource <- list(c('Agriculture','Food_sector'),
                      c('International_tourism','Food_sector'),
                      c('School_age','Elders'),
                      c('Test_rate','Test_start'),
                      c('Social_distancing_rate','Social_distancing_min'),
                      c('R0','Hospital_capacity','Elders'),
                      c('R0','Social_distancing_rate','Test_rate','Test_start','Social_distancing_min'),
                      c('R0','beta'))
  multisource_names <- c('Sectors','Tourism','Age groups','Testing','Social distancing','Hosp, elders, R0','Test, SD, R0','R0, beta')
  for (il in 1:length(income_levels)){
    klistvoi <- list()
    klistmi <- list()
    for (ks in 1:length(strategies)){
      inp2 <- diseases[jd];
      inp3 <- strategies[ks];
      income_level <- income_levels[il];
      results <- read.csv(paste0('results/VOI_',inp2,'_',inp3,'_',income_level,'.csv'),header=T);
      # labsh <- read.csv(paste0('results/labsh_',inp2,'_',inp3,'_',income_level,'.csv'),header=T);
      # colnames(labsh) <- c('Labour share of GVA','School contacts')
      # results <- cbind(labsh,results)
      sourcelist <- list()
      for(src in 1:length(multisource)) sourcelist[[src]] <- results[,colnames(results)%in%multisource[[src]]]
      
      firstreultcol <- which(colnames(results)=='Cost')
      
      colnames(results) <- gsub('_',' ',colnames(results))
      colnames(results)[colnames(results)=='School age'] <- 'School-age'
      colnames(results)[colnames(results)=='Social distancing min'] <- 'Social distancing max'
      colnames(results)[colnames(results)=='workp'] <- 'Work contacts'
      
      sourcemat <- results[,1:(firstreultcol-1)]
      
      outcomes <- results[,firstreultcol:(ncol(results))]
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
          fittedvalues <- evppifit(y,sourcesj,pars=colnames(sourcesj))
          mi[j+ncol(sourcemat)] <- infotheo::mutinformation(infotheo::discretize(fittedvalues),infotheo::discretize(y))
          # voi[j+ncol(sourcemat)] <- voi::evppivar(y,sourcesj,pars=colnames(sourcesj))[2]/vary*100
          voi[j+ncol(sourcemat)] <- (vary - mean((y - fittedvalues) ^ 2)) / vary * 100
        }
        milist[[i]] <- mi
        voilist[[i]] <- voi
      }
      
      voitab <- do.call(rbind,voilist)
      colnames(voitab) <- c(colnames(sourcemat),multisource_names)
      rownames(voitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3)
      mitab <- do.call(rbind,milist)
      colnames(mitab) <- c(colnames(sourcemat),multisource_names)
      rownames(mitab) <- paste0(colnames(outcomes),': ',income_level,', ',inp3)
      
      klistvoi[[ks]] <- voitab
      klistmi[[ks]] <- mitab
    }
    ilistvoi[[il]] <- do.call(rbind,klistvoi)
    ilistmi[[il]] <- do.call(rbind,klistmi)
  }
  
  voiall <- do.call(rbind,ilistvoi)
  roworder <- unlist(lapply(colnames(outcomes),function(x)which(grepl(paste0(x,':'),rownames(voiall)))))
  voiall <- voiall[roworder,]
  miall <- do.call(rbind,ilistmi)
  miall <- miall[roworder,]
  
  saveRDS(voiall,paste0('results/voi_',inp2,'.Rds'))
  saveRDS(miall,paste0('results/mi_',inp2,'.Rds'))
}

library(scales)
for (jd in 1:length(diseases)) {
  ivoioutcomelist <- list()
  for (il in 1:length(income_levels)){
    voioutcomelist <- list()
    for (ks in 1:length(strategies)){
      inp2 <- diseases[jd];
      inp3 <- strategies[ks];
      income_level <- income_levels[il];
      results <- read.csv(paste0('results/VOI_',inp2,'_',inp3,'_',income_level,'.csv'),header=T);
      
      firstreultcol <- which(colnames(results)=='Cost')
      outcomes <- results[,firstreultcol:(ncol(results))]

      # meltedres <- melt(outcomes[,2:4]/outcomes[,1]*100,value.name = 'value')
      for(i in 1:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/results$GDP*100
      for(i in 2:ncol(outcomes)) outcomes[,i] <- outcomes[,i]/outcomes[,1]
      
      
      setDT(outcomes)
      outcomes[,log10cost:=log10(Cost)]
      x <- hist(outcomes$log10cost,plot = F, breaks=seq(-10,10,by=.15))
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
    annotate(geom = "rect",xmin = 10,xmax = 50,ymin = -Inf,ymax = +Inf,alpha = 0.2) +
    annotate(geom = "rect",xmin = 50,xmax = 100,ymin = -Inf,ymax = +Inf,alpha = 0.4) +
    annotate(geom = "rect",xmin = 100,xmax = 500,ymin = -Inf,ymax = +Inf,alpha = 0.6) +
    annotate(geom = "rect",xmin = 500,xmax = 1000,ymin = -Inf,ymax = +Inf,alpha = 0.8) +
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
  ggsave(p,filename = paste0('logcostshare_',inp2,'.pdf'),width=10,height=10)
}
    
    

nscen <- nrow(voiall)
nparam <- ncol(voiall)
brks <- c(-1,1,10,30,50,70,101)
# brks <- c(-1,1,5,10,20,40,80)/100
discxs <- matrix(cut(as.numeric(voiall),breaks=brks,labels=1:(length(brks)-1)),
                 nrow=nscen,ncol=nparam)
xvals <- 1:nparam
xvals[xvals>17] <- xvals[xvals>17] + .1
xvals[xvals>6] <- xvals[xvals>6] + .1
xvals[xvals>4] <- xvals[xvals>4] + .1
yvals <- 1:nscen
yvals[yvals>12] <- yvals[yvals>12] + .1*(ceiling(yvals[yvals>12]/12)-1)
df <- data.frame(x=rep(xvals,each=nscen),y=rep(yvals,nparam),xs=c(discxs[nscen:1,]))

options(repr.plot.width = 1, repr.plot.height = 0.75)
lbls <- c('<1','1-10','10-30','30-50','50-70','70-100')
p <- ggplot(df, aes(x, y, fill = factor(xs,levels=1:(length(brks)-1),labels=lbls))) +
  geom_tile() +
  scale_fill_manual(values = colorRampPalette(c("midnightblue","royalblue","lightskyblue","white"))(length(brks)-1)) +   
  theme_bw(base_size=12) +                                   # Minimal theme
  labs(title = "") +
  scale_x_continuous(name='',breaks=xvals,labels=colnames(voiall),expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1,vjust=1)) +
  scale_y_continuous(name='',breaks=rev(yvals),labels=rownames(voiall),expand=c(0,0))+
  theme(
    plot.title = element_text(hjust = 1),        
    legend.position="right", legend.spacing.x = unit(0.15, 'cm'),legend.box.spacing = unit(0.0, 'cm')) +               
  guides(fill = guide_legend(#nrow = 1,
    title.theme = element_text(size=10),title.vjust = .8,title.hjust = 1,#title.position = "left",      
    #label.position="bottom", 
    label.hjust = 0, #label.vjust = 0.3, 
    #label.theme = element_text(size=8,angle = 45),
    title = 'EVPPI, % ',override.aes=list(colour='black')))


ggsave(p,filename='evppi.png',width=7,height=4)


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


