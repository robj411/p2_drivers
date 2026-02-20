
args <- commandArgs(TRUE)
# print(args)
if(length(args)==0){
  lb_file = '../data/20251104 updated scenario delivery and costing.xlsx'
}else{
  lb_file = args[1]
}

## start ################################

library(earth)
library(ggplot2)
library(infotheo)
library(data.table)
library(viridis)
library(GGally)
library(splines)
library(doParallel)
library(MASS)
library(openintro)
library(evmix)
library(scales)
library(patchwork)
library(xlsx)
library(rJava)
library(prismatic)
library(ggpubr)
library(geomtextpath)

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java \\jre1.8.0_451/') # for 64-bit version
Sys.setenv(JAVA_HOME='C:\\Program Files (x86)\\Java\\jre1.8.0_451/') # for 32-bit version

cl <- makeCluster(4)
registerDoParallel(cl)

scenario_tab <- readxl::read_xlsx(lb_file, sheet = 'Delivery')
lbsheets = readxl::excel_sheets(lb_file)
# onecost = readxl::read_xlsx(lb_file,sheet = "Lump Sum")
# acosts = readxl::read_xlsx(lb_file,sheet = "Annual Cost")
# rcosts = readxl::read_xlsx(lb_file,sheet = "Response Cost")
strategies <- c('No Closures','School Closures','Economic Closures','Elimination')
income_levels <- c('LLMIC','UMIC','HIC')

result_cols <- c('Cost','YLL','School','GDP_loss')


cost_samples = readRDS('results/cost_samples.Rds')
cost_diffs = list()
for (i in 1:3){
  cost_diffs[[i]] = cost_samples[[i]][,-1]
  for(j in 1:ncol(cost_diffs[[i]]))
    cost_diffs[[i]][,j] = cost_diffs[[i]][,j] - cost_samples[[i]][,1]
}

scenario_names <- colnames(cost_samples[[1]])
nScen = length(scenario_names)
scenario_levels <- 1:nScen
bau_scens = which(grepl('BAU',scenario_names))
scen_scens = which(!grepl('BAU',scenario_names))
bau_names = scenario_names[bau_scens]
just_scen_names = scenario_names[scen_scens]

nbscens = length(bau_names)
ncscens = length(scen_scens)


## plot scenarios ##################################


scen_to_keep <- 1:nScen
scenario_tab[1,c(is.na(scenario_tab[1,]))] = ''

df_list <- list()
for(scen in 1:nScen){
  for(il in 1:3){
    start_index <- (scen-1)*6 + 1 + 2
    get_index <- start_index - 1 + il
    ilevel <- gsub('s','',scenario_tab[2,get_index])
    scenario <- colnames(scenario_tab)[start_index]
    vals <- 0
    for(i in c(0,3)){
      vaccine <- scenario_tab[1,start_index+i][[1]]
      vindex <- get_index + i
      vals <- vals +  as.numeric(scenario_tab[-c(1:2),vindex][[1]])
    }
    df_list[[length(df_list)+1]] <- data.frame(Day=1:length(vals), Percent=cumsum(as.numeric(vals))*100, #Vaccine=vaccine,
                                               Scenario=scenario, Income=ilevel)
  }
}

plotscens <- do.call(rbind,df_list)
ggplot(subset(plotscens,Scenario%in%scenario_names)) + 
  geom_line(aes(x=Day,y=Percent,colour=Scenario),linewidth=1.25,alpha=.75) +
  facet_grid(~factor(Income,levels=c('LLMIC','UMIC','HIC')),scales='free_y') + 
  theme_bw(base_size=15) + labs(colour='') +
  theme(legend.position = 'top') + 
  guides(colour=guide_legend(nrow=2,byrow=F))



## savings and choices ###############################

# read in all scenarios
vsls <- gdps <- topresults <- reslist <- list()
choicestab <- c()
# income_levels <- 'HIC'
for(sl in scenario_levels){
  topresults[[sl]] <- list()
  scenario <- scenario_names[sl]
  allresults <- data.frame()
  for (i in 1:length(income_levels)){
    income_level <- income_levels[i];
    inputtab <- read.csv(paste0('results/inputs_',income_level,'.csv'),header=T);
    vsls[[i]] <- inputtab$vsl * inputtab$gdp
    gdps[[i]] <- 1e6 * inputtab$gdp
    for (k in 1:length(strategies)){
      inp3 <- strategies[k];
      results <- read.csv(paste0('results/outputs/outputs_',inp3,'_',income_level,'_scen',sl,'.csv'),header=T);
      results$policy <- inp3
      results$igroup <- income_level
      results$samplei <- 1:nrow(results)
      results$Costsl <- results$Cost/(inputtab$vsl * inputtab$gdp)
      results$Costpc <- 100 * results$Cost/inputtab$gdp
      results$gdp <- inputtab$gdp
      results$vsl <- inputtab$vsl
      # results <- cbind(results,inputtab)
      allresults <- rbind(allresults,results)
    }
  }
  setDT(allresults)
  allresults[,Costrnd:=Cost]
  allresults[policy!='No Closures',Costrnd:=Cost*(1+runif(sum(allresults$policy!='No Closures'))*1e-8)]
  allresults[,mincost:=Costrnd==min(Costrnd),by=.(igroup,samplei)]
  choices <- allresults[,.(N=sum(mincost)),by=.(igroup,policy)]
  choices$scenario <- scenario
  choicestab <- rbind(choicestab,choices)
  allresults[mincost==T,mean(Cost/gdp),by=.(igroup)]
  
  evalues <- allresults[,mean(Cost/gdp),by=.(igroup,policy)]
  topchoices <- evalues[,list(policy[which.min(V1)],round(V1[which.min(V1)]*100)),by=igroup]
  print(topchoices)
  print(allresults[policy=='No Closures',sum(Breach_before>0),by=igroup])
  # allresults[,keeprow:=policy==topchoices$V1[which(topchoices$igroup==igroup)],by=igroup]
  
  ##!! decision under uncertainty, or decision under certainty?
  # topresults[[bl]][[v]] <- subset(allresults,keeprow)
  topresults[[sl]] <- subset(allresults,mincost)
  reslist[[sl]] <- allresults
  setorder(topresults[[sl]],igroup,samplei)
}
saveRDS(topresults,'results/topresults.Rds')
saveRDS(reslist,'results/reslist.Rds')
saveRDS(list(vsls=vsls,gdps=gdps),'results/vsl_gdp.Rds')

alldiffs <- difftabs <- list()
# get all differences with bau
for(sl in scen_scens){
  alldiffs[[length(alldiffs)+1]] <- list()
  for(refsl in bau_scens){
    
    difftab <- copy(topresults[[refsl]])
    # difftab[,policy:=NULL]
    difftab[,scenCostsl:=topresults[[sl]]$Costsl]
    difftab[,Cost:=Cost-topresults[[sl]]$Cost]
    difftab[,Costorder:=Costpc]
    difftab[,Costpc:=Costpc-topresults[[sl]]$Costpc]
    difftab[,Costsl:=Costsl-topresults[[sl]]$Costsl]
    difftab[,YLL:=YLL-topresults[[sl]]$YLL]
    difftab[,School:=School-topresults[[sl]]$School]
    difftab[,GDP_loss:=GDP_loss-topresults[[sl]]$GDP_loss]
    difftab[,scen_Exit_wave:=topresults[[sl]]$Exit_wave]
    difftab[,scen_Mitigated_deaths:=topresults[[sl]]$Mitigated_deaths]
    difftab[,scen_deaths:=topresults[[sl]]$Deaths]
    difftab[,scen_deaths1:=topresults[[sl]]$Deaths1]
    difftab[,scen_deaths2:=topresults[[sl]]$Deaths2]
    difftab[,scen_deaths3:=topresults[[sl]]$Deaths3]
    difftab[,scen_deaths4:=topresults[[sl]]$Deaths4]
    difftab[,scen_policy:=topresults[[sl]]$policy]
    difftab[,scen_End_simulation:=topresults[[sl]]$End_simulation]
    
    # difflist[[length(difflist)+1]] <- difftab
    
    saveRDS(difftab,paste0('results/differences/difftab_',refsl,'_',sl,'.Rds'))
    if(refsl==1) difftabs[[length(difftabs)+1]] <- difftab
    
    alldiffs[[length(alldiffs)]][[which(bau_scens==refsl)]] <- difftab
  }
}
saveRDS(alldiffs,'results/alldiffs.Rds')


# summary(subset(reslist[[sl]],policy=='No Closures'&R0<2.6&R0>2.4)$Deaths/50000)

saveRDS(choicestab,paste0('results/choicestab.Rds'))

with(do.call(rbind,topresults),table(policy,igroup))

plotdur <- list()
for(i in scenario_levels){
  plotdur[[i]] <- copy(topresults[[i]])
  plotdur[[i]]$scenario <- scenario_names[i]
}
plotdur <- do.call(rbind,plotdur)
mplotdur <- melt(plotdur[,.(scenario,igroup,End_mitigation,End_simulation)],measure.vars=c('End_mitigation','End_simulation'))
ggplot(mplotdur) + 
  geom_density(aes(x=value/365,y=..density..),fill='darkorange',colour='midnightblue',alpha=.75,linewidth=2) +
  facet_grid(~factor(variable,labels=c('Mitigation end','Simulation end')),scales='free') +
  theme_bw(base_size=15) +
  scale_x_continuous(limits=c(0,NA)) +
  scale_y_continuous(expand=c(0,NA)) +
  labs(x='Day',y='Density')

dispcols <- c('Cost','GDP_loss','YLL','School',
              'policy','samplei','Deaths','igroup')
unique(subset(topresults[[1]],End_simulation>3600))[,..dispcols]

endhosp <- list()
for(j in scenario_levels){
  endhosp[[j]] <- list()
  endhosp[[j]] <- copy(topresults[[j]])
  endhosp[[j]]$scenario <- scenario_names[i]
}
endhosp <- do.call(rbind,endhosp)

# unique(subset(endhosp,End_hosp>1000)[,.(igroup,policy,samplei)])

## negatives ###########################

dispcols <- c('Cost','GDP_loss','YLL','School',
              'policy','scen_policy','samplei','Deaths4','scen_deaths4','Deaths','scen_deaths','displaced')
storecols <- c('Cost','GDP_loss','YLL','School',#'yll_per_death1','yll_per_death2','yll_per_death3','yll_per_death4',
               'Deaths1','scen_deaths1','Deaths2','scen_deaths2','Deaths3','scen_deaths3','Deaths4','scen_deaths4','Deaths','scen_deaths','policy','scen_policy','igroup','scenario','gdp','samplei')
nneg <- 0
ncdisp <- exitwaves <- others <- drag <- c()
for(sl in 1:ncscens){
  difftab <- difftabs[[sl]]
  for(income_level in income_levels){
    subtab <- subset(difftab,igroup==income_level)
    subtab$sample <- 1:nrow(subtab)
    subtab$scenario <- sl
    subtab[,drags:=End_simulation<scen_End_simulation]
    subtab[,noclosuredisplacement:=Deaths4>scen_deaths4 & policy==scen_policy]
    subtab[,exitwave:=scen_Mitigated_deaths < Mitigated_deaths & Deaths < scen_deaths]
    drag <- rbind(drag,subset(subtab,Cost < 0 & drags &! noclosuredisplacement)[,..storecols])
    exitwaves <- rbind(exitwaves,subset(subtab,Cost < 0 & exitwave & !drags &!noclosuredisplacement)[,..storecols])
    ncdisp <- rbind(ncdisp, subset(subtab,Cost < 0 & noclosuredisplacement )[,..storecols])
    others <- rbind(others, subset(subtab,Cost < 0 & !noclosuredisplacement & !exitwave& !drags)[,..storecols])
    subtab[,displaced:=Deaths4>scen_deaths4 & Deaths < scen_deaths]
    befores <- subset(subtab,Cost < 0 & scen_Mitigated_deaths <= Mitigated_deaths &!noclosuredisplacement&!exitwave)[,..dispcols]
    if(nrow(befores)>0){
      print(c('before',sl,income_level))
      print(befores)
    }
    afters <- subset(subtab,Cost < 0 & scen_Mitigated_deaths > Mitigated_deaths &!noclosuredisplacement&!exitwave)[,..dispcols]
    if(nrow(afters)>0){
      print(c('after',sl,income_level))
      print(afters)
    }
    nneg <- nneg + nrow(subset(subtab,Cost< 0 ))
  }
} 
# nrow(drag)
# summary(with(drag,Cost/gdp*100))
# nrow(exitwaves)
# summary(with(exitwaves,Cost/gdp*100))
# nrow(ncdisp)
# summary(with(ncdisp,Cost/gdp*100))
# nrow(others)
# summary(with(others,Cost/gdp*100))

for (or in 1:nrow(others)){
  si <- others[or,samplei]
  scen <- others[or,scenario]
  ig <- others[or,igroup]
  pol <- others[or,policy]
  subset(reslist[[scen]],igroup==ig&policy==pol&samplei==si)[,c('Cost','GDP_loss','YLL','School','Deaths1','Deaths2','Deaths3','Deaths4','Deaths')]
  subset(reslist[[1]],igroup==ig&policy==pol&samplei==si)[,c('Cost','GDP_loss','YLL','School','Deaths1','Deaths2','Deaths3','Deaths4','Deaths')]
}

# 
# agevals <- c('0-4','5-19','20-64','65+','Total')
# barpt <- setDT(reshape2::melt(subset(ncdisp,samplei==1859&scenario==4)[,grepl('death',colnames(ncdisp),ignore.case = T),with=F]))
# barpt$var <- 'BAU'
# barpt$var[grepl('yll',barpt$variable)] <- 'YLLperdeath'
# barpt$var[grepl('scen',barpt$variable)] <- 'Scenario'
# barpt[,agegroup:=substr(variable,nchar(as.character(variable)),nchar(as.character(variable))),by=variable]
# barpt[,ages:=agevals[as.numeric(agegroup)],by=agegroup]
# barpt$ages <- factor(barpt$ages,levels=agevals,labels=agevals)
# barpt2 <- dcast(barpt,formula=ages~var)
# barpt2[,deathsaverted:=BAU-Scenario]
# barpt2[,YLLaverted:=YLLperdeath*deathsaverted]
# barpt2$ages[is.na(barpt2$ages)] <- 'Total'
# barpt2$YLLaverted[is.na(barpt2$YLLaverted)] <- sum(barpt2$YLLaverted,na.rm=T)
# 
# 
# library(patchwork)
# ylims = range(barpt2$deathsaverted/1000)
# (p1 <- ggplot(barpt2) + geom_bar(aes(x=ages,y=deathsaverted/1000,fill=ages),stat='identity',show.legend = F) +
#     theme_bw(base_size = 15) +
#     scale_fill_manual(values=c(`0-4`='chocolate3',`5-19`='chocolate3',`20-64`='chocolate3',`65+`='chocolate3','Total'='midnightblue')) +
#     labs(x='',y='Deaths averted, thousand') +
#     scale_y_continuous(limits=ylims)
# )
# p2 <- ggplot(barpt2) + geom_bar(aes(x=ages,y=YLLaverted/1e5,fill=ages),stat='identity',show.legend = F) +
#   theme_bw(base_size = 15) +
#   scale_fill_manual(values=c(`0-4`='chocolate3',`5-19`='chocolate3',`20-64`='chocolate3',`65+`='chocolate3','Total'='midnightblue')) +
#   labs(x='',y='YLL averted, hundred thousand') +
#   scale_y_continuous(limits=ylims)
# p1 + p2
# ggsave(p1 + p2,filename='results/displacement.png',width=7)



get_pcs <- function(df, x = "x", y = "y") {
  
  stopifnot(all(c(x, y) %in% names(df)))
  d <- na.omit(df[, c(x, y)])
  names(d) <- c("x", "y")
  
  # PCA to get axes
  pc <- prcomp(d, center = TRUE, scale. = FALSE)
  v1 <- pc$rotation[, 1]          # axis of max variance (unit vector)
  v2 <- pc$rotation[, 2]          # perpendicular axis (unit vector)
  ctr <- pc$center                # mean used for centering
  sc  <- pc$x                     # scores in PC coordinates
  
  # IQR on PC1 (along v1)
  q1_pc1  <- as.numeric(quantile(sc[,1], 0.25))
  q3_pc1  <- as.numeric(quantile(sc[,1], 0.75))
  med_pc1 <- as.numeric(median(sc[,1]))
  
  # IQR on PC2 (along v2)
  q1_pc2  <- as.numeric(quantile(sc[,2], 0.25))
  q3_pc2  <- as.numeric(quantile(sc[,2], 0.75))
  
  # Helper to go from center + a*v1 + b*v2 back to (x, y)
  to_xy <- function(a, b) {
    # (x, y) = ctr + a*v1 + b*v2
    c(ctr[1] + a*v1[1] + b*v2[1],
      ctr[2] + a*v1[2] + b*v2[2])
  }
  
  # Segment along PC1: from Q1 to Q3 (b = 0)
  p1 <- to_xy(q1_pc1, 0)
  p2 <- to_xy(q3_pc1, 0)
  
  # Segment along PC2 through the PC1 median point: a = med_pc1
  p3 <- to_xy(med_pc1, q1_pc2)
  p4 <- to_xy(med_pc1, q3_pc2)
  
  segs <- rbind(
    data.frame(x = p1[1], y = p1[2], xend = p2[1], yend = p2[2],
               axis = "PC1 IQR")#,
    # data.frame(x = p3[1], y = p3[2], xend = p4[1], yend = p4[2],
    #            axis = "PC2 IQR @ PC1 median")
  )
  
  segs
  
  
}
