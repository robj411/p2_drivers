
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


source('functions.R')

library(tidyr)
library(dplyr)
library(stringi)
library(data.table)
library(splines)
library(haven)
library(hrbrthemes)
library(MASS)
library(statmod)
library(TruncExpFam)
library(EnvStats)
library(extraDistr)
library(PearsonDS)
library(foreach)
library(parallel)

# cl <- makeCluster(4)
# registerDoParallel(cl)

set.seed(0)

## levels for scenario options

basic_rates = c(.07, .07, .02, .02)
crs = c(0, 0.7, 2)
dms = c(365, 200, 100)



## read and sample parameters #####################################
# also sets model variables

caps <- readODS::read_ods('cost_parameters.ods',sheet = 1)

rawpar <- get_parameters(nsamples=10000)[[1]]



## scenario variables  ######################################################



## scenario inputs

scennames = c('BAU',paste0('S',sprintf("%02d", 1:12)))
nscen = length(scennames)
pop_proportional = c(1,1,2,2,1,1,1,2,2,1,2,2,2)
scen_dm = c(1, 1,1,1, 1,1, 2,2,2, 3,3,3, 1)
scen_capres = c(1, 1,2,3, 2,3, 1,2,3, 1,2,3, 1)
bpsv_scen = rep(F,length(scennames))
bpsv_scen[2:4] = T

scenrates = matrix(rep(basic_rates, nscen), ,ncol=length(INCOMELEVELS),byrow=T)
scenrates[13,3:4] = 0.04

scenario_df = data.frame(Scenario=scennames, Proportional=pop_proportional, DM=dms[scen_dm], CR=crs[scen_capres], BPSV=bpsv_scen, scenrates)




##################################################################
##################################################################

# functions for allocation
function_list = list(res=allocate_res_doses, ex=allocate_exbui_doses, bui=allocate_exbui_doses)
rep_res = list(res=allocate_res_doses, ex=allocate_res_doses, bui=allocate_res_doses)

# dm=365; bpsvflag=F; cr=0; propflag=1; rates=c(.07, .07, .02)
# dm=365; bpsvflag=T; cr=2; propflag=1; rates=c(.07, .07, .02)
get_scenario = function(dm=365, bpsvflag=F, cr=0, propflag=1, rates=c(.07, .07, .02, .02)){

  # timings of the phases
  phase_duration <<- sapply(0:3,function(x) pfixed[[paste0('weeks_P',x,'_',dm)]])
  
  ##!! should be sum of phases
  time_to_approval = round(dm/7)

  ## delivery ###################################################################
  # delivery variables are not sampled in this implementation because each rollout realisation 
  # is written to csv (for input into another pipeline)

  ## ssv supply
  ssv_supply = get_ssv_supply(dm=dm, capres=cr, bpsv=bpsvflag, weeks_init = WEEKS_INIT, weeks_scale = WEEKS_SCALE)
  supplies = ssv_supply[['supplies']]
  
  ## ssv allocation 
  allocation_functions = list(function_list, rep_res)[[propflag]]
  alloc_and_del = allocate_and_deliver_doses(allocation_functions, supplies, del_rates=rates,
                                             time_to_approval=time_to_approval, 
                                             hic_only = pfixed$hic_cap_res/(cr+pfixed$hic_cap_res))
  alloc = alloc_and_del[['alloc']]
  second_dose_delivery = alloc_and_del[['second_dose_delivery']]
  all_dose_delivery <<- alloc_and_del[['all_dose_delivery']]
  second_dose_delivery$il = factor(second_dose_delivery$il, levels=INCOMELEVELS)
  second_dose_delivery = reshape2::dcast(second_dose_delivery,formula=Week~il,value.var='doses')
  
  ssv_py_inc_booster <- annualise_ssv_procurement(supplies, alloc, cap_res = pfixed$hic_cap_res + cr)
  
  # BPSV
  if(bpsvflag==T){
    bpsv_supply = get_bpsv_supply()
    b_supply = bpsv_supply$b_supply
    b_alloc = bpsv_supply$b_alloc
    bpsv_dose_delivery = bpsv_supply$bpsv_dose_delivery
    total_bpsv = max(b_supply$doses)
  }else{
    bpsv_dose_delivery <- NULL
  }
  
  ## costs ###########################################################################
  # enabling costs
  # enabling_costs <- rep(0,NSAMPLES)
  # capres_costs_per_year = rep(0,NSAMPLES)
  # 
  # ssv_proccost_discounted <- ssv_proccost_undiscounted <- rep(0, NSAMPLES)
  # ssv_rd_costsamples = rep(0,NSAMPLES)
  # ssv_delivery_costs_list = rep(0,NSAMPLES)
  # 
  # bpsv_del_cost <- bpsv_rd_costsamples_dyd <- bpsv_rd_costsamples_no_d <- bpsvresrd <- 
  #   inv_cost_per_year <- inv_cost <- bpsvproc <- rep(0,NSAMPLES)
  
  foreach(i = 1:NSAMPLES, .combine='rbind')%do%{
    time_to_rd = pardf$time_to_rd[i]
    
    # scaling for phase costs
    durations = OLD_DURATIONS[i,]
    timescale = phase_duration/durations #get_duration_scalar(phase_duration,OLD_DURATIONS[i,])
    
    # parameters for functions
    pos = POS_SSV[i,]
    pto = PTO_SSV[i,]
    inflation = pardf$inflation[i]
    ex = EX[i,]
    exi = (1+inflation)*ex
    cost_lic = pardf$cost_lic[i]
    icost_lic = (1+inflation)*cost_lic
    n_ssv_successes = pardf$n_ssv_successes[i]
    q_ssv_successes = pardf$q_ssv_successes[i]
    cost_un = pardf$cost_un[i]
    cost_res = pardf$cost_cogs[i]*(1+pardf$profit[i])*(1+pardf$cost_travel[i])
    discount = pardf$discount[i]
    cost_per_dose = COST_PER_DOSE_LIST[[i]]
    cost_enab = pardf$cost_enab[i]
    cap_cost_per_dose = pardf$cost_capres[i]
    
    # enabling costs
    enabling_costs = get_enabling_costs(cost_enab = cost_enab, 
                                        discount = discount, 
                                        time_to_target_dm = time_to_rd)
    
    # ssv r&d
    ssv_rd_costsamples = get_ssv_randd_costs(pos = pos,
                                             pto = pto, 
                                             exi = exi, 
                                             timescale = timescale, 
                                             cost_lic = cost_lic, 
                                             n_ssv_successes, 
                                             q_ssv_successes)
  
    # capacity reservations
    capres_costs_per_year = get_cap_res_cost(cap_cost_per_dose = cap_cost_per_dose, 
                                             cap_res_vol = cr)
    
    # procurement costs
    ssv_procurement_costs = get_ssv_procurement_costs(ssv_py_inc_booster = ssv_py_inc_booster, 
                                                      cost_res = cost_res, 
                                                      cost_un = cost_un)
    ssv_proccost_discounted = sum(ssv_procurement_costs/(1+pardf$discount[i])^c(1:NYEARS))
    ssv_proccost_undiscounted = sum(ssv_procurement_costs)
  
    # apply costs per dose to quantiles
    ssv_delivery_costs_list = get_ssv_delivery_costs(all_ssv_delivered = all_dose_delivery, 
                                                     cost_per_dose = cost_per_dose,
                                                     discount = discount)
    dis_ssv_delivery_costs = ssv_delivery_costs_list[['dis_countrycosts']]
    ssv_delivery_costs = ssv_delivery_costs_list[['countrycosts']]
    
    # BPSV costs
    if(bpsvflag==T){
      
      pos = POS_BPSV[i,]
      pto = PTO_BPSV[i,]
      inex = INEX[i,]
      inexi = (1+inflation)*inex
      inex_weight = pardf$inex_weight[i]
      n_bpsv_candidates = pardf$n_bpsv_candidates[i]
      cost_cogs = pardf$cost_cogs[i]
      cost_ff = pardf$cost_ff[i]
      cost_travel = pardf$cost_travel[i]
      profit = pardf$profit[i]
      cost_bpsv_res = pardf$cost_bpsv_res[i]
      bpsv_replenishment = pardf$bpsv_replenishment[i]
      n_bpsv_p1 = pardf$n_bpsv_p1[i]
      bpsv_res_upfront = pardf$bpsv_res_upfront[i]
      
      bpsvcosts = get_bpsv_costs(total_bpsv = total_bpsv,
                                 pos = pos, 
                                 pto = pto, 
                                 ex = ex, 
                                 exi= exi,
                                 inexi = inexi, 
                                 inex_weight = inex_weight,
                                 n_bpsv_candidates = n_bpsv_candidates,
                                 n_bpsv_p1 = n_bpsv_p1,
                                 bpsv_res_upfront = bpsv_res_upfront,
                                 y_durations = durations[1:3]/52,
                                 old_duration = durations[4], # years
                                 discount = discount,
                                 icost_lic = icost_lic,
                                 bpsv_replenishment = bpsv_replenishment,
                                 cost_cogs = cost_cogs,
                                 cost_bpsv_res = cost_bpsv_res,
                                 cost_res = cost_res,
                                 profit = profit,
                                 cost_travel = cost_travel,
                                 cost_ff = cost_ff)
      
      bpsv_rd_costsamples_dyd = bpsvcosts[['bpsv_rd_costsamples_dyd']]
      bpsv_rd_costsamples_no_d = bpsvcosts[['bpsv_rd_costsamples_no_d']]
      bpsvresrd = bpsvcosts[['bpsvresrd']]
      inv_cost_per_year = bpsvcosts[['inv_cost_per_year']]
      bpsvproc = bpsvcosts[['bpsvproc']]
      dis_upfront_bpsv = bpsvcosts[['dis_bpsv_upfront']]
      upfront_bpsv = bpsvcosts[['bpsv_upfront']]
      
      bpsv_del_cost = get_bpsv_del_costs(bpsv_dose_delivery, cost_per_dose)
      # colnames(bpsv_dose_delivery) = INCOMELEVELS
    }else{
      bpsv_rd_costsamples_no_d <- bpsv_rd_costsamples_dyd <- inv_cost_per_year <- bpsvresrd <- bpsvproc <- bpsv_del_cost <- upfront_bpsv <- dis_upfront_bpsv <- 0
    }
    outvec =     c(enabling_costs, bpsv_rd_costsamples_no_d,bpsv_rd_costsamples_dyd,
                   capres_costs_per_year, inv_cost_per_year,
                   ssv_rd_costsamples,ssv_proccost_discounted,ssv_proccost_undiscounted,
                   dis_ssv_delivery_costs,ssv_delivery_costs,
                   bpsvresrd,bpsvproc,bpsv_del_cost)
    outvec
  } -> x
  xmat = matrix(x, nrow=NSAMPLES, byrow=F)
  # print(xmat)
  outcosts = list()
  for(i in 1:ncol(xmat)) outcosts[[i]] = unname(xmat[,i])
  names(outcosts) = c('enabling_costs', 'bpsv_rd_costsamples_no_d','bpsv_rd_costsamples_dyd',
                      'capres_costs_per_year', 'inv_cost_per_year',
                      'ssv_rd_costsamples','ssv_proccost_discounted','ssv_proccost_undiscounted',
                      'dis_ssv_delivery_costs','ssv_delivery_costs',
                      'bpsvresrd','bpsvproc','bpsv_del_cost')
  for(i in names(outcosts)) assign(i, outcosts[[i]])
  
  ## return
  return(list(costs=list(upfront=list(enabling=enabling_costs,
                                      upfront_bpsv=upfront_bpsv,
                                      dis_upfront_bpsv=dis_upfront_bpsv,
                                      bpsv_rd_undiscounted=bpsv_rd_costsamples_no_d,
                                      bpsv_rd_discounted=bpsv_rd_costsamples_dyd),
                         annual=list(capres=capres_costs_per_year,
                                     investigational_reserve=inv_cost_per_year),
                         response=list(ssv_rd=ssv_rd_costsamples,
                                      ssv_proc_discounted=ssv_proccost_discounted,
                                      ssv_proc_undiscounted=ssv_proccost_undiscounted,
                                      ssv_delivery_discounted=dis_ssv_delivery_costs,
                                      ssv_delivery_undiscounted=ssv_delivery_costs,
                                      bpsv_response_rd=bpsvresrd,
                                      bpsv_proc=bpsvproc,
                                      bpsv_delivery=bpsv_del_cost)),
              delivery=list(bpsv=combine_llmic(POPS0[4], POPS0[3], bpsv_dose_delivery),
                            ssv=combine_llmic(POPS0[4], POPS0[3], second_dose_delivery),
                            ssv4=second_dose_delivery)))
}
  
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(format(x,scientific = F))), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

format_to_print2 <- function(x,z=2){
  v <- signif(x,z)
  sapply(v, function(i)  formatC(i, format="f", digits=decimalplaces(i), big.mark=","))
}

# CEPI assumptions
# times to R&D for the three missions
dm_rd_times <- cbind(0, pardf$years_200, pardf$years_100)

scenario_results = list()
tosave <- list()
for(i in 1:3) tosave[[i]] <- list()


sheet2 = read.csv('../../cepi_results/Delta_LIR_IQR_pc_GDP_BAU.csv',check.names = F)


for(s in 1:nscen){ #c(1,10)){# 
  
  dm = scenario_df$DM[s]
  whichdm <- which(dms==dm)
  pardf$time_to_rd = dm_rd_times[,whichdm]
  bpsvflag = scenario_df$BPSV[s]
  cr <- scenario_df$CR[s]
  propflag = scenario_df$Proportional[s]
  
    scenario_results[[s]] = get_scenario(dm, bpsvflag, cr, propflag, rates=scenrates[s,])
    # })
  
  thisscen = scenario_results[[s]]
  allupfront = with(thisscen$costs$upfront, enabling + dis_upfront_bpsv + bpsv_rd_discounted)
  allannual = with(thisscen$costs$annual, capres + investigational_reserve)
  allresp = with(thisscen$costs$response, ssv_rd + ssv_proc_discounted + ssv_delivery_discounted + 
                   bpsv_response_rd + bpsv_proc + bpsv_delivery)
  
  if(s==1){
    bauallupfront = allupfront
    bauallannual = allannual
    bauallresp = allresp
    baucosts = list(allupfront, allannual, allresp)
    
    print(summary(with(thisscen$costs$response, ssv_rd + ssv_proc_undiscounted + ssv_delivery_undiscounted + 
                         bpsv_response_rd + bpsv_proc + bpsv_delivery)))
  }else{
    diffallupfront = allupfront - bauallupfront
    diffallannual = allannual - bauallannual
    diffallresp = allresp - bauallresp
    
    delta_lir = sheet2[match(scennames[s],sheet2$to), 3:4]
    
    cat(paste(scennames[s], ' & ', paste0(format_to_print2(quantile(diffallupfront, c(1,3)/4)),collapse='--{}')
                , ' & ', paste0(format_to_print2(quantile(diffallannual, c(1,3)/4)),collapse='--{}')
              , ' & ', paste0(format_to_print2(quantile(diffallresp, c(1,3)/4)),collapse='--{}')
              , ' & ', paste0(delta_lir,collapse='--{}')
              , '\\\\ \n'
    ))
    
    tosave[[1]][[scennames[s]]] <- allupfront
    tosave[[2]][[scennames[s]]] <- allannual
    tosave[[3]][[scennames[s]]] <- allresp
  }
  
}

poptarget = POPS15/POPS0
reachday = matrix(0,nrow=nscen,ncol=4)
for(s in 1:nscen){
  thisscen = scenario_results[[s]]$delivery$ssv4
  for(j in 2:ncol(thisscen)){
    index = which(thisscen[,j]>0.4*poptarget[j-1])[1]
    reachday[s,j-1] = thisscen[index,1]
  }
}
sapply(1:nscen,function(x) cat(paste(paste0(c(scennames[x], reachday[x,]), collapse=' & '),'\\\\\n'))) -> x

cat(paste0('\\renewcommand*{\\MinNumber}{',min(reachday),'} \n\\renewcommand*{\\MaxNumber}{',max(reachday),'}'))


if(NSAMPLES>=10000){
  bounds = lapply(1:3,function(x) cbind(BAU=baucosts[[x]], do.call(cbind, tosave[[x]])))
  saveRDS(bounds,file = '../results/cost_samples.Rds')
}
    
      

scenario_results[[s]]$delivery$ssv -> deliveries


diurnise_deliveries = function(deliveries){
  # dimension of daily rollout
  maxdays = max(deliveries$Week)*7
  # matrix for results
  daily_doses = matrix(0,ncol=ncol(deliveries),nrow=maxdays)
  # one income level at a time
  for(j in 2:ncol(deliveries)){
    # deliveries are cumulative
    doses_per_week = diff(c(0,deliveries[,j]))
    # final week with deliveries
    last_week = length(doses_per_week) - which(rev(doses_per_week)>0)[1] + 1
    # take the recent max
    penultimate_week_doses = max(doses_per_week[last_week-c(1:8)]) 
    for(i in 1:last_week){
      # week to day
      daily_doses[(i-1)*7+1:7,j] = rep(doses_per_week[i],7)/7*100
    }
    # extend indefinitely - allows for different population sizes
    daily_doses[last_week:maxdays, j] = penultimate_week_doses/7
  }
  daily_doses[,1] = 1:maxdays
  colnames(daily_doses) = colnames(deliveries)
  colnames(daily_doses)[1] = 'Day' 
  daily_doses
}

diurnise_deliveries(deliveries)[(69*6):(70*7),]


      
