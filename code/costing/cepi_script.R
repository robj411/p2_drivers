
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

rawpar <- get_parameters(nsamples=10000, nweeks=5*52)[[1]]



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
  ssv_py_inc_booster <- annualise_ssv_procurement(supplies)
  
  ## ssv allocation 
  allocation_functions = list(function_list, rep_res)[[propflag]]
  alloc_and_del = allocate_and_deliver_doses(allocation_functions, supplies, del_rates=rates,
                                             time_to_approval=time_to_approval)
  alloc = alloc_and_del[['alloc']]
  second_dose_delivery = alloc_and_del[['second_dose_delivery']]
  all_dose_delivery <<- alloc_and_del[['all_dose_delivery']]
  second_dose_delivery$il = factor(second_dose_delivery$il, levels=INCOMELEVELS)
  second_dose_delivery = reshape2::dcast(second_dose_delivery,formula=Week~il,value.var='doses')
  
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
    timescale = timescale
    pos = POS[i,]
    pto = PTO[i,]
    ex = EX[i,]
    icost_lic = (1+pardf$inflation[i])*pardf$cost_lic[i]
    n_ssv_candidates = pardf$n_ssv_candidates[i]
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
                                             ex = ex, 
                                             timescale = timescale, 
                                             icost_lic = icost_lic, 
                                             n_ssv_candidates = n_ssv_candidates)
  
    # capacity reservations
    capres_costs_per_year = get_cap_res_cost(cap_cost_per_dose = cap_cost_per_dose, 
                                             cap_res_vol = cr)
    
    # procurement costs
    ssv_procurement_costs = get_ssv_procurement_costs(ssv_py_inc_booster = ssv_py_inc_booster, 
                                                      cost_res = cost_res, 
                                                      cost_un = cost_un, 
                                                      cap_res_vol = cr)
    ssv_proccost_discounted = sum(ssv_procurement_costs/(1+pardf$discount[i])^c(1:NYEARS))
    ssv_proccost_undiscounted = sum(ssv_procurement_costs)
  
    # apply costs per dose to quantiles
    ssv_delivery_costs_list = get_ssv_delivery_costs(all_ssv_delivered = all_dose_delivery, 
                                                     cost_per_dose = cost_per_dose,
                                                     discount = discount)
  
    # BPSV costs
    if(bpsvflag==T){
      
      inex = INEX[i,]
      inex_weight = pardf$inex_weight[i]
      n_bpsv_candidates = pardf$n_bpsv_candidates[i]
      cost_cogs = pardf$cost_cogs[i]
      cost_ff = pardf$cost_ff[i]
      cost_travel = pardf$cost_travel[i]
      profit = pardf$profit[i]
      cost_bpsv_res = pardf$cost_bpsv_res[i]
      bpsv_replenishment = pardf$bpsv_replenishment[i]
      
      bpsvcosts = get_bpsv_costs(total_bpsv = total_bpsv,
                                 pos = pos, 
                                 pto = pto, 
                                 ex = ex, 
                                 inex = inex, 
                                 inex_weight = inex_weight,
                                 n_bpsv_candidates = n_bpsv_candidates,
                                 durations = durations[1:3]/52, # weeks
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
      inv_cost_per_year = bpsvcosts[['inv_cost_per_year']]*1e3
      bpsvproc = bpsvcosts[['bpsvproc']]
      
      bpsv_del_cost = get_bpsv_del_costs(bpsv_dose_delivery, cost_per_dose)
      # colnames(bpsv_dose_delivery) = INCOMELEVELS
    }else{
      bpsv_rd_costsamples_no_d <- bpsv_rd_costsamples_dyd <- inv_cost_per_year <- bpsvresrd <- bpsvproc <- bpsv_del_cost <- 0
    }
    outvec =     c(enabling_costs, bpsv_rd_costsamples_no_d,bpsv_rd_costsamples_dyd,
      capres_costs_per_year, inv_cost_per_year,
      ssv_rd_costsamples,ssv_proccost_discounted,ssv_proccost_undiscounted,
      ssv_delivery_costs_list,bpsvresrd,bpsvproc,bpsv_del_cost)
    outvec
  } -> x
  outcosts = list()
  for(i in 1:ncol(x)) outcosts[[i]] = unname(x[,i])
  names(outcosts) = c('enabling_costs', 'bpsv_rd_costsamples_no_d','bpsv_rd_costsamples_dyd',
                      'capres_costs_per_year', 'inv_cost_per_year',
                      'ssv_rd_costsamples','ssv_proccost_discounted','ssv_proccost_undiscounted',
                      'ssv_delivery_costs_list','bpsvresrd','bpsvproc','bpsv_del_cost')
  for(i in names(outcosts)) assign(i, outcosts[[i]])
  
  ## return
  return(list(costs=list(upfront=list(enabling=enabling_costs,
                                      bpsv_rd_undiscounted=bpsv_rd_costsamples_no_d,
                                      bpsv_rd_discounted=bpsv_rd_costsamples_dyd),
                         annual=list(capres=capres_costs_per_year,
                                     investigational_reserve=inv_cost_per_year),
                         response=list(ssv_rd=ssv_rd_costsamples,
                                      ssv_proc_discounted=ssv_proccost_discounted,
                                      ssv_proc_undiscounted=ssv_proccost_undiscounted,
                                      ssv_delivery=ssv_delivery_costs_list,
                                      bpsv_response_rd=bpsvresrd,
                                      bpsv_proc=bpsvproc,
                                      bpsv_delivery=bpsv_del_cost)),
              delivery=list(bpsv=combine_llmic(POPS0[4], POPS0[3], bpsv_dose_delivery),
                            ssv=combine_llmic(POPS0[4], POPS0[3], second_dose_delivery))))
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
for(s in 2:4){#1:nscen){
  
  dm = scenario_df$DM[s]
  whichdm <- which(dms==dm)
  pardf$time_to_rd = dm_rd_times[,whichdm]
  bpsvflag = scenario_df$BPSV[s]
  cr <- scenario_df$CR[s]
  propflag = scenario_df$Proportional[s]
  
  # profvis::profvis({
    scenario_results[[s]] = get_scenario(dm, bpsvflag, cr, propflag, rates=scenrates[s,])
    # })
  
  # print(scennames[s])
  # print(summary( scenario_results[[s]]$costs$upfront$bpsv_rd_discounted))
  # for(j in 1:length(scenario_results[[s]]$costs)){
  #   for(k in 1:length(scenario_results[[s]]$costs[[j]])){
  #     print(names(scenario_results[[s]]$costs[[j]])[k])
  #     print(format_to_print2(mean(scenario_results[[s]]$costs[[j]][[k]]),3))
  #   }
  # }
    
}



    
      


      
