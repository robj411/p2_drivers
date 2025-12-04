
## setting up ######################################

get_parameters = function(nsamples = 100,
                          nweeks = 52){
  
  ## global variables ############################
  
  NSAMPLES <<- nsamples
  NWEEKS <<- nweeks
  NYEARS <<- ceiling(NWEEKS/52)
  MAN_TYPES <<- c('res','ex','bui')
  PC_THRESHOLDS <<- c(0,11,31)
  INCOMELEVELS <<- c('HIC','UMIC','LMIC','LIC')
  NLEVELS <<- length(INCOMELEVELS)
  
  
  ## parameters 
  
  
  params <- rawpar <- list()
  
  for(i in 1:nrow(caps)){
    parameters <- as.numeric(strsplit(caps$Parameters[i],',')[[1]])
    bounds <- as.numeric(strsplit(caps$Bounds[i],',')[[1]])
    if(caps$Distribution[i]=='Uniform'){
      distcall = function(nsamples) runif(nsamples,parameters[1],parameters[2])
    }else if(caps$Distribution[i]=='Exponential'){
      distcall = function(nsamples) rexp(nsamples, rate = 1/parameters[1])
    }else if(caps$Distribution[i]=='Inverse Gaussian'){ # s = dispersion
      distcall = function(nsamples) rinvgauss(nsamples, mean = parameters[1], shape = parameters[2])
    }else if(caps$Distribution[i]=='PearsonV'){ 
      distcall = function(nsamples) rinvgamma(nsamples, alpha = parameters[1], beta = parameters[2])
      # distcall = function(nsamples) rpearsonV(nsamples, shape = parameters[1], scale = parameters[2], location = 0)
    }else if(caps$Distribution[i]=='Log normal'){ 
      m <- parameters[1]
      s <- parameters[2]
      location <- log(m^2 / sqrt(s^2 + m^2))
      shape <- sqrt(log(1 + (s^2 / m^2)))
      distcall = function(nsamples) rlnorm(nsamples, meanlog = location, sdlog = shape)
    }else if(caps$Distribution[i]=='PearsonVI'){ 
      # distcall = function(nsamples) rbetapr(nsamples, shape1 = parameters[1], shape2 = parameters[2], scale = parameters[3])
      distcall = function(nsamples) rpearsonVI(nsamples, a = parameters[1], b = parameters[2], scale = parameters[3], location=0)
    }else if(caps$Distribution[i]=='Triangular'){
      distcall = function(nsamples) rtri(nsamples,min=parameters[1],mode=parameters[2],max=parameters[3])
    }else if(caps$Distribution[i]=='Multinomial'){
      values = na.omit(c(parameters))
      distcall = function(nsamples) sample(values,size=nsamples,replace=T)
    }
    
    if(caps$Distribution[i]!='Constant'){
      samples = distcall(NSAMPLES)
      if(length(bounds)>1){
        fails = samples < bounds[1] | samples > bounds[2]
        while(sum(fails)>0){
          samples[fails] = distcall(sum(fails))
          fails = samples < bounds[1] | samples > bounds[2]
        }
      }
      params[[caps$Description[i]]] = samples
    }
    
    # params has only the samples; rawpar has also the constants
    if(caps$Distribution[i]=='Constant'){
      rawpar[[caps$Code[i]]] = parameters
    }else{
      rawpar[[caps$Code[i]]] = params[[caps$Description[i]]]
    }
  }
  
  ics = c('hic','umic','lmic','lic')
  
  # some variables need to be fixed atm as they feed into delivery trajectories
  pop0pars = paste0('pop_',ics,'_0')
  pop15pars = paste0('pop_',ics,'_15')
  pop65pars = paste0('pop_',ics,'_65')
  
  dms <<- c(365, 200, 100)
  phase_dur_par = apply(expand.grid(0:3,dms), 1, function(x) paste0('weeks_P',x[1],'_',as.numeric(x[2])))
  weeks_init_pars = paste0('weeks_init_',MAN_TYPES)
  weeks_init_ex_pars = paste0('weeks_init_ex_',c('nb','bp'))
  weeks_scale_pars = paste0('weeks_scale_',MAN_TYPES)
  fixedpars = c(pop0pars, pop15pars, pop65pars,'final_vaccine_coverage','vaccine_wastage', # volumes made
                phase_dur_par,'duration_3_resp', # phase durations
                weeks_init_pars,weeks_init_ex_pars,weeks_scale_pars,'week_trans_start', # weeks to manufacturing
                'man_curr','man_glo', # manufacturing capacities
                'bpsv_inv_res','hic_cap_res') # reserve amounts
  pfixed = list()
  for(p in fixedpars) 
    pfixed[[p]] = rawpar[[p]]
  pfixed <<- pfixed
  
  pardf <<- as.data.frame(do.call(cbind,rawpar[!names(rawpar)%in%names(pfixed)]))
  
  ## populations and vaccine demand
  
  pops0 = sapply(pop0pars,function(x)pfixed[[x]])
  # pops0 = pops[1:3]
  # pops0[3] = sum(pops[3:4])
  POPS0 <<- pops0
  
  pops3 = sapply(pop15pars,function(x)pfixed[[x]])
  # pops3 = pops[1:3]
  # pops3[3] = sum(pops[3:4])
  POPS15 <<- pops3
  
  pops65 = sapply(pop65pars,function(x)pfixed[[x]])
  # pops365 = pops65[1:3]
  # pops365[3] = sum(pops65[3:4])
  POPS65 <<- pops65
  
  DEMAND65 <<- pfixed$final_vaccine_coverage*POPS65 / (1-pfixed$vaccine_wastage) # total bpsv doses per income level
  POPFRACS65 <<- POPS65/sum(POPS65)
  TOTAL_BPSV_DEMAND <<- 1 # billion. not: sum(DEMAND65/1e9)
  
  DEMAND15 <<- pfixed$final_vaccine_coverage*POPS15*2 / (1-pfixed$vaccine_wastage)
  POPFRACS15 <<- POPS15/sum(POPS15)
  # two doses plus two boosters
  MAX_DEMAND <<- sum(DEMAND15)*2
  
  # points at which costs shift per income level (ignoring the 80% mark)
  DEL_COST_THRESHOLDS <<- sapply(c(1,3)/10,function(x) x*POPS15) 
  
  
  # trial durations in normal times (change from years to weeks)
  OLD_DURATIONS <<- do.call(cbind,lapply(paste0('duration_',0:3),function(x)52*pardf[[x]]))
  
  
  # probability of success per phase
  POS <<- sapply(paste0('pos_',0:3), function(x) pardf[[x]])
  # probability to occur is the product of prior phases
  PTO = POS
  for(i in 2:4) PTO[,i] = apply(POS[,1:i],1,prod)
  PTO <<- PTO
  
  
  # costs of experienced and inexperienced manufacturers per phase
  EX <<- sapply(paste0('cost_',0:3,'_ex'),function(x) pardf[[x]])
  INEX <<- sapply(paste0('cost_',0:3,'_inex'),function(x) pardf[[x]])
  
  
  
  # times to scale up for different manufacturing types
  weeks_init = sapply(weeks_init_pars, function(x) pfixed[[x]])
  weeks_scale = sapply(weeks_scale_pars, function(x) pfixed[[x]])
  names(weeks_init) <- names(weeks_scale) <- MAN_TYPES
  WEEKS_INIT <<- weeks_init
  WEEKS_SCALE <<- weeks_scale
  
  # costs per dose delivery
  costperdose_names = apply(expand.grid(ics,PC_THRESHOLDS), 1, function(x) paste0('cost_',x[1],'_',as.numeric(x[2])))
  cost_per_dose = as.data.frame(sapply(costperdose_names,function(x) pardf[[x]]))
  # popsllmic = pops[3:4]/sum(pops[3:4])
  # for(i in PC_THRESHOLDS){
  #   lic = paste0('cost_lic_',as.numeric(i))
  #   lmic = paste0('cost_lmic_',as.numeric(i))
  #   llmic = paste0('cost_llmic_',as.numeric(i))
  #   cost_per_dose[[llmic]] = popsllmic[1]*cost_per_dose[[lmic]] + popsllmic[2]*cost_per_dose[[lic]]
  #   cost_per_dose[[lmic]] <- cost_per_dose[[lic]] <- NULL
  # }
  
  cost_per_dose_list = list()
  for(i in 1:NSAMPLES){
    cost_per_dose_list[[i]] = matrix(0,nrow=NLEVELS, ncol=length(PC_THRESHOLDS))
    for(il in 1:NLEVELS){
      lab1 = paste0('cost_',tolower(INCOMELEVELS)[il],'_')
      for(j in 1:length(PC_THRESHOLDS)){
        lab = paste0(lab1,as.numeric(PC_THRESHOLDS[j]))
        cost_per_dose_list[[i]][il,j] = cost_per_dose[[lab]][i]
      }
    }
  }
  COST_PER_DOSE_LIST <<- cost_per_dose_list
  
  
  list(rawpar, params)
  
}

## utility functions ##############################

# combine lic and lmic into llmic
combine_llmic = function(N_lic, N_lmic, doses){
  doses$LLMIC = (doses$LIC * N_lic + doses$LMIC * N_lmic)/(N_lmic + N_lic) 
  doses$LMIC <- doses$LIC <- NULL
  doses
}

# accumulating costs over time
accumulate = function(discount=0, from=1, to=15) {
  if(to<from){
    0
  }else if(to==from){
    1/(1+discount)^(to-1)
  }else{
    sum(sapply(from:to, function(x)1/(1+discount)^(x-1)))
  }
}


get_duration_scalar = function(phase_duration, old_durations){
  # phases are costs per week
  # reactive phases are costed by dividing the phase cost by its duration in normal times and
  # multiplying by the number of weeks in the expedited period
  scalars = rep(0,length(phase_duration))
  for(i in 1:length(phase_duration))
    scalars[i] = phase_duration[i]/old_durations[i]
  scalars
}


## delivery functions ###############################


dose_supply = function(week, w0, weeks_init, weeks_scale, man_cap){
  if((week - w0) <= weeks_init){
    doses = 0
  }else if((week - w0) <= (weeks_init + weeks_scale)){
    doses = man_cap/52 * (week - w0 - weeks_init ) / weeks_scale
  }else{
    doses = man_cap/52
  }
  return(doses)
}

get_ssv_supply = function(dm=365, capres=0, bpsv=F, weeks_init, weeks_scale){
  
  if(bpsv==T){
    weeks_init[['ex']] = pfixed$weeks_init_ex_bp
  }else{
    weeks_init[['ex']] = pfixed$weeks_init_ex_nb
  }
  
  time_to_approval = round(dm/7)
  w0 = time_to_approval - pfixed$week_trans_start
  man_cap_res = capres
  
  man_cap_ex = pfixed$man_curr - man_cap_res
  man_cap_bui = pfixed$man_glo - pfixed$man_curr
  
  man_cap = c(man_cap_res, man_cap_ex, man_cap_bui)
  names(man_cap) = MAN_TYPES
  supplies = list()
  for(x in MAN_TYPES) supplies[[x]] = rep(0,NWEEKS)
  for(y in 1:NYEARS){
    
      for(w in ((y-1)*52+1):min(y*52,NWEEKS)){
        if(sum(unlist(supplies))<MAX_DEMAND/1e9){
          for(x in MAN_TYPES){
            # if(w==1) print(c(x, weeks_init[[x]], weeks_scale[[x]], man_cap[[x]]))
            supplies[[x]][w] = dose_supply(week=w, 
                                           w0=w0, 
                                           weeks_init=weeks_init[[x]], 
                                           weeks_scale=weeks_scale[[x]], 
                                           man_cap=man_cap[[x]])
        }
      }
    }
  }
  cumulative = cumsum(Reduce("+",supplies))
  
  
  ## return
  list(supplies=supplies,
       cumulative=cumulative)
}


annualise_ssv_procurement = function(sup){
  # supplied doses
  # sup = supplies[[s]]
  res = sup[[1]]
  exbui = sup[[2]] + sup[[3]]
  # total doses per year
  annual_res_doses = sapply(1:NYEARS,function(y) sum(na.omit(res[(1:52) + (y-1)*52])))
  annual_exbui_doses = sapply(1:NYEARS,function(y) sum(na.omit(exbui[(1:52) + (y-1)*52])))
  annual_doses = annual_res_doses + annual_exbui_doses
  booster1 <- booster2 <- ssv_doses <- rep(0,NYEARS)
  for(y in 1:NYEARS){
    ssv_doses[y] = min(annual_doses[y], sum(DEMAND15)/1e9-sum(ssv_doses))
  }
  booster_doses = annual_doses - ssv_doses
  for(y in 2:NYEARS){
    booster1[y] = min(ssv_doses[y-1]/2, cumsum(booster_doses)[y]-sum(booster1[1:(y-1)]))
  }
  for(y in 3:NYEARS){
    booster2[y] = min(booster1[y-1], cumsum(booster_doses)[y]-sum(booster1[1:y])-sum(booster2[1:(y-1)]))
  }
  ssv_py_inc_booster = ssv_doses + booster1 + booster2
}


get_bpsv_supply = function(){
  
  ## manufacture
  
  man_cap_bpsv = pfixed$man_curr - pfixed$hic_cap_res
  bpsv_supplies <- c()
  w = 0
  while(sum(bpsv_supplies)<TOTAL_BPSV_DEMAND){ 
    w = w+1
    bpsv_supplies[w] = dose_supply(week=w, w0=0, weeks_init=WEEKS_INIT[['res']], weeks_scale=WEEKS_SCALE[['res']], man_cap=man_cap_bpsv)
  }
  
  b_cumulative = cumsum(bpsv_supplies)
  bpsv_weeks = w
  
  b_supply <- data.frame(doses=b_cumulative,week=1:bpsv_weeks)
  
  ## allocation
  b_cumulative = b_cumulative + pfixed$bpsv_inv_res/1e9
  
  b_allocations = t(POPFRACS65 %*% t(b_cumulative))
  
  b_alloc = as.data.frame(b_allocations)
  b_alloc$Week = 1:bpsv_weeks
  colnames(b_alloc)[1:NLEVELS] = INCOMELEVELS
  
  ## delivery
  
  bpsv_dose_delivery <- list()
  for(il in 1:NLEVELS){
    max_doses_each_wk = basic_rates[il] * POPS0[il]/1e9
    # new doses in per week
    doses_per_week = diff(c(0,b_alloc[,il]))
    doses <- rep(0,bpsv_weeks)
    doses_left = 0
    for(w in 1:bpsv_weeks){
      # doses this week = stock plus flow
      doses_left = doses_left + doses_per_week[w]
      if(w>pfixed$duration_3_resp)
        # doses given are the minimum of: the delivery rate; the doses left (per population); the fraction of the population still unvaccinated
        doses[w] =  min(max_doses_each_wk, doses_left, max(0,DEMAND65[il] - sum(doses)) )
      # subtract doses given to update the stock
      doses_left = doses_left - doses[w]
    }
    bpsv_dose_delivery[[il]] = cumsum(doses)/(POPS0[il]/1e9)
  }
  bpsv_dose_delivery = cbind(Week = 1:bpsv_weeks, 
                             as.data.frame(do.call(cbind,bpsv_dose_delivery)))
  colnames(bpsv_dose_delivery)[1+1:NLEVELS] = INCOMELEVELS
  
  ## return
  list(b_supply=b_supply,
       b_alloc=b_alloc,
       bpsv_dose_delivery=bpsv_dose_delivery)
}

## ssv allocation: reserved doses
allocate_res_doses = function(cum_received, new_doses){
  if(cum_received[1] < pfixed$hic_cap_res){
    allocation = c(new_doses, 0, 0, 0)
  }else if(cum_received[1] < DEMAND15[1]/1e9){
    allocation = new_doses * POPFRACS15
  }else if(cum_received[2] < DEMAND15[2]/1e9){
    allocation = c(0, new_doses * POPS15[2:NLEVELS]/sum(POPS15[2:NLEVELS]))
  }else if(cum_received[3] < DEMAND15[3]/1e9){
    allocation = c(0, 0, new_doses * POPS15[3:NLEVELS]/sum(POPS15[3:NLEVELS]))
  }else if(cum_received[4] < DEMAND15[4]/1e9){
    allocation = c(0,0,0,new_doses)
  }else{
    allocation = rep(0,NLEVELS)
  }
  return(allocation)
}

allocate_exbui_doses = function(cum_received, new_doses){
  if(cum_received[1] < DEMAND15[1]/1e9){
    allocation = c(new_doses, 0, 0, 0)
  }else if(cum_received[2] < DEMAND15[2]/1e9){
    allocation = c(0, new_doses, 0, 0)
  }else if(cum_received[3] < DEMAND15[3]/1e9){
    allocation = c(0,0,new_doses, 0)
  }else if(cum_received[4] < DEMAND15[4]/1e9){
    allocation = c(0,0, 0,new_doses)
  }else{
    allocation = rep(0,NLEVELS)
  }
  return(allocation)
}

allocate_and_deliver_doses = function(allocation_functions, supplies, del_rates, time_to_approval){
  
  ## allocation
  # allocate doses one week at a time
  cumulative_doses = matrix(0,nrow=NWEEKS,ncol=NLEVELS)
  for(w in 2:NWEEKS){
    current_count = cumulative_doses[w-1,]
    
    for(x in MAN_TYPES){
      new_doses = supplies[[x]][w]
      allocation_res = allocation_functions[[x]](cum_received = current_count, new_doses=new_doses)
      current_count = current_count + allocation_res
    }
    
    cumulative_doses[w,] = current_count
  }
  alloc = as.data.frame(cumulative_doses)
  alloc$Week = 1:NWEEKS
  colnames(alloc)[1:NLEVELS] = INCOMELEVELS
  
  ## delivery
  
  second_dose_delivery <- all_dose_delivery <- list()
  # one income level at a time
  for(il in 1:length(INCOMELEVELS)){
    max_doses_each_wk = del_rates * POPS0[il]/1e9
    # new doses in per week
    doses_per_week = c(0, diff(alloc[,il]))
    first_doses <- second_doses <- alldoses <- rep(0,NWEEKS)
    doses_left = 0
    for(w in 1:NWEEKS){
      max_doses_this_week = max_doses_each_wk
      # doses this week = stock plus flow
      doses_left = doses_left + doses_per_week[w]
      if(w > time_to_approval){
        if(w>4){
          # second doses given are the minimum of: the delivery rate; the doses left (per population); the fraction of the population vaccinated once, at least three weeks ago, but not twice.
          second_doses[w] = min(max_doses_this_week, doses_left, max(0, sum(first_doses[1:(w-4)])-sum(second_doses)))
          # subtract doses given to update the stock
          doses_left = doses_left - second_doses[w]
          # update max this week
          max_doses_this_week = max_doses_each_wk - second_doses[w]
        }
        # first doses given are the minimum of: the (remaining) delivery rate; the doses left (per population); the fraction of the population still unvaccinated
        first_doses[w] =  min(max_doses_this_week, doses_left, max(0,DEMAND15[il]/1e9/2 - sum(first_doses)) )
        # subtract doses given to update the stock
        doses_left = doses_left - first_doses[w]
        alldoses[w] = second_doses[w]+first_doses[w]
      }
    }
    # all_dose_delivery in billions
    all_dose_delivery[[il]] = cumsum(alldoses)
    # second_dose_delivery in percentage
    second_dose_delivery[[il]] = data.frame(doses=cumsum(second_doses)/(POPS0[il]/1e9), il=INCOMELEVELS[il], Week=1:NWEEKS)
  }
  second_dose_delivery = do.call(rbind,second_dose_delivery)
  all_dose_delivery = as.data.frame(do.call(cbind,all_dose_delivery))
  all_dose_delivery$Week = 1:NWEEKS
  
  ## return
  list(alloc=alloc,
       second_dose_delivery=second_dose_delivery,
       all_dose_delivery=all_dose_delivery)
}


## all costs #####################################

get_all_costs = function(time_to_rd, timescale, 
                         cap_cost_per_dose, cap_res_vol, 
                         ssv_py_inc_booster,
                         all_ssv_delivered,
                         bpsv_dose_delivery){
  
  
}


## ssv costs ######################################

get_enabling_costs = function(cost_enab, discount=0, time_to_target_dm = 0){
  rd_years = time_to_target_dm
  cost_enab * accumulate(discount, from=1, to=rd_years)
}


get_cap_res_cost = function(cap_cost_per_dose, cap_res_vol = 0){
  (cap_res_vol + pfixed$hic_cap_res)*cap_cost_per_dose
}


get_ssv_randd_costs = function(pos, pto=NULL, ex, timescale, icost_lic, n_ssv_candidates){
  
  ## trial success-weighted costs
  
  if(is.null(pto)){
    pto = pos
    for(i in 2:length(pos)) pto[i] = prod(pos[1:i])
  }
  
  # add 1 column for probability of phase 0 to occur
  pto = c(1, pto)
  
  ssv_phasecost = ex * timescale 
  
  # add licence cost column (adjusted for inflation)
  ssv_phase_liccost = c(ssv_phasecost, icost_lic)
  
  # adjust for pto
  ssv_phase_pto = ssv_phase_liccost*pto
  
  ssv_rd_costsamples = n_ssv_candidates * sum(ssv_phase_pto)
  
  ssv_rd_costsamples
}


get_ssv_procurement_costs = function(ssv_py_inc_booster, cost_res, cost_un, cap_res_vol=0){
  # capacity distribution
  man_cap_res = pfixed$hic_cap_res + cap_res_vol
  # price according to utilisation of reserved capacity
  cost_per_year = sapply(1:NYEARS, function(y) 
    min(man_cap_res,ssv_py_inc_booster[y]) * cost_res 
    + max(ssv_py_inc_booster[y]-man_cap_res, 0) * cost_un
  )
  cost_per_year
}
## cost delivery

get_ssv_delivery_costs = function(all_ssv_delivered, cost_per_dose, discount){
  ssv_delivery_costs = 0
  cost_per_dose = cost_per_dose
  costvec = rep(1, nrow(all_ssv_delivered))
  for(j in 1:length(INCOMELEVELS)){
    scendosesj = all_ssv_delivered[,j]
    # the jth column is HIC, UMIC or LLMIC
    popj = POPS15[j]/1e9
    ##!! should use DEL_COST_THRESHOLDS, in total pop number
    ind1 = which(scendosesj/popj > 0.1)
    ind3 = which(scendosesj/popj > 0.3)
    costvec[ind1[1]:ind3[1]] <- 2
    costvec[ind3[1]:length(costvec)] <- 3 
    costsamples = rep(0,length(costvec))
    cost_per_dose_row = cost_per_dose[j, ]
    costsamples = cost_per_dose_row[costvec]
    
    weeklycountrycosts = costsamples*c(0,diff(scendosesj))
    annualcountrycosts = c()
    for(i in 1:NYEARS) annualcountrycosts[i] = sum(weeklycountrycosts[c(1:52)+(i-1)*52])
    # discount to starting year
    dcountrycosts = annualcountrycosts/(1+discount)^c(0:(NYEARS-1))
    countrycosts = sum(dcountrycosts)
    # add country levels together, multiplying by population size
    ssv_delivery_costs = ssv_delivery_costs + countrycosts
  }
  ssv_delivery_costs
}

## bpsv costs ########################

get_bpsv_costs = function(total_bpsv, pos, pto, ex, inex, 
                          inex_weight = 0.875,
                          n_bpsv_candidates = 8,
                          durations = c(2,2,2),
                          old_duration = 104,
                          discount = 0,
                          icost_lic = 368320,
                          cost_bpsv_res = 1012,
                          cost_res = 6.29,
                          bpsv_replenishment = 3,
                          cost_cogs = 4.83,
                          profit = 0.2,
                          cost_travel = 0.12,
                          cost_ff = 0.14){
  
  if(is.null(pto)){
    pto = pos
    for(i in 2:length(pos)) pto[i] = prod(pos[1:i])
  }
  
  # add 1 column for probability of phase 0 to occur
  pto = c(1, pto)
  
  # cost per phase is a weighted sum
  inex_costs = inex_weight * inex + (1-inex_weight) * ex
  
  d0 = durations[1]
  d1 = durations[2]
  dc = 1/(1+discount)
  time_scalar = c(1, dc^d0, dc^(d0+d1))
  firstthreephases = inex_costs[1:3]*pto[1:3]
  bpsv_rd_costsamples_dyd = n_bpsv_candidates*sum(firstthreephases*time_scalar)/1e9
  bpsv_rd_costsamples_no_d = n_bpsv_candidates*sum(firstthreephases)/1e9
  
  ## bpsv reactive r&d costs
  bpsvresrd = n_bpsv_candidates * pto[4] * (pfixed$duration_3_resp/old_duration * ex[4] + pos[4]*icost_lic)
  
  ## investigational reserve costs
  
  # cost for bpsv reserve per dose: subtract fill--finish cost, and add on profit
  cost_bpsvinv = cost_cogs*(1-cost_ff)*(1+profit)
  
  # cost per year for the inventory: replenished every three years.
  bpsv_inv_res = pfixed$bpsv_inv_res
  inv_cost_per_year = (cost_bpsvinv/bpsv_replenishment*bpsv_inv_res + cost_bpsv_res)/1e9
  # time to completion of phase 2
  time_to_bpsv = d0 + d1 + durations[3]
  
  ## dose procurement cost
  # fill and finish
  bpsv_ff_trans = (cost_ff+cost_travel)*cost_cogs*(1+profit)
  # procurement cost
  bpsvproc = cost_res*total_bpsv + bpsv_inv_res*bpsv_ff_trans/1e9
  
  # return
  list(bpsv_rd_costsamples_dyd=bpsv_rd_costsamples_dyd,
       bpsv_rd_costsamples_no_d=bpsv_rd_costsamples_no_d,
       bpsvresrd=bpsvresrd,
       inv_cost_per_year=inv_cost_per_year,
       bpsvproc=bpsvproc,
       time_to_bpsv=time_to_bpsv)
  
}

get_bpsv_del_costs = function(bpsv_dose_delivery, cost_per_dose_mat){
  bpsv_del_cost = 0
  for(il in 1:length(INCOMELEVELS)){
    icthresh = DEL_COST_THRESHOLDS[il,]
    ildoses = bpsv_dose_delivery[[INCOMELEVELS[il]]]
    # bpsv_dose_delivery is a fraction of population
    icdemand = max(ildoses)*POPS0[il] # DEMAND65[il]
    cost_per_dose_il = cost_per_dose_mat[il,]
    firstcost = cost_per_dose_il[1]
    if(icdemand < icthresh[1]){
      cost = icdemand * firstcost
    }else if(icdemand < icthresh[2]){
      cost = icthresh[1] * firstcost + (icdemand - icthresh[1]) * cost_per_dose_il[2]
    }else{
      cost = icthresh[1] * firstcost + (icthresh[2] - icthresh[1]) * cost_per_dose_il[2]  + (icdemand - icthresh[2]) * cost_per_dose_il[3]
    }
    # print(round(summary(cost/1e9)))
    bpsv_del_cost = bpsv_del_cost + cost
  }
  bpsv_del_cost
}





generate_scenario = function(days_mission=365, 
                             capacity_reservation=0, 
                             BPSV=F, 
                             pop_proportional_allocation=F, 
                             vaccination_rates=c(.07,.07,.02)){
  
  if(length(days_mission)>1){
    warning('Only first value used for variable days_mission.',immediate.=T)
    days_mission = days_mission[1]
  }
  if(length(capacity_reservation)>1){
    warning('Only first value used for variable capacity_reservation',immediate.=T)
    capacity_reservation = capacity_reservation[1]
  }
  if(length(BPSV)>1){
    warning('Only first value used for variable BPSV',immediate.=T)
    BPSV = BPSV[1]
  }
  if(length(pop_proportional_allocation)>1){
    warning('Only first value used for variable pop_proportional_allocation',immediate.=T)
    pop_proportional_allocation = pop_proportional_allocation[1]
  }
  if(length(vaccination_rates)>3){
    warning('Only first three values used for variable vaccination_rates',immediate.=T)
    vaccination_rates = vaccination_rates[1:3]
  }
  
  dms = c(365,200,100)
  if(!days_mission%in%dms){
    stop('days_mission should be 365, 200, or 100')
  }else{
    dm_index = which(dms==days_mission)
  }
  
  if(!BPSV%in%c(T,F))
    stop('BPSV should be T or F')
  if(!pop_proportional_allocation%in%c(T,F))
    stop('pop_proportional_allocation should be T or F')
  
  
  
  list(inputs=list(days_mission,capacity_reservation,BPSV,pop_proportional_allocation,vaccination_rates),
       parameters=parameters, delivery=delivery, costs = costs)
  
}
# generate_scenario(days_mission =2){

# if(any(crs + pfixed$$hic_cap_res<0)) stop('Capacity reservation must be at least 0')
# if(any(crs + pfixed$$hic_cap_res>9)) error('Capacity reservation must be at most 9')
# }
# 
generate_cepi_scenarios = function(){
  
  
}#






