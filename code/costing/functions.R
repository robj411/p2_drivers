
## setting up ######################################

get_parameters = function(nsamples = 100){
  
  ## global variables ############################
  
  NSAMPLES <<- nsamples
  MAN_TYPES <<- c('res','ex','bui')
  PC_THRESHOLDS <<- c(0,11,31)
  INCOMELEVELS <<- c('HIC','UMIC','LMIC','LIC')
  NLEVELS <<- length(INCOMELEVELS)
  
  
  ## parameters 
  
  
  params <- rawpar <- list()
  
  for(i in 1:nrow(caps)){
    parameters <- as.numeric(strsplit(caps$Parameters[i],',')[[1]])
    # bounds <- as.numeric(strsplit(caps$Bounds[i],',')[[1]])
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
      distcall = function(nsamples) rbetapr(nsamples, shape1 = parameters[1], shape2 = parameters[2], scale = parameters[3])
      # distcall = function(nsamples) rpearsonVI(nsamples, a = parameters[1], b = parameters[2], scale = parameters[3], location=0)
    }else if(caps$Distribution[i]=='Triangular'){
      distcall = function(nsamples) rtri(nsamples,min=parameters[1],mode=parameters[2],max=parameters[3])
    }else if(caps$Distribution[i]=='Multinomial'){
      values = na.omit(c(parameters))
      distcall = function(nsamples) sample(values,size=nsamples,replace=T)
    }
    
    if(caps$Distribution[i]!='Constant'){
      samples = distcall(NSAMPLES)
      # if(length(bounds)>1){
      #   fails = samples < bounds[1] | samples > bounds[2]
      #   while(sum(fails)>0){
      #     samples[fails] = distcall(sum(fails))
      #     fails = samples < bounds[1] | samples > bounds[2]
      #   }
      # }
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
  covid_candidates <- sapply(paste0('covid_ssv_',0:4), function(x) rawpar[[x]])
  phase_dur_par = apply(expand.grid(0:3,dms), 1, function(x) paste0('weeks_P',x[1],'_',as.numeric(x[2])))
  weeks_init_pars = paste0('weeks_init_',MAN_TYPES)
  weeks_init_ex_pars = paste0('weeks_init_ex_',c('nb','bp'))
  weeks_scale_pars = paste0('weeks_scale_',MAN_TYPES)
  fixedpars = c(pop0pars, pop15pars, pop65pars,'final_vaccine_coverage','vaccine_wastage', # volumes made
                phase_dur_par,'duration_3_resp', # phase durations
                weeks_init_pars,weeks_init_ex_pars,weeks_scale_pars,'week_trans_start', # weeks to manufacturing
                'man_curr','man_glo', # manufacturing capacities
                'bpsv_inv_res','hic_cap_res', # reserve amounts
                'n_boosters', covid_candidates) 
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
  TOTAL_BPSV_DEMAND <<- sum(DEMAND65/1e9) # 1 # billion. not: 
  # vaccine wastage applied to BPSV only
  # vaccine_wastage = 1 - pfixed$final_vaccine_coverage *sum(POPS65)/1e9
  
  DEMAND15 <<- pfixed$final_vaccine_coverage*POPS15*2 #/ (1-pfixed$vaccine_wastage)
  POPFRACS15 <<- POPS15/sum(POPS15)
  # two doses plus two boosters
  MAX_DEMAND <<- sum(DEMAND15)*2
  
  # points at which costs shift per income level (ignoring the 80% mark)
  DEL_COST_THRESHOLDS <<- sapply(c(1,3)/10,function(x) x*POPS15*2) 
  
  
  # trial durations in normal times (change from years to weeks)
  OLD_DURATIONS <<- do.call(cbind,lapply(paste0('duration_',0:3),function(x)52*pardf[[x]]))
  
  
  # probability of success per phase for bpsv
  POS_BPSV <<- matrix(sapply(paste0('pos_',0:3), function(x) pardf[[x]]), ncol=4, byrow=F)
  # probability to occur is the product of prior phases
  PTO_BPSV = POS_BPSV
  for(i in 2:4) PTO_BPSV[,i] = apply(POS_BPSV[,1:i, drop = FALSE],1,prod)
  PTO_BPSV <<- PTO_BPSV
  
  # probability of success per phase for ssv
  pos_ssv = matrix(0,ncol=4,nrow=NSAMPLES)
  for(i in 1:4){
    successes = sum(covid_candidates[i:5]) - covid_candidates[i]
    failures = covid_candidates[i]
    pos_ssv[,i] = rbeta(NSAMPLES, successes+1, failures+1)
  }
  POS_SSV <<- pos_ssv
  # probability to occur is the product of prior phases
  PTO_SSV = POS_SSV
  for(i in 2:4) PTO_SSV[,i] = apply(POS_SSV[,1:i, drop = FALSE],1,prod)
  PTO_SSV <<- PTO_SSV
  
  
  
  
  # costs of experienced and inexperienced manufacturers per phase
  EX <<- matrix(sapply(paste0('cost_',0:3,'_ex'),function(x) pardf[[x]]*pardf$inflation),ncol=4,byrow=F)
  INEX <<- matrix(sapply(paste0('cost_',0:3,'_inex'),function(x) pardf[[x]]*pardf$inflation),ncol=4,byrow=F)
  
  
  
  # times to scale up for different manufacturing types
  weeks_init = sapply(weeks_init_pars, function(x) pfixed[[x]])
  weeks_scale = sapply(weeks_scale_pars, function(x) pfixed[[x]])
  names(weeks_init) <- names(weeks_scale) <- MAN_TYPES
  WEEKS_INIT <<- weeks_init
  WEEKS_SCALE <<- weeks_scale
  
  # assume one year to start, time to cover all populations (assuming lowest rate 0.02), and number of boosters:
  NYEARS <<- ceiling(1 + NLEVELS*pfixed$final_vaccine_coverage/0.02/52 + pfixed$n_boosters)
  
  NWEEKS <<- ceiling(NYEARS*52)
  
  
  # costs per dose delivery
  costperdose_names = apply(expand.grid(ics,PC_THRESHOLDS), 1, function(x) paste0('cost_',x[1],'_',as.numeric(x[2])))
  cost_per_dose = as.data.frame(sapply(costperdose_names,function(x) pardf[[x]], simplify = FALSE))
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
  man_cap_res = pfixed$hic_cap_res + capres
  
  man_cap_ex = pfixed$man_curr - man_cap_res
  man_cap_bui = pfixed$man_glo - pfixed$man_curr
  
  man_cap = c(man_cap_res, man_cap_ex, man_cap_bui)
  names(man_cap) = MAN_TYPES
  supplies = list()
  for(x in MAN_TYPES) supplies[[x]] = rep(0,NWEEKS)
  for(w in 1:NWEEKS){
    # keep production going; allocation will be limited by demand
    # if(sum(unlist(supplies))<MAX_DEMAND*2/1e9){
      for(x in MAN_TYPES){
        supplies[[x]][w] = dose_supply(week=w, 
                                       w0=w0, 
                                       weeks_init=weeks_init[[x]], 
                                       weeks_scale=weeks_scale[[x]], 
                                       man_cap=man_cap[[x]])
      }
    # }
  }
  cumulative = cumsum(Reduce("+",supplies))
  
  
  ## return
  list(supplies=supplies,
       cumulative=cumulative)
}


annualise_ssv_procurement = function(supplies, alloc, cap_res){
  # supplied doses are per week
  # browser()
  res = supplies[[1]]
  exbui = supplies[[2]] + supplies[[3]]
  # allocated doses are cumulative
  total_alloc = rowSums(alloc[,1:4])
  
  ssv_first_schedule = diff(c(0, sapply(1:NYEARS,function(y) total_alloc[y*52])))
  # print(sum(ssv_first_schedule))
  first_schedule_completes = which(total_alloc >= sum(ssv_first_schedule))[1]
  # print(DEMAND15/1e9)
  # print(alloc[first_schedule_completes,])
  
  res_fs = res[1:first_schedule_completes]
  un_fs = exbui[1:first_schedule_completes]
  
  # max production per year for fs
  max_annual_res_fs = sapply(1:NYEARS,function(y) sum(na.omit(res_fs[(1:52) + (y-1)*52])))
  max_annual_exbui_fs = sapply(1:NYEARS,function(y) sum(na.omit(un_fs[(1:52) + (y-1)*52])))
  
  fs_res = pmin(ssv_first_schedule, max_annual_res_fs)
  fs_un = pmin(max_annual_exbui_fs, ssv_first_schedule - fs_res)
  
  # max production per year
  max_annual_res_doses = sapply(1:NYEARS,function(y) sum(na.omit(res[(1:52) + (y-1)*52])))
  max_annual_exbui_doses = sapply(1:NYEARS,function(y) sum(na.omit(exbui[(1:52) + (y-1)*52])))
  
  res_left_for_booster = max_annual_res_doses - fs_res
  un_left_for_booster = max_annual_exbui_doses - fs_un
  
  booster_res <- booster_un <- rep(0,NYEARS)
  for(y in 2:NYEARS){
    boosters_due = sum(ssv_first_schedule[max(1,y-2):(y-1)])/2
    booster_res[y] = min(res_left_for_booster[y], boosters_due)
    booster_un[y] = min(un_left_for_booster[y], boosters_due - booster_res[y])
  }
  
  annual_res_doses = booster_res + fs_res
  annual_exbui_doses = booster_un + fs_un
  # print(sum(annual_res_doses + annual_exbui_doses))
  # print(sum(fs_res+fs_un))
  # print(sum(fs_res))
  # print(sum(fs_un))
  # print(sum(booster_res))
  # print(sum(booster_un))
  list(res=annual_res_doses, un=annual_exbui_doses)
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
  #  t(POPFRACS65 %*% t(bpsv_supplies))[13:23,]*1000
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
allocate_res_doses = function(cum_received, new_doses, hic_only=1){
  demand15_bn = DEMAND15/1e9
  allocation = rep(0,NLEVELS)
  
  first_recipient = which(cum_received < demand15_bn)[1]
  
  if(!is.na(first_recipient)){
    shared_doses = new_doses
    hic_only_doses = 0 
    
    if(first_recipient==1){ # allocate fraction hic_only to hic
      shared_doses = (1 - hic_only) * new_doses
      hic_only_doses = hic_only * new_doses
    }
    
    allocation[first_recipient:NLEVELS] = shared_doses * POPS15[first_recipient:NLEVELS]/sum(POPS15[first_recipient:NLEVELS])
    allocation[1] = allocation[1] + hic_only_doses
  }
  
  return(allocation)
}

allocate_exbui_doses = function(cum_received, new_doses, hic_only=0){
  demand15_bn = DEMAND15/1e9
  allocation = rep(0,NLEVELS)
  
  sole_recipient = which(cum_received < demand15_bn)[1]
  if(!is.na(sole_recipient))
    allocation[sole_recipient] = new_doses
  
  return(allocation)
}

allocate_and_deliver_doses = function(allocation_functions, supplies, del_rates, time_to_approval, hic_only){
  
  ## allocation
  # allocate doses one week at a time
  cumulative_doses = matrix(0,nrow=NWEEKS,ncol=NLEVELS)
  for(w in 2:NWEEKS){
    current_count = cumulative_doses[w-1,]
    week_count = rep(0, NLEVELS)
    for(x in MAN_TYPES){
      new_doses = supplies[[x]][w]
      allocation_res = allocation_functions[[x]](cum_received = current_count, 
                                                 new_doses=new_doses, 
                                                 hic_only = ifelse(x=='res',hic_only,0))
      # week_count = week_count + allocation_res
      current_count = current_count + allocation_res
    }
    # cumulative_doses[w,] = current_count + week_count
    cumulative_doses[w,] = current_count
  }
  alloc = as.data.frame(cumulative_doses)
  alloc$Week = 1:NWEEKS
  colnames(alloc)[1:NLEVELS] = INCOMELEVELS
  
  ## delivery
  
  second_dose_delivery <- all_dose_delivery <- list()
  # one income level at a time
  for(il in 1:length(INCOMELEVELS)){
    max_doses_each_wk = del_rates[il] * POPS0[il]/1e9
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
          second_doses[w] = min(max_doses_this_week, doses_left, first_doses[w-4])
                                # max(0, sum(first_doses[1:(w-4)])-sum(second_doses)))
          # subtract doses given to update the stock
          doses_left = doses_left - second_doses[w]
          # update max this week
          max_doses_this_week = max_doses_each_wk - second_doses[w]
        }
        # first doses given if
        # they are the minimum of: the (remaining) delivery rate; the doses left (per population); the fraction of the population still unvaccinated
        if(sum(first_doses) < DEMAND15[il]/1e9/2){
          first_doses[w] =  min(max_doses_this_week, doses_left, max(0, DEMAND15[il]/1e9/2-sum(first_doses))) 
        }else{
          first_doses[w] = 0
        }
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
  (cap_res_vol + pfixed$hic_cap_res)*cap_cost_per_dose * 1e3 # billions to millions
}


get_ssv_randd_costs = function(pos, pto=NULL, exi, timescale, cost_lic, n_ssv_successes, q_ssv_successes){
  
  ## trial success-weighted costs
  
  if(is.null(pto)){
    pto = pos
    for(i in 2:length(pos)) pto[i] = prod(pos[1:i])
  }
  
  # number of candidates to achieve target number of successes with given probability
  n_ssv_candidates = n_ssv_successes + qnbinom(q_ssv_successes,size = n_ssv_successes, prob = pto[4])
  
  # add 1 column for probability of phase 0 to occur
  pto = c(1, pto)
  
  ssv_phasecost = exi * timescale 
  
  # add licence cost column 
  ssv_phase_liccost = c(ssv_phasecost, cost_lic)
  
  # adjust for pto
  ssv_phase_pto = ssv_phase_liccost*pto
  
  ssv_rd_costsamples = n_ssv_candidates * sum(ssv_phase_pto)
  
  ssv_rd_costsamples/1e6
}


get_ssv_procurement_costs = function(ssv_py_inc_booster, cost_res, cost_un){
  # capacity distribution
  # man_cap_res = pfixed$hic_cap_res + cap_res_vol
  # price according to utilisation of reserved capacity
  # cost_per_year = sapply(1:NYEARS, function(y) 
  #   min(man_cap_res,ssv_py_inc_booster[y]) * cost_res 
  #   + max(ssv_py_inc_booster[y]-man_cap_res, 0) * cost_un
  # )
  cost_per_year = ssv_py_inc_booster$res * cost_res + ssv_py_inc_booster$un * cost_un
  cost_per_year*1e3 # from billions to millions
}
## cost delivery

thresholded_costs = function(costs, thresholds, doses){
  if(doses <= 0){
    cost = 0
  }else{
    toolow = thresholds<0
    toohigh = thresholds>doses
    weights = diff(c(0, thresholds[!toolow&!toohigh], doses))
    if(any(toolow)) costs = costs[-which(toolow)]
    cost = weights * costs[1:length(weights)]
    # print(cost)
  }
  sum(cost)
}

get_ssv_delivery_costs = function(all_ssv_delivered, cost_per_dose, discount){
  ssv_delivery_cost_mat = matrix(0,nrow=NLEVELS,ncol=NYEARS)
  for(j in 1:length(INCOMELEVELS)){
    costvec = rep(1, nrow(all_ssv_delivered))
    scendosesj = pmin(all_ssv_delivered[,j], DEMAND15[j]/1e9)
    pop_thresh = DEL_COST_THRESHOLDS[j,]/1e9
    cost_per_dose_row = cost_per_dose[j, ]
    
    annualdoses <- annualcosts <- c()
    for(i in 1:NYEARS) {
      annualdoses[i] = scendosesj[i*52]
      doses_so_far = ifelse(i>1,annualdoses[i-1],0)
      annualcosts[i] <- thresholded_costs(costs=cost_per_dose_row, 
                                          thresholds=pop_thresh-doses_so_far, 
                                          doses=annualdoses[i]-doses_so_far)
    }
    # print(diff(annualdoses))
    # print(annualcosts)
    
    ssvperyear = diff(c(0, annualdoses))    
    booster_doses <- booster_costs <- rep(0,NYEARS)
    for(i in 2:NYEARS){
      booster_doses[i] = sum(ssvperyear[max(i-2,1):(i-1)])/2
      booster_costs[i] = booster_doses[i]*cost_per_dose_row[3]
    }
    
    ssv_delivery_cost_mat[j,] = annualcosts+booster_costs
  }
  # print(ssv_delivery_cost_mat)
  # discount to starting year # from billions to millions
  dcountrycosts = (colSums(ssv_delivery_cost_mat*1e3))/(1+discount)^c(0:(NYEARS-1))
  dis_countrycosts = sum(dcountrycosts)
  list(countrycosts=sum(ssv_delivery_cost_mat*1e3),
       dis_countrycosts=dis_countrycosts
  )
}

## bpsv costs ########################

get_bpsv_costs = function(total_bpsv, pos, pto, ex, exi, inexi, 
                          inex_weight = 0.875,
                          n_bpsv_candidates = 8,
                          n_bpsv_p1 = 0,
                          bpsv_res_upfront = 0.115,
                          y_durations = c(2,2,2),
                          old_duration = 104,
                          discount = 0,
                          cost_lic = 287750,
                          cost_bpsv_res = 0.01012,
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
  inex_costs = inex_weight * inexi + (1-inex_weight) * exi
  
  rep_year = function(cost,pto,dur) rep(cost*pto/dur,dur)
  
  dc = 1/(1+discount)
  inex_costs_per_year = c(rep_year(inex_costs[1],pto[1], y_durations[1]), rep_year(inex_costs[2],pto[2], y_durations[2]), rep_year(inex_costs[3],pto[3], y_durations[3]))
  d_inex_costs_per_year = inex_costs_per_year * dc^(0:(length(inex_costs_per_year)-1))
  inex_costs_per_year_from1 = c(rep_year(exi[2],1, y_durations[2]), rep_year(exi[3],pos[2], y_durations[3]))
  d_inex_costs_per_year_from1 = inex_costs_per_year_from1 * dc^(0:(length(inex_costs_per_year_from1)-1))
  
  bpsv_rd_costsamples_dyd_py = ((n_bpsv_candidates - n_bpsv_p1)*d_inex_costs_per_year + 
                               n_bpsv_p1*d_inex_costs_per_year_from1)/1e6 # million
  bpsv_rd_costsamples_no_d_py = ((n_bpsv_candidates - n_bpsv_p1)*inex_costs_per_year + 
                                n_bpsv_p1*inex_costs_per_year_from1)/1e6 # million
  bpsv_rd_costsamples_dyd = sum(bpsv_rd_costsamples_dyd_py)
  bpsv_rd_costsamples_no_d = sum(bpsv_rd_costsamples_no_d_py)
  
  ## bpsv reactive r&d costs
  bpsvresrd = n_bpsv_candidates * pto[4] * (pfixed$duration_3_resp/old_duration * exi[4] + pos[4]*cost_lic)/1e6
  # bpsvresrd = (pfixed$duration_3_resp/old_duration * ex[4] + cost_lic)/1e6
  
  ## investigational reserve costs
  
  # cost for bpsv reserve per dose: subtract fill--finish cost, and add on profit
  cost_bpsvinv = cost_cogs*(1-cost_ff)*(1+profit)
  
  # cost per year for the inventory: replenished every three years.
  bpsv_inv_res = pfixed$bpsv_inv_res
  # in millions
  inv_cost_per_year = (cost_bpsvinv/bpsv_replenishment + cost_bpsv_res)*bpsv_inv_res/1e6
  # time to completion of phase 2
  time_to_bpsv = sum(y_durations)
  # discounted upfront cost
  bpsv_upfront = bpsv_res_upfront*bpsv_inv_res/1e6 # in million USD
  dis_bpsv_upfront = bpsv_upfront/(1+discount)^time_to_bpsv
  
  ## dose procurement cost
  # fill and finish
  bpsv_ff_trans = (cost_ff+cost_travel)*cost_cogs*(1+profit)
  # procurement cost
  bpsvproc = cost_res*total_bpsv*1e3 + bpsv_inv_res*bpsv_ff_trans/1e6
  
  # return
  list(bpsv_rd_costsamples_dyd=bpsv_rd_costsamples_dyd,
       bpsv_rd_costsamples_no_d=bpsv_rd_costsamples_no_d,
       bpsv_rd_costsamples_dyd_py=bpsv_rd_costsamples_dyd_py,
       bpsv_rd_costsamples_no_d_py=bpsv_rd_costsamples_no_d_py,
       bpsvresrd=bpsvresrd,
       inv_cost_per_year=inv_cost_per_year,
       bpsvproc=bpsvproc,
       time_to_bpsv=time_to_bpsv,
       bpsv_upfront=bpsv_upfront,
       dis_bpsv_upfront=dis_bpsv_upfront)
  
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
  unname(bpsv_del_cost)/1e6
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






