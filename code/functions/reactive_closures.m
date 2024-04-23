% function implementing reactive closure mitigation strategies
%
% t: current time
% y: vector of state variables
% data: struct of general model parameters
% nStrata: number of strata within y
% dis: struct of pathogen parameters
% i: current mandate index, corresponding to a configuration in
% get_strategy_design.m
% p2: struct of p2 intervention parameters
%
% value: vector of values associated with events, where we are interested
% to know if one crosses 0
% isterminal: whether to interrupt the solver
% direction: whether the sign change is detected going in one direction
% only

function [value,isterminal,direction] = reactive_closures(t,y,data,nStrata,dis,i,p2)
    
    y_mat = reshape(y,nStrata,[]);
    compindex = data.compindex;
    S_index = compindex.S_index;
    I_index = compindex.I_index;
    H_index = compindex.H_index;
    
    S    = y_mat(:,S_index(1));
    H    = y_mat(:,H_index(1));
    Hv1    = y_mat(:,H_index(2));
    Hv2    = y_mat(:,H_index(3));
    Sn   = y_mat(:,S_index(2));
    S01   = y_mat(:,S_index(3));
    Sv1   = y_mat(:,S_index(4));
    Sv2   = y_mat(:,S_index(5));
    S02   = y_mat(:,S_index(6));
    S12   = y_mat(:,S_index(7));
    
    Is   = y_mat(:,I_index(2));
    Isv1   = y_mat(:,I_index(4));
    Isv2   = y_mat(:,I_index(6));
    
    occ   = max(1,sum(H+Hv1+Hv2)); 
    
    dis2 = update_hosp_dis_parameters(occ, p2, dis);
    
    g3    = dis2.g3;
    mu    = dis2.mu;

    %% get waning for R values and H derivatives
    
    dis2 = update_vax_dis_parameters(dis2, S, Sn, compindex, y_mat);

    Hdot   =         dis2.h.*Is   -(g3+mu).*H;
    Hv1dot = dis2.h_v1.*Isv1 -(g3+mu).*Hv1;
    Hv2dot = dis2.h_v2.*Isv2 -(g3+mu).*Hv2;
    
    occdot = sum(Hdot+Hv1dot+Hv2dot);
    r      = occdot/occ;
    Tcap   = t + log(p2.Hmax/occ)/r - 7;
    
    %% isolating
    sumN = sum(data.NNs);
    [p3, p4] = fraction_averted_self_isolating(sum(sum(y_mat(:,I_index))), sumN, p2, t, i);
          
    %% distancing
    ddk    = 10^6*sum(mu.*(H + Hv1 + Hv2))/sumN;
%     betamod = social_distancing(data.sd_baseline, data.sd_death_coef, data.sd_mandate_coef,ddk, data.rel_stringency(i));
     
    %% Event 1: Response Time
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0) ;
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: Early Lockdown
    
    value(2)     = - abs((i-2)*(i-4)) + min(t-(data.tvec(end-1)+0.1),0) + min(t-Tcap,0);
    if value(2)==0
%         disp([t r])
    end
    direction(2) = 1;
    if r>0.025
        isterminal(2) = 1;
    else
        isterminal(2) = 0;
    end
    
    %% Event 3: Late Lockdown
    
    value(3)      = - abs((i-1)*(i-2)*(i-4)) + min(t-(data.tvec(end-1)+0.1),0) + min(occ-0.95*p2.Hmax,0);
    direction(3)  = 1;
    isterminal(3) = 1;
    
    %% Event 4: Reopening
    
    value(4)      = abs(i-3) + abs(min(t-(data.tvec(end-1)+7),0)) + max(0,occ-p2.hosp_release_trigger);
    direction(4)  = -1;
    isterminal(4) = 1;
    
    %% Event 5: End
    
    % not in hard lockdown: i = 1, 2 or 4; ivals = 0
    ivals = -abs((i-1)*(i-2)*(i-4));
    % have passed penultimate timepoint: tval = 0
    tval = min(t-(data.tvec(end-1)+7),0) + min(t-7-max(p2.tpoints),0);
    % (low growth rate OR occupancy is low) AND have reached end of vaccine
    % rollout: otherval = 0
    otherval = -abs(min(0.025-r,0)*max(0,occ-p2.hosp_release_trigger));
    R2flag = otherval + ivals + tval;
    % only check open economy if current R<1
    R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
    if ivals==0 && tval==0 && R_est<0.95
        if otherval~=0
        % only compute R if R2flag is not already 0 and ivals and tval
        % conditions are both met
            Rt2 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, ddk, data, 5);
            R2flag = min(0.95-Rt2,0);
        end
    end
    
    value(5)      = R2flag; 
    %measures can be removed if (not in hard lockdown) and ((Rt<1) or (after end of vaccination campaign and below 25% occupancy or low growth rate))
    direction(5)  = 0;
    isterminal(5) = 1;
    
    %% Event 6: importation
    value(6)      =  min(t-data.t_import,0);
    direction(6)  = 1;
    isterminal(6) = 1;
    
    %% Event 7: end
    % mitigation is over
    ival = -abs((i-5));
    % t is greater than the end of the vaccine rollout: tval = 0
    tval = min(t-(max(p2.tpoints)+7),0) + min(t-(data.tvec(end-1)+7),0);
    % hval: no patients
    sumH = sum(H + Hv1 + Hv2);
    hval = min(p2.hosp_final_threshold-sumH,0);
    R3flag = tval + hval + ival;
    if tval==0 && hval==0 && ival==0
%         Rt3 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,dis.beta,p3,p4, ddk, data, 5);
        Rthresh = exp(dis.generation_time*log(2) / p2.final_doubling_time_threshold); % R value for which doubling time is 30 days
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        R3flag = min(Rthresh - R_est,0);
%         disp([Rt3, R_est])
%         disp([t/1000 Rt3 sumH])
%         if t<600 
%             disp([t/100 hval tval r Rt3 Rthresh ])
%         end
    end
    value(7)      =  R3flag;
    direction(7)  = 0;
    isterminal(7) = 1;
    
end
