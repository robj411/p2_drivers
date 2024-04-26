% function implementing elimination mitigation strategy
%
% t: current time
% y: vector of state variables
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

function [value,isterminal,direction] = elimination(t,y,data,nStrata,dis,i,p2)
    
    y_mat = reshape(y,nStrata,[]);
    compindex = data.compindex;
    sumN = sum(data.NNs);
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
    
    hospital_occupancy = H + Hv1 + Hv2;
    sumH = sum(hospital_occupancy);
    occ   = max(1,sumH);     
    
    %% isolating
    [p3, p4] = fraction_averted_self_isolating(sum(sum(y_mat(:,I_index))), sumN, p2, t, i);
    
    %% correct ph for previous infection
    
    dis2 = update_vax_dis_parameters(dis, S, Sn, compindex, y_mat);
     % correct for hosp occupancy
    dis2 = update_hosp_dis_parameters(occ, p2, dis2);

    g3    = dis2.g3;
    mu    = dis2.mu;
    
    Hdot   =         dis2.h.*Is   -(g3+mu).*H;
    Hv1dot = dis2.h_v1.*Isv1 -(g3+mu).*Hv1;
    Hv2dot = dis2.h_v2.*Isv2 -(g3+mu).*Hv2;
    
    occdot = sum(Hdot+Hv1dot+Hv2dot);
    
    %% distancing
    ddk    = 10^6*sum(dis2.mu.*hospital_occupancy)/sumN;    
    
    %% Event 1: Early Lockdown
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0);
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: Late Lockdown
    
    value(2)      = - abs(i-1) + min(occ-0.95*p2.Hmax,0);
    direction(2)  = 1;
    isterminal(2) = 1;
    
    %% Event 3: Reopening
    
    % Rt values for 3rd config
    R1flag3 = -1;
    minttvec3 = min(t-(data.tvec(end-1)+7),0);
    R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
    if i==2 && minttvec3==0  && R_est<.95
        Rt1 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, ddk, data, 3);
        R1flag3 = min(0.95-Rt1,0);
    end
    
    value(3)      = - abs(i-2) + minttvec3 + R1flag3;
    direction(3)  = 0;
    isterminal(3) = 1;
    
    %% Event 4: Relockdown
    
%     R1flag4 = -1;
    minttvec4 = min(t-(data.tvec(end-1)+7),0);
%     if  (i==3  && minttvec4==0)      
%         Rt1 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,dis.beta,p3,p4, ddk, data, 3);
    R1flag4 = min(R_est-1.2000,0);
%     end
    
    value(4)      = - abs(i-3) + minttvec4 + R1flag4;
    direction(4)  = 0;
    isterminal(4) = 1;
    
    %% Event 5: End
    % i is in 1:4: ival = 0
    ival = -abs((i-1)*(i-2)*(i-3)*(i-4));
    % t is greater than the penultimate timepoint: tval = 0
    tval = min(t-(data.tvec(end-1)+0.1),0);
    % t is greater than the end of the vaccine rollout: otherval = 0
    otherval = min(t-max(p2.tpoints),0);
    R2flag = otherval + ival + tval;
    if ival==0 && tval==0
        if otherval~=0 && R_est<.95
            Rt2 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, ddk, data, 5);
            R2flag = min(0.95-Rt2,0);
%             disp([t Rt2])
        end
    end
    
    value(5)      = R2flag; % min(t-max(p2.tpoints),0)*min(1.00-Rt2,0); 
    %measures can be removed at any stage if (Rt<1) or (after end of vaccination campaign)
    direction(5)  = 0;
    isterminal(5) = 1;
    
    %% Event 6: importation
    value(6)      =  min(t-data.t_import,0);
    direction(6)  = 1;
    isterminal(6) = 1;
    
    %% Event 7: end simulation
    % i is 5
    ival = -abs((i-5));
    % t is greater than the end of the vaccine rollout and a week since the last changepoint (which was end mitigation): tval = 0
    tval = min(t-(max([p2.tpoints data.tvec(end-1)])+7),0);
    % tlong: one month since end mitigation
    tlong = min(t-(data.tvec(end-1)+365),0);
    % patient numbers declining
    hdotval = min(-occdot,0);
    % hval: few patients
    hval = min(p2.hosp_final_threshold - sumH,0);
    % either: more than one year has passed
    % or: H is low and coming down
    R3flag = ival + (hval + hdotval + ival + tval);
    if ival==0 && tlong==0 && R3flag < 0
%         Rt3 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,dis.beta,p3,p4, ddk, data, 5);
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        doubling_time = log(2)*dis.generation_time/max(R_est-1,1e-5);
        R3flag = min(doubling_time - p2.final_doubling_time_threshold,0);
%         disp([t/1000 Rt3 betamod p3 p4 sumH ])
    end
    
    value(7)      =  R3flag;
    direction(7)  = 1;
    isterminal(7) = 1;
    
end 


