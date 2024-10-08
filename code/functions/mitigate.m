% function implementing strategies
%
% t: current time
% y: vector of state variables
% data: struct of general model parameters
% nStrata: number of strata within y
% dis: struct of pathogen parameters
% i: current mandate index, corresponding to a configuration in
% get_strategy_design.m
% p2: struct of p2 intervention parameters
% strategy: string, which strategy type
%
% value: vector of values associated with events, where we are interested
% to know if one crosses 0
% isterminal: whether to interrupt the solver
% direction: whether the sign change is detected going in one direction
% only

function [value,isterminal,direction] = mitigate(t,y,data,nStrata,dis,i,p2,strategy)
    
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
    still_susc = sum(S+S01+Sv1+Sv2+S02+S12);
    
    Is   = y_mat(:,I_index(2));
    Isv1   = y_mat(:,I_index(4));
    Isv2   = y_mat(:,I_index(6));
    
    Isum = sum(y_mat(:,I_index)')';
    
    hospital_occupancy = H + Hv1 + Hv2;
    sumH = sum(hospital_occupancy);
    occ   = max(1e-10,sumH); 
    
    %% isolating
    
    [p3, p4] = fraction_averted_self_isolating(sum(sum(y_mat(:,I_index))), sumN, p2, t, i);
    
    %% correct for vaccine waning and previous infection
    
    dis2 = update_vax_dis_parameters(dis, S, Sn, compindex, y_mat);
    
    %% correct for hosp occupancy
    
    dis2 = update_hosp_dis_parameters(occ, p2, dis2, t);

    g3    = dis2.g3;
    mu    = dis2.mu;
    
    Hdot   =         dis2.h.*Is   -(g3+mu).*H;
    Hv1dot = dis2.h_v1.*Isv1 -(g3+mu).*Hv1;
    Hv2dot = dis2.h_v2.*Isv2 -(g3+mu).*Hv2;
    
    occdot = sum(Hdot+Hv1dot+Hv2dot);
    
    % variables for reactive closures
    r      = occdot/max(occ,1);
    
    
    %% uncosted transmission reduction
    
    deaths_per_mill    = 10^6*sum(mu.*hospital_occupancy)/sumN;    
    
    %% Event 1: Response time
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0);
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: 
    if strcmp(strategy,"Unmitigated") || strcmp(strategy,"Elimination")
    	% Late Lockdown
	    value(2)      = - abs(i-1) + min(occ-0.95*p2.Hmax,0);
	    direction(2)  = 1;
	    isterminal(2) = 1;
    else
    	% Early Lockdown; reactive closure
        Rflag2 = - abs((i-2)*(i-4)) + min(t-(data.tvec(end-1)+0.1),0) ;
        if i==2 && Rflag2 == 0
            R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
            Rt = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, deaths_per_mill, data, 2, t, Isum);
            growth_rate = (Rt-1)/dis.generation_time;
            time_to_breach   = log(p2.Hmax/occ)/growth_rate;
%             disp([t r Rt growth_rate occ time_to_breach log(p2.Hmax/occ)])
            
            Rflag2 = min(4-time_to_breach,0);
        end
	    value(2)     =  Rflag2;
	    direction(2) = 1;
	    if r>0.025
            isterminal(2) = 1;
	    else
            isterminal(2) = 0;
        end
    end
	    
    
    %% Event 3: 
    % for unmitigated:
    value(3:4)      = 0; 
    direction(3:4)  = 0;
    isterminal(3:4) = 0;
    
    if strcmp(strategy,"Elimination") % Reopening
	    % Rt values for 3rd config
	    R1flag3 = -1;
	    minttvec3 = min(t-(data.tvec(end-1)+7),0);
	    R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        
	    if i==2 && minttvec3==0  && R_est<1
            Rt1 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, deaths_per_mill, data, 3, t, Isum);
            R1flag3 = min(1-Rt1,0);
	    end
	    
	    value(3)      = - abs(i-2) + minttvec3 + R1flag3;
	    direction(3)  = 0;
	    isterminal(3) = 1;
    elseif strcmp(strategy,"Reactive closures") % Late Lockdown
    
	    value(3)      = - abs((i-1)*(i-2)*(i-4)) + min(t-(data.tvec(end-1)+0.1),0) + min(occ-0.95*p2.Hmax,0);
	    direction(3)  = 1;
	    isterminal(3) = 1;
    end
    
    %% Event 4: 
    
    if strcmp(strategy,"Elimination") % Relockdown
        % it's been a week, or response time was less than a week ago:
	    minttvec4 = - min(t-(data.tvec(end-1)+7),0) * min(p2.Tres+7 - t,0) + min(t-(data.tvec(end-1)+.1),0);
	    R1flag4 = min(R_est-1.2000,0);
%         if i==3 && minttvec4==0  && t<50
%             Rt4 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, deaths_per_mill, data, 3, t, Isum);
%             R4flag = min(Rt4-1.2,0);
%             disp([t Rt4 R_est R1flag4])
% 	    end
% 	    
	    value(4)      = - abs(i-3) + minttvec4 + R1flag4;
	    direction(4)  = 0;
	    isterminal(4) = 1;
    elseif strcmp(strategy,"Reactive closures") % Reopening
%         if t<120
%             disp([i-3,t,(p2.hosp_release_trigger - occ)/1000, r])
%         end
	    value(4)      = - abs(i-3) + min(t-(data.tvec(end-1)+7),0) + min(0, p2.hosp_release_trigger - occ) + min(-r,0);
	    direction(4)  = 0;
	    isterminal(4) = 1;
    end
    
    %% Event 5: End mitigation
    
    if strcmp(strategy,"Elimination") || strcmp(strategy,"Unmitigated")   
        %measures can be removed at any stage if (Rt<1) or (after end of vaccination campaign)        
	    % i is in 1:4: ival = 0
	    ivals = -abs((i-1)*(i-2)*(i-3)*(i-4));
	    % t is greater than the penultimate timepoint: tval = 0
	    tval = min(t-(data.tvec(end-1)+7),0) + min(t-7-max(p2.tpoints),0);
	    % have reached end of vaccine rollout: otherval = 0
	    otherval = min(t-max(p2.tpoints+7),0);
    elseif strcmp(strategy,"Reactive closures")
	    %measures can be removed if (not in hard lockdown) and ((Rt<1) or (after end of vaccination campaign and below 25% occupancy or low growth rate))
	    % not in hard lockdown: i = 1, 2 or 4; ivals = 0
	    ivals = -abs((i-1)*(i-2)*(i-4)*(i-3));
	    % have passed penultimate timepoint: tval = 0
	    tval = min(t-(data.tvec(end-1)+7),0) + min(t-7-max(p2.tpoints),0);
	    % (low growth rate OR occupancy is low) AND have reached end of vaccine
	    % rollout: otherval = 0
	    otherval = -abs(min(0.025-r,0)*max(0,occ-p2.hosp_release_trigger));
    end
    
    R2flag = otherval + ivals + tval;
    if ivals==0 && tval==0 
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        if otherval~=0 && R_est<1
            % only compute R if R2flag is not already 0 and ivals and tval
            % conditions are both met
            Rt2 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, deaths_per_mill, data, 5, t, Isum);
            R2flag = min(1-Rt2,0);
        end
    end
	    
    value(5)      = R2flag; 
    direction(5)  = 0;
    isterminal(5) = 1;

    %% Event 6: end exit wave

    % we are in the exit wave phase
    ival = -abs(i-6);
    % t is one week greater than the last changepoint (which was end mitigation): tval = 0
    tval = min(t-(data.tvec(end-1)+30),0);
    % current R is less than 1
%     R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
    % tlong: one year since end mitigation
    tlong = min(t-(data.tvec(end-1)+365),0);
    R6flag = tlong + ival + tval;
    
    if ival==0 && tval==0 && R6flag ~= 0
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4);
        Rt6 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, deaths_per_mill, data, 6, t, Isum);
        R6flag = min(1.0 - R_est, 0);
%         disp([t/1000 R_est Rt6])
%         R6flag = min(doubling_time - p2.final_doubling_time_threshold,0);
    end
    
    value(6)      =  R6flag;
    direction(6)  = 1;
    isterminal(6) = 1;
    
    %% Event 7: end simulation
    % mitigation is over
    ival = -abs(i-5);
    % t is one month greater than the end of the vaccine rollout and the last changepoint (which was end mitigation): tval = 0
    tval = min(t-(data.tvec(end-1)+30),0);
    % tlong: one year since end mitigation
    tlong = min(t-(data.tvec(end-1)+365),0);
    % patient numbers declining
    hdotval = min(-occdot,0);
    % hval: few patients
    hval = min(p2.hosp_final_threshold - sumH,0);
    % either: more than one year has passed
    % or: H is low and coming down
    R7flag = ival + hval + tval + hdotval;
%     R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
%     disp([t/1000 hdotval still_susc/50000000 R_est])
%     disp([ival tlong R7flag still_susc/50000000])
    if ival==0 && tlong==0 && R7flag ~= 0 %&& still_susc/50000000 < .85
%         Rt3 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,dis.beta,p3,p4, deaths_per_mill, data, 5, t);
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        doubling_time = log(2)*dis.generation_time/max(R_est-1,1e-5);
        R7flag = min(doubling_time - p2.final_doubling_time_threshold,0);
%         disp(doubling_time)
    end
    
    
    value(7)      =  R7flag;
    direction(7)  = 1;
    isterminal(7) = 1;
    
    %% Event 8: importation
    value(8)      =  min(t-data.t_import,0);
    direction(8)  = 1;
    isterminal(8) = 1;
    
end


