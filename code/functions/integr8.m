% function that calls the ode solver
%
% data: struct containing fixed values
% contact_matrix: 49 by 49 contact matrix used this time
% i: current event state
% t0: start time of this time segment
% tend: maximum time for projection
% dis: struct of disease parameters
% y0: initial conditions for solver
% p2: struct of response parameters
%
% tout: vector of times
% Iclass: array of I values
% Iaclass: array of asymptomatic I values
% Isclass: array of symptomatic I values
% Hclass: array of I values
% Dclass: array of I values
% p3: vector of asymptomatic self isolating
% p4: vector of symptomatic self isolating
% betamod: vector of social distancing values
% y0new: end conditions / initial conditions for next time
% inext: next event state
% still_susc: vector of number still susceptible
% data: struct containing fixed values, with potential for updated state-6
% variables

function [tout,Iclass,Iaclass,Isclass,Hclass,Dclass,p3,p4,betamod,y0new,inext,still_susc,data]=...
          integr8(data,contact_matrix,i,t0,tend,dis,y0,p2)
      
    %% copy over only required items
    
    rundata = struct;
    rundata.NNs = data.NNs; 
    rundata.yll = data.yll; 
    rundata.Dvec = data.Dvec;
    rundata.tvec = data.tvec;
    rundata.compindex = data.compindex; 
    rundata.t_import = data.t_import; 
    rundata.imand = data.imand; 
    rundata.seedsize = data.seedsize; 
    rundata.sd_baseline = data.sd_baseline;
    rundata.sd_death_coef = data.sd_death_coef;
    rundata.sd_mandate_coef = data.sd_mandate_coef;
    rundata.sd_decay_rate = data.sd_decay_rate;
    
    NN0 = data.NNs; 
    %% CALL

    compindex = data.compindex;
    S_index = compindex.S_index;
    E_index = compindex.E_index;
    I_index = compindex.I_index;
    H_index = compindex.H_index;
    D_index = compindex.D_index;
    
    nStrata = size(NN0,1);
    fun  = @(t,y)ODEs(rundata,contact_matrix,i,t,dis,y,p2);
    strategy = data.strategy;

    if strcmp(strategy,'Elimination')
        options = odeset('Events',@(t,y)mitigate(t,y,data,nStrata,dis,i,p2,'Elimination'));
    elseif strcmp(strategy,'Economic Closures') || strcmp(strategy,'School Closures')
        options = odeset('Events',@(t,y)mitigate(t,y,data,nStrata,dis,i,p2,'Reactive closures'));
    elseif strcmp(strategy,'No Closures')
        options = odeset('Events',@(t,y)mitigate(t,y,data,nStrata,dis,i,p2,'Unmitigated'));
    end

    % solve ODEs
    [tout,yout,~,~,ie] = ode45(fun,[t0 tend],y0,options);

    % process event
    y0new     = yout(end,:)'; 
    y_mat    = reshape(y0new,nStrata,[]);
    
    % ie is the index ("value") returned
    % inext is the value it maps to (from get_strategy_design)
    % disp([max(tout) i ie' ie'])
    if tout(end)<tend
        if ie <= 6
            inext = data.inext(ie(end));
%             disp(data.rel_stringency(inext))
        elseif ie == 8 % importation event
            % keep the same i
            inext = i;
            % move 5 people from S to E
            current_S = y_mat(:,S_index(1));
            imported = 5/sum(current_S)*current_S;
            new_E = y_mat(:,E_index(1)) + imported;
            new_S = current_S - imported;
            y_mat(:,S_index(1)) = new_S;
            y_mat(:,E_index(1)) = new_E;
            y0new = reshape(y_mat,[],1);
        else %end
            inext = 0;
        end
    else
        inext = NaN;
    end
    

    %% OUTPUT VARIABLES

    indices = 1:nStrata;
    Ia   = yout(:,(I_index(1)-1)*nStrata + indices);
    Is   = yout(:,(I_index(2)-1)*nStrata + indices);
    Iav1   = yout(:,(I_index(3)-1)*nStrata + indices);
    Isv1   = yout(:,(I_index(4)-1)*nStrata + indices);
    Iav2   = yout(:,(I_index(5)-1)*nStrata + indices);
    Isv2   = yout(:,(I_index(6)-1)*nStrata + indices);
    H     = yout(:,(H_index(1)-1)*nStrata + indices);
    Hv1     = yout(:,(H_index(2)-1)*nStrata + indices);
    Hv2     = yout(:,(H_index(3)-1)*nStrata + indices);
    Died     = yout(:,(D_index(1)-1)*nStrata + indices);
    still_susc = sum(yout(:,(repelem(S_index([1,3:end]),1,length(indices))-1)*nStrata + repmat(indices,1,length(S_index)-1)),2);

    Iclass   = Ia + Is + Iav1 + Isv1 + Iav2 + Isv2; 
    Hclass   = H + Hv1 + Hv2; 
    Dclass   = Died;

    %% TIME-DEPENDENT PARAMETERS
    
    occ   = max(1,sum(Hclass,2));
    mu = zeros(size(Hclass));
    for ii = 1:size(mu,1)
        dis2 = update_hosp_dis_parameters(occ(ii), p2, dis);
        mu(ii,:)    = dis2.mu;
    end
    
    sumNN0 = sum(NN0);
    ddk = 10^6*sum(mu.*Hclass,2)/sumNN0;
    foi0 = data.basic_foi;
    foi1 = contact_matrix*(1./NN0);
    sd_so_far = ((foi1'*NN0+1e-10)./(foi0'*NN0+1e-10));
    betamod = betamod_wrapped(ddk, data, i, 1-sd_so_far, tout);
    
    [p3, p4] = fraction_averted_self_isolating(sum(Iclass,2), sumNN0, p2, tout, i);
    
    Iaclass = zeros(size(Ia,1),size(Ia,2),3);
    Isclass = zeros(size(Ia,1),size(Ia,2),3);
    Iaclass(:,:,1) = Ia; 
    Iaclass(:,:,2) = Iav1; 
    Iaclass(:,:,3) = Iav2; 
    Isclass(:,:,1) = Is; 
    Isclass(:,:,2) = Isv1; 
    Isclass(:,:,3) = Isv2; 
    
    %% get exit status variables

    if inext==5 
        
        S    = y_mat(:,compindex.S_index(1));
        Sn   = y_mat(:,compindex.S_index(2));
        S01   = y_mat(:,compindex.S_index(3));
        Sv1   = y_mat(:,compindex.S_index(4));
        Sv2   = y_mat(:,compindex.S_index(5));
        S02   = y_mat(:,compindex.S_index(6));
        S12   = y_mat(:,compindex.S_index(7));

        dis2 = update_vax_dis_parameters(dis2, S, Sn, compindex, y_mat);
    
        t = tout(end);
        closure = 1 - data.workerConfigMat(:,i);   
        trial_vals = 0:0.25:1;
        Rthresh = log(2)*dis.generation_time/p2.final_doubling_time_threshold + 1;
        Rs = zeros(1, length(trial_vals));
        for j =1:length(trial_vals)
            trial_val = trial_vals(j);
            openness = 1 - trial_val*closure;
            Rs(j) = get_trial_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3(end),p4(end), ddk(end), data,openness,data.hw(2,:),5, t);
            if Rs(j) < Rthresh
                break;
            end
        end
            
        if min(Rs>Rthresh)
            index = length(trial_vals);
        else
            index = find(Rs<Rthresh,1);
        end
        if data.exittype==2
            index=1;
        end
        openness = 1 - trial_vals(index)*closure;
        data.workerConfigMat(:,6) = openness;
        data.hw(6,:) = data.hw(2,:); % this can go into setup
        Dtemp = p2MakeDs(data,data.NNs,openness,data.hw(6,:));
        % store matrix in list
        data.Dvec(:,:,6) = Dtemp;
        if data.exittype > 0 % if not in exit mode, save data for restart
            inext = 6;
        end
    end 

end

%%

function [f] = ODEs(data,contact_matrix,i,t,dis,y,p2)

    NN0 = data.NNs;
    nStrata = size(NN0,1);
    y_mat = reshape(y,nStrata,[]); 

    %% initial conditions
    
    compindex = data.compindex;
    S_index = compindex.S_index;
    E_index = compindex.E_index;
    I_index = compindex.I_index;
    H_index = compindex.H_index;
    R_index = compindex.R_index;
    D_index = compindex.D_index;
    V_index = compindex.V_index;

    S =      y_mat(:,S_index(1));
    E =      y_mat(:,E_index(1));
    Ia =    y_mat(:,I_index(1));
    Is =    y_mat(:,I_index(2));
    H =      y_mat(:,H_index(1));
    R =      y_mat(:,R_index(1));
    Sn =     y_mat(:,S_index(2));
    S01 =   y_mat(:,S_index(3));
    Sv1 =    y_mat(:,S_index(4));
    Sv2 =    y_mat(:,S_index(5));
    S02 =   y_mat(:,S_index(6));
    S12 =   y_mat(:,S_index(7));
    Ev1 =    y_mat(:,E_index(2));
    Iav1 =    y_mat(:,I_index(3));
    Isv1 =    y_mat(:,I_index(4));
    Hv1 =    y_mat(:,H_index(2));
    Rv1 =    y_mat(:,R_index(2));
    Ev2 =    y_mat(:,E_index(3));
    Iav2 =    y_mat(:,I_index(5));
    Isv2 =    y_mat(:,I_index(6));
    Hv2 =    y_mat(:,H_index(3));
    Rv2 =    y_mat(:,R_index(3));
    DE =    y_mat(:,D_index(1));
    
    %% SELF-ISOLATION

    [p3, p4] = fraction_averted_self_isolating(sum(Ia+Is + Iav1+Isv1 + Iav2+Isv2), sum(NN0), p2, t, i);
    
    Ina = (1-p3) .* Ia;
    Inav1 = (1-p3) .* Iav1;
    Inav2 = (1-p3) .* Iav2;
    Ins = (1-p4) .* Is;
    Insv1 = (1-p4) .* Isv1;
    Insv2 = (1-p4) .* Isv2;

    %% TIME-DEPENDENT DISEASE PARAMETERS

    % unchanged
    sig1 = dis.sig1;
    sig2 = dis.sig2;
    g1   = dis.g1;
    nu   = dis.nu;
    hrv1 = dis.hrv1;
    
    % correct for vaccine
    [dis2, V, B, vaccine_pp, booster_pp] = update_vax_dis_parameters(dis, S, Sn, compindex, y_mat);
    
    % correct for hosp occupancy
    hospital_occupancy = H + Hv1 + Hv2;
    dis2 = update_hosp_dis_parameters(max(1,sum(hospital_occupancy)), p2, dis2);
    
    scv1 = dis2.scv1;
    scv2 = dis2.scv2;
    h    = dis2.h;
    g2   = dis2.g2;
    g2_v1 = dis2.g2_v1;
    g2_v2 = dis2.g2_v2;
    h_v1  = dis2.h_v1;
    h_v2  = dis2.h_v2;
    g3 = dis2.g3;
    mu = dis2.mu;

    %% FOI
    
    foi = get_foi(dis2, hospital_occupancy, data, i,...
        Ina,Ins,Inav1,Insv1,Inav2,Insv2,contact_matrix, t);
    

    %% VACCINATION
    
    [v1rates, v1rater, v2rates, v2rater, v12rates, v12rater] = ...
        get_vax_rates(p2, t, nStrata, R,S,DE, Rv1,Sv1);

    %% EQUATIONS

    Sdot=      -v1rates - v2rates      -S.*foi  +nu.*R   ; 
    S01dot=    v1rates               - hrv1*S01   -S01.*foi;
    Sv1dot=    -v12rates +              hrv1*S01   -Sv1.*(1-scv1).*foi  + nu.*Rv1;  
    Sndot=      -Sn.*foi    -(v1rates+v2rates).*Sn./S;  
    S02dot=    v2rates               - hrv1*S02   - S02.*foi;
    S12dot=    v12rates              - hrv1*S12   - S12.*(1-scv1).*foi;
    Sv2dot=     hrv1*S02 + hrv1*S12     -Sv2.*(1-scv2).*foi  + nu.*Rv2;  
    
    Edot=             (S + S01 + S02).*foi  - (sig1+sig2).*E ; 
    Ev1dot=     (S12 + Sv1).*(1-scv1).*foi  - (sig1+sig2).*Ev1;
    Ev2dot=             Sv2.*(1-scv2).*foi  - (sig1+sig2).*Ev2;

    [Iadot, Isdot, Iav1dot, Isv1dot, Iav2dot, Isv2dot] =  get_Idot(dis2, compindex, y_mat);

    Hdot=       h.*Is         -(g3+mu).*H;
    Hv1dot=     h_v1.*Isv1    -(g3+mu).*Hv1;
    Hv2dot=     h_v2.*Isv2    -(g3+mu).*Hv2;

    Rdot=       g1.*Ia       +g2.*Is          +g3.*H      - v1rater  - v2rater  - nu.*R;
    Rv1dot=     g1.*Iav1     +g2_v1.*Isv1     +g3.*Hv1    + v1rater  - v12rater - nu.*Rv1;   
    Rv2dot=     g1.*Iav2     +g2_v2.*Isv2     +g3.*Hv2    + v2rater  + v12rater - nu.*Rv2;   

    DEdot=      mu.*H     + mu.*Hv1    + mu.*Hv2   ;  
    
    Vdot =  - dis.nuv1*V - ... % decay
        (mu.*Hv1 + v12rater + hrv1*S12).*vaccine_pp + ... % exit
        v1rater + hrv1*S01; % entry
    Bdot =  - dis.nuv2*B - ...% decay
        mu.*Hv2.*booster_pp + ... % exit
        v2rater + v12rater + hrv1*S12 + hrv1*S02; % entry

    %% OUTPUT
    
    f_mat = zeros(size(y_mat));
    f_mat(:,S_index(1)) = Sdot;
    f_mat(:,E_index(1)) = Edot;
    f_mat(:,I_index(1)) = Iadot;
    f_mat(:,I_index(2)) = Isdot;
    f_mat(:,H_index(1)) = Hdot;
    f_mat(:,R_index(1)) = Rdot;
    f_mat(:,S_index(2)) = Sndot;
    f_mat(:,S_index(3)) = S01dot;
    f_mat(:,S_index(4)) = Sv1dot;
    f_mat(:,S_index(5)) = Sv2dot;
    f_mat(:,S_index(6)) = S02dot;
    f_mat(:,S_index(7)) = S12dot;
    f_mat(:,E_index(2)) = Ev1dot;
    f_mat(:,I_index(3)) = Iav1dot;
    f_mat(:,I_index(4)) = Isv1dot;
    f_mat(:,H_index(2)) = Hv1dot;
    f_mat(:,R_index(2)) = Rv1dot;
    f_mat(:,E_index(3)) = Ev2dot;
    f_mat(:,I_index(5)) = Iav2dot;
    f_mat(:,I_index(6)) = Isv2dot;
    f_mat(:,H_index(3)) = Hv2dot;
    f_mat(:,R_index(3)) = Rv2dot;
    f_mat(:,D_index(1)) = max(DEdot,0);
    f_mat(:,V_index(1)) = Vdot;
    f_mat(:,V_index(2)) = Bdot;
    
    f = reshape(f_mat,[],1);
    eps10 = eps*1e12;
    f(y<eps10) = max(0,f(y<eps10)); %%! exit wave was lost

end


