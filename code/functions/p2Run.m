% run model
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% strategy: string specifying which mitigation strategy is followed
% p2: struct of p2 intervention parameters
%
% data: struct of general model parameters
% returnobject: struct of epi outcomes 

function [data, returnobject] = p2Run(data, dis, strategy, p2)

    %% get configurations
    
    [data, configuration] = get_strategy_design(data, strategy, p2);
    
    %% get contact matrices
    
    nSectors = data.nSectors;
    nStrata = length(data.NNs);
    nConfigs = length(configuration)/data.nSectors;

    NNbar = data.NNs;
    configMat = reshape(configuration,nSectors,nConfigs);
    workerConfigMat = configMat;
    zs = zeros(size(data.NNs));
%     NNvec = repmat(NNbar,1,nConfigs);
%     NNvec                = repmat(NNbar(1:nSectors),1,int).*workerConfigMat;
%     NNworkSum            = sum(NNvec,1);
%     NNvec(nSectors+1:nStrata,:)     = repmat(NNbar(nSectors+1:nStrata),1,int);
%     NNvec(nSectors+adInd,:)    = sum(NNbar([1:nSectors,nSectors+adInd]))-NNworkSum;
%     data.NNvec = NNvec;
    
    % get low-contact matrices
    config_min_econ = data.x_econ(:,2)*.95;
    Dminecon = p2MakeDs(data,NNbar,config_min_econ,data.wfh(2,:));
    candidate_infectees_econ = get_candidate_infectees(nStrata, dis, NNbar, zs, zs, 0, 0, NNbar, Dminecon);
    config_min_sch = data.x_schc(:,2)*.95;
    Dminschc = p2MakeDs(data,NNbar,config_min_sch,data.wfh(2,:));
    candidate_infectees_schc = get_candidate_infectees(nStrata, dis, NNbar, zs, zs, 0, 0, NNbar, Dminschc);
    candidate_infectees_min = min(candidate_infectees_schc,candidate_infectees_econ);
   
    candidate_infectees = zeros(nConfigs,1);
    Dvec = zeros(nStrata,nStrata,nConfigs);
    for i = 1:nConfigs
        % get matrix per configuration
        Dtemp = p2MakeDs(data,NNbar,configMat(:,i),data.hw(i,:));
        % store matrix in list
        Dvec(:,:,i) = Dtemp;

        candidate_infectees(i) = get_candidate_infectees(nStrata,dis,NNbar, zs, zs, 0, 0, NNbar, Dtemp);
%         candidate_infectees(i) = get_R(nStrata, dis, NNbar, zs, zs,...
%            NNbar, Dtemp, 1, 1, 0, 0);
    end
    candidate_infectees_max = candidate_infectees(1);
    
    % store amount by which configurations reduce contacts
%     data.rel_mobility = (candidate_infectees_max - candidate_infectees)./(candidate_infectees_max-candidate_infectees_min);
    rel_mobility = (candidate_infectees)./(candidate_infectees_max);
    rel_mobility_min = candidate_infectees_min/candidate_infectees_max;
    data.rel_stringency = (1-rel_mobility) / max(1-[rel_mobility_min; rel_mobility]);
    data.Dvec = Dvec;
    data.strategy = strategy;

    data.basic_foi = Dvec(:,:,1)*(1./NNbar);
    
    returnobject = p2SimVax(data, dis, workerConfigMat, p2);

end

%%

function returnobject = p2SimVax(data, dis, workerConfigMat, p2)    
    
    
    %% PARAMETERS
    nStrata = size(data.NNs,1);
    nSectors = data.nSectors;
    NNbar = data.NNs;
    sumWorkingAge = sum(NNbar([1:nSectors,nSectors+3]));
    compindex = data.compindex;
    nODEs = max(struct2array(compindex));
    S0 = NNbar;
    Dvec = data.Dvec;
    
    % initial conditions
    t0 = data.tvec(1);
    y0_mat = zeros(nStrata,nODEs);
    y0_mat(:,compindex.S_index(1)) = S0;
    y0_mat(:,compindex.S_index(2)) = S0;
    y0 = reshape(y0_mat,[],1);

    % store outputs
    zn = zeros(1,nStrata);
    zn3 = zeros(1,nStrata,3);
    Tout       = t0;
    Iout       = zn;
    Iaout     = zn3;
    Isout     = zn3;
    Hout       = zn;
    Dout       = zn;
    Wout       = [];
    hwout      = [];
    p3out    = 0;
    p4out    = 0;
    betamodout = 1;
    Sout       = sum(S0);
    rout       = 0;

    %% LOOP

    i = 1;
    isequence = [];

    tend = data.tvec(end);

    while Tout(end)<tend 

        Wit               = workerConfigMat(:,i);    
        NNfeed            = NNbar;
        NNfeed(NNfeed==0) = 1;
        D                 = Dvec(:,:,i);

        %Vaccination Rollout by Sector
        NNnext              = NNbar;
        NNnext(nSectors+[1,2])    = 1;
        NNnext([1:nSectors,nSectors+3]) = NNnext([1:nSectors,nSectors+3])/sumWorkingAge;
        NNnext(end)         = 1;
        p2.NNnext = NNnext;

        isequence = [isequence; [t0 i]];
        
        [tout,Iclass,Iaclass,Isclass,Hclass,Dclass,p3,p4,betamod,y0,inext,still_susc]=...
         integr8(data,D,i,t0,tend,dis,y0,p2);
        if inext==0
            tend = tout(end);
        end
        
        Tout       = [Tout;tout(2:end)];  
        Iout       = [Iout;Iclass(2:end,:)];
        Iaout     = [Iaout;Iaclass(2:end,:,:)];
        Isout     = [Isout;Isclass(2:end,:,:)];
        Hout       = [Hout;Hclass(2:end,:)];
        Dout       = [Dout;Dclass(2:end,:)]; 
        W   = Wit'.*ones(length(tout),nSectors);
        Wout       = [Wout;W(1:end-1,:)];    
        hw  = data.hw(i,:).*ones(length(tout),nSectors);
        hwout      = [hwout;hw(1:end-1,:)];
        p3out    = [p3out;p3(2:end)];
        p4out    = [p4out;p4(2:end)];
        betamodout = [betamodout;betamod(2:end)];
        Sout       = [Sout;still_susc(2:end,:)];
        
        if Tout(end)<tend
            data.tvec = [data.tvec(1:end-1),Tout(end),tend];
            t0 = Tout(end);
            i                     = inext;
        end   

    end

    %% OUTPUTS:  

    Wout  = [Wout;Wout(end,:)];
    hwout = [hwout;hwout(end,:)];
    returnobject = struct;
    returnobject.Tout = Tout;
    returnobject.Stotal = sum(Sout,2);
    returnobject.workers = Wout;
    returnobject.homeworkers = hwout;
    returnobject.Itot = sum(Iout,2);
    returnobject.Iamat = Iaout;
    returnobject.Ismat = Isout;
    returnobject.Htot = sum(Hout,2);
    returnobject.hospmat = Hout;
    returnobject.deathtot = sum(Dout,2);
    returnobject.death1 = sum(Dout(:,nSectors+1),2);
    returnobject.death2 = sum(Dout(:,nSectors+2),2);
    returnobject.death3 = sum(Dout(:,[1:nSectors,nSectors+3]),2);
    returnobject.death4 = sum(Dout(:,nSectors+4),2);
    returnobject.deathmat = Dout;
    returnobject.betamod = betamodout;
    pout = struct;
    pout.p3 = p3out;
    pout.p4 = p4out;
    returnobject.selfisolation = pout;
    returnobject.isequence = isequence;
  
end

%%

function [tout,Iclass,Iaclass,Isclass,Hclass,Dclass,p3,p4,betamod,y0new,inext,still_susc]=...
          integr8(data,D,i,t0,tend,dis,y0,p2)
      
    %% copy over only required items
    
    rundata = struct;
    
%     rundata.nSectors = data.nSectors;
%     rundata.hw = data.hw; 
%     rundata.inext = data.inext; 
    rundata.NNs = data.NNs; 
    rundata.Dvec = data.Dvec;
    rundata.tvec = data.tvec;
    rundata.compindex = data.compindex; 
    rundata.t_import = data.t_import; 
    rundata.imand = data.imand; 
    rundata.rel_stringency = data.rel_stringency;

    rundata.sd_baseline = data.sd_baseline;
    rundata.sd_death_coef = data.sd_death_coef;
    rundata.sd_mandate_coef = data.sd_mandate_coef;
    
    NN0 = data.NNs; 
    %% CALL

    compindex = data.compindex;
    S_index = compindex.S_index;
    E_index = compindex.E_index;
    I_index = compindex.I_index;
    H_index = compindex.H_index;
    D_index = compindex.D_index;
    
    nStrata = size(NN0,1);
    fun  = @(t,y)ODEs(rundata,D,i,t,dis,y,p2);
    strategy = data.strategy;

    if strcmp(strategy,'Elimination')
        options = odeset('Events',@(t,y)elimination(t,y,data,nStrata,dis,i,p2));
    elseif strcmp(strategy,'Economic Closures') || strcmp(strategy,'School Closures')
        options = odeset('Events',@(t,y)reactive_closures(t,y,data,nStrata,dis,i,p2));
    elseif strcmp(strategy,'No Closures')
        options = odeset('Events',@(t,y)unmitigated(t,y,data,nStrata,dis,i,p2));
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
        if ie <= length(data.inext)
            inext = data.inext(ie(end));
        elseif ie == length(data.inext)+1 % importation event
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
    D     = yout(:,(D_index(1)-1)*nStrata + indices);
    still_susc = sum(yout(:,(repelem(S_index([1,3:end]),1,length(indices))-1)*nStrata + repmat(indices,1,length(S_index)-1)),2);

    Iclass   = Ia + Is + Iav1 + Isv1 + Iav2 + Isv2; 
    Hclass   = H + Hv1 + Hv2; 
    Dclass   = D;

    %% TIME-DEPENDENT PARAMETERS
    
    occ   = max(1,sum(Hclass,2));
    mu = zeros(size(Hclass));
    for ii = 1:size(mu,1)
        dis2 = update_hosp_dis_parameters(occ(ii), p2, dis);
        mu(ii,:)    = dis2.mu;
    end
    
    sumNN0 = sum(NN0);
    ddk = 10^6*sum(mu.*Hclass,2)/sumNN0;
    betamod = betamod_wrapped(ddk, data, i);
    
    [p3, p4] = fraction_averted_self_isolating(sum(Iclass,2), sumNN0, p2, tout, i);
    
    Iaclass = zeros(size(Ia,1),size(Ia,2),3);
    Isclass = zeros(size(Ia,1),size(Ia,2),3);
    Iaclass(:,:,1) = Ia; 
    Iaclass(:,:,2) = Iav1; 
    Iaclass(:,:,3) = Iav2; 
    Isclass(:,:,1) = Is; 
    Isclass(:,:,2) = Isv1; 
    Isclass(:,:,3) = Isv2; 
    
    %% compute Rt at end of period
    
%     S    = y_mat(:,compindex.S_index(1));
%     Sn    = y_mat(:,compindex.S_index(2));
%     S01   = y_mat(:,compindex.S_index(3));
%     Sv1   = y_mat(:,compindex.S_index(4));
%     Sv2   = y_mat(:,compindex.S_index(5));
%     
%     dis2 = update_vax_dis_parameters(dis2, S, Sn, compindex, y_mat);
%     
%     betam = social_distancing(p2.sdl,p2.sdb,ddk,data.rel_mobility(i));
%     
%     Rtret = get_R(nStrata,dis2,Stest,Stest1,Stest2,...
%             data.NNvec(:,3),data.Dvec(:,:,3),dis.beta,betam(end),pout(end),pout(end));
%         disp([2 sum(Iclass(end,:))/10^6 Ip(end) pout(end) Rtret])

end

%%

function [f] = ODEs(data,D,i,t,dis,y,p2)

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
        Ina,Ins,Inav1,Insv1,Inav2,Insv2,D);
    

    %% VACCINATION
    
    [v1rates, v1rater, v2rates, v2rater, v12rates, v12rater] = ...
        get_vax_rates(p2, t, nStrata, R,S,DE, Rv1,Sv1);

    if t<83
%         betamod = betamod_wrapped(10^6*sum(dis2.mu.*hospital_occupancy)/sum(NN0),p2, data, i);
%         Rt1 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,NN0,data.Dvec(:,:,1),dis2.beta,1,0,0);
%         disp([t foi'])
    end

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

%%



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
    
    hospital_occupancy = H + Hv1 + Hv2;
    sumH = sum(hospital_occupancy);
    occ   = max(1,sumH);     
    
    %% isolating
    [p3, p4] = fraction_averted_self_isolating(sum(sum(y_mat(:,I_index))), sumN, p2, t, i);
    
    %% correct ph for previous infection
    
    dis2 = update_vax_dis_parameters(dis, S, Sn, compindex, y_mat);
     % correct for hosp occupancy
    dis2 = update_hosp_dis_parameters(occ, p2, dis2);
    
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
    % t is greater than the end of the vaccine rollout: otherval = 0
    tval = min(t-(max(p2.tpoints)+7),0);
    % hval: no patients
    hval = min(p2.hosp_final_threshold - sumH,0);
    R3flag = ival + tval + hval;
    if ival==0 && tval==0 && hval==0
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
%         Rt3 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,dis.beta,p3,p4, ddk, data, 5);
        Rthresh = exp(dis.generation_time*log(2) / p2.final_doubling_time_threshold); % R value for which doubling time is 30 days
        R3flag = min(Rthresh - R_est,0);
%         disp([t/1000 Rt3 betamod p3 p4 sumH ])
    end
    value(7)      =  R3flag;
    direction(7)  = 1;
    isterminal(7) = 1;
    
end 

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
    
    value(4)      = abs(i-3) + abs(min(t-(data.tvec(end-1)+7),0)) + max(0,occ-p2.thl);
    direction(4)  = -1;
    isterminal(4) = 1;
    
    %% Event 5: End
    
    % not in hard lockdown: i = 1, 2 or 4; ivals = 0
    ivals = -abs((i-1)*(i-2)*(i-4));
    % have passed penultimate timepoint: tval = 0
    tval = min(t-(data.tvec(end-1)+0.1),0);
    % (low growth rate OR occupancy is low) AND have reached end of vaccine
    % rollout: otherval = 0
    otherval = -abs(min(0.025-r,0)*max(0,occ-p2.thl)) + min(t-7-max(p2.tpoints),0);
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
    tval = min(t-(max(p2.tpoints)+7),0);
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
    direction(7)  = 1;
    isterminal(7) = 1;
    
end

function [value,isterminal,direction] = unmitigated(t,y,data,nStrata,dis,i,p2)
    

    y_mat = reshape(y,nStrata,[]);
    compindex = data.compindex;
    S_index = compindex.S_index;
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
    
    occ  = max(1,sum(H+Hv1+Hv2)); 
        
    %% get waning for R
    
    dis2 = update_vax_dis_parameters(dis, S, Sn, compindex, y_mat);
    
    %% isolating
    sumN = sum(data.NNs);
    [p3, p4] = fraction_averted_self_isolating(sum(sum(y_mat(:,compindex.I_index))), sumN, p2, t, i);
          
    %% distancing
    ddk    = 10^6*sum(dis2.mu.*(H + Hv1 + Hv2))/sumN;
%     betamod = social_distancing(data.sd_baseline, data.sd_death_coef, data.sd_mandate_coef,ddk, data.rel_stringency(i));
     

    %% Event 1: Response Time
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0);
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: Late Lockdown
    
    value(2)      = - abs(i-1) + min(occ-0.95*p2.Hmax,0);
    direction(2)  = 1;
    isterminal(2) = 1;

    %% Event 3: End
    
    % not in i>=5; ivals = 0
    ivals = -floor(i/5);
    % have passed penultimate timepoint: tval = 0
    tval = min(t-(data.tvec(end-1)+0.1),0);
    % have reached end of vaccine rollout: otherval = 0
    otherval = min(t-max(p2.tpoints),0);
    R2flag = otherval + ivals + tval;
    if ivals==0 && tval==0 
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        if otherval~=0 && R_est<1
        % only compute R if R2flag is not already 0 and ivals and tval
        % conditions are both met
            Rt2 = get_R(nStrata,dis2,S+S01+S02,Sv1+S12,Sv2,dis.beta,p3,p4, ddk, data, 5);
            R2flag = min(1.00-Rt2,0);
        end
    end
    
    value(3)      = R2flag; 
    %measures can be removed at any stage if (Rt<1) or (after end of vaccination campaign)
    direction(3)  = 0;
    isterminal(3) = 1;
    
    %% Event 4: importation
    
    value(4)      =  min(t-data.t_import,0);
    direction(4)  = 1;
    isterminal(4) = 1;
    
    %% Event 5: end
    % i is 5
    ival = -abs((i-5));
    % t is greater than the end of the vaccine rollout: otherval = 0
    tval = min(t-(max(p2.tpoints)+7),0);
    % hval: no patients
    sumH = sum(H + Hv1 + Hv2);
    hval = min(p2.hosp_final_threshold-sumH,0);
    R3flag = ival + tval + hval;
    if ival==0 && tval==0 && hval==0
%         Rt3 = get_R(nStrata,dis2,S+S01,Sv1,Sv2,dis.beta,p3,p4, ddk, data, 5);
        Rthresh = exp(dis.generation_time*log(2) / p2.final_doubling_time_threshold); % R value for which doubling time is 30 days
        R_est = get_R_est(dis2, compindex, y_mat, p3, p4); 
        R3flag = min(Rthresh - R_est,0);
%         disp([t/1000 Rt3 sumH])
    end
    value(5)      =  R3flag;
    direction(5)  = 1;
    isterminal(5) = 1;
    
    
end

