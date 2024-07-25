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
    
    data = get_strategy_design(data, strategy, p2);
    
    %% get contact matrices
    
    nSectors = data.nSectors;
    nStrata = length(data.NNs);
    nConfigs = length(data.configuration)/data.nSectors;

    NNbar = data.NNs;
    configMat = reshape(data.configuration,nSectors,nConfigs);
    data.workerConfigMat = configMat;
    
    Dvec = zeros(nStrata,nStrata,nConfigs);
    for i = 1:nConfigs
        % get matrix per configuration
        Dtemp = p2MakeDs(data,NNbar,configMat(:,i),data.hw(i,:));
        % store matrix in list
        Dvec(:,:,i) = Dtemp;
    end
    
    data.Dvec = Dvec;
    data.strategy = strategy;

    data.basic_foi = Dvec(:,:,1)*(1./NNbar);
    
    [data, returnobject] = p2SimVax(data, dis, p2);

end

%%

function [data, returnobject] = p2SimVax(data, dis, p2)    
    
    
    %% PARAMETERS
    nStrata = size(data.NNs,1);
    nSectors = data.nSectors;
    NNbar = data.NNs;
    sumWorkingAge = sum(NNbar([1:nSectors,nSectors+3]));
    compindex = data.compindex;
    nODEs = max(struct2array(compindex));
    S0 = NNbar;

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

        Wit               = data.workerConfigMat(:,i);    
        NNfeed            = NNbar;
        NNfeed(NNfeed==0) = 1;
        contact_matrix                 = data.Dvec(:,:,i);

        %Vaccination Rollout by Sector
        NNnext              = NNbar;
        NNnext(nSectors+[1,2])    = 1;
        NNnext([1:nSectors,nSectors+3]) = NNnext([1:nSectors,nSectors+3])/sumWorkingAge;
        NNnext(end)         = 1;
        p2.NNnext = NNnext;

        isequence = [isequence; [t0 i]];
        
        [tout,Iclass,Iaclass,Isclass,Hclass,Dclass,p3,p4,betamod,y0,inext,still_susc,data]=...
         integr8(data,contact_matrix,i,t0,tend,dis,y0,p2);
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
            if i==5
                ysave = y0;
            end
        end   

    end

    %% OUTPUTS:  

    Wout  = [Wout;Wout(end,:)];
    hwout = [hwout;hwout(end,:)];
    returnobject = struct;
    returnobject.Tout = Tout;
    returnobject.Sout = Sout;
    returnobject.Stotal = sum(Sout,2);
    returnobject.workers = Wout;
    returnobject.homeworkers = hwout;
    returnobject.Iout = Iout;
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
%     returnobject.y0 = ysave;
  
end



