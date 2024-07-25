% rerun model starting from end of mitigation, trying one of two exit
% strategies
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% p2: struct of p2 intervention parameters
% returned: struct of outcomes from first run, which had no exit strategy 
%
% returnobject: struct of epi outcomes 


function returnobject = reRun(data, dis, p2, returned)    
    
    
    %% PARAMETERS
    isequence = returned.isequence;
    restart_time = isequence(find(isequence(:,2)>4,1),1);
    i = isequence(find(isequence(:,2)>4,1)-1,2);
    isequence = isequence(find(isequence(:,2)<5),:);
    
    nSectors = data.nSectors;

    % initial conditions
    t0 = restart_time;
    last_index = find(returned.Tout >= t0, 1);
    Tout = returned.Tout(1:last_index);
    Sout = returned.Sout(1:last_index,:);
    Wout = returned.workers(1:(last_index-1),:);
    hwout = returned.homeworkers(1:(last_index-1),:);
    Iout = returned.Iout(1:last_index,:);
    Iaout = returned.Iamat(1:last_index,:,:);
    Isout = returned.Ismat(1:last_index,:,:);
    Hout = returned.hospmat(1:last_index,:);
    Dout = returned.deathmat(1:last_index,:);
    betamodout = returned.betamod(1:last_index);
    pout = returned.selfisolation;
    p3out = pout.p3(1:last_index);
    p4out = pout.p4(1:last_index);
    y0 = returned.y0;
    
    %% LOOP

    i = 6;
    if data.exittype==2
        data.Dvec(:,:,6) = data.Dvec(:,:,1);
        data.workerConfigMat(:,6) = data.workerConfigMat(:,1);
    end

    tend = data.tvec(end);

    while Tout(end)<tend 

        Wit              = data.workerConfigMat(:,i);    
        contact_matrix   = data.Dvec(:,:,i);
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
            i  = inext;
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
%     returnobject.y0 = y0;
  
end


