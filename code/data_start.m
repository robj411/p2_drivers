function data = data_start()

    load('../country_mats/Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields    = fieldnames(data);
    ikeep     = [6,7,8,13,14,16,17,18];
    data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    data.adInd = 3;
    data.lx = length(data.B);
    data.tvec = [-75 365*2];
    data.alp = 1;
    data.EdInd    = 41;%education sector index
    data.HospInd  = [32,43,44];%hospitality sector indices
    
    compindex = struct;
    compindex.S_index = [1,10,11,12,13]; % S, Sn, Sp, Sv, Sb
    compindex.E_index = [2,14,15];
    compindex.I_index = [3:6,16:19,20:23];
%     Ina=    y_mat(:,compindex.I_index(1));
%     Isa=    y_mat(:,compindex.I_index(2));
%     Ins=    y_mat(:,compindex.I_index(3));
%     Iss=    y_mat(:,compindex.I_index(4));
    compindex.H_index = [7,24:25];
    compindex.R_index = [8,26:27];
    compindex.D_index = [9];
    compindex.V_index = [28];
    data.compindex = compindex;


end
