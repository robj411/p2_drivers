% load in basic data and make general definitions
%
% data: struct of general model parameters

function data = data_start()

    load('../country_mats/Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields    = fieldnames(data);
    ikeep     = [6,7,8,13,14,16,17,18];
    data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    data.adInd = 3;
    data.nSectors = length(data.B);
    data.tvec = [0 365*10];
    data.EdInd    = 41;%education sector index
    data.HospInd  = [32,43,44];%hospitality sector indices
    
    contacts = struct;
    contacts.B = data.B;
    contacts.C = data.C;
    data = rmfield(data,'B');
    data = rmfield(data,'C');
    data = rmfield(data,'hospA2');
    data = rmfield(data,'hospA3');
    data = rmfield(data,'hospA4');
    data.contacts = contacts;
    
    compindex = struct;
    
    compindex.S_index = [1, 8:11,24:25]; % S, Sn, S01, Sv, Sb, S02, S12
    compindex.E_index = [2, 12:13];
    compindex.I_index = [3:4, 14:15,16:17];
%     Ia=    y_mat(:,compindex.I_index(1));
%     Is=    y_mat(:,compindex.I_index(2));
    compindex.H_index = [5, 18:19];
    compindex.R_index = [6, 20:21];
    compindex.D_index = [7];
    compindex.V_index = [22:23];
    
    compindex.vaccine = [compindex.S_index(4),compindex.S_index(7),...
        compindex.E_index(2),...
        compindex.I_index(3),...
        compindex.I_index(4),...
        compindex.H_index(2),...
        compindex.R_index(2)];
    
    compindex.booster = [compindex.S_index(5),...
        compindex.E_index(3),...
        compindex.I_index(5),...
        compindex.I_index(6),...
        compindex.H_index(3),...
        compindex.R_index(3)];    

    data.compindex = compindex;


end
