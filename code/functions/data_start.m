% load in basic data and make general definitions
%
% data: struct of general model parameters

function data = data_start()

    data      = struct; %rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    closurefile = '../data/closures.xlsx';
    sheets = sheetnames(closurefile);
    for i = 1:numel(sheets)
        thissheet = sheets{i};
        closurei = table2array(readtable(closurefile,'FileType','spreadsheet','Sheet',sheets{i}));
        data.(thissheet) = closurei;
    end
    
    data.adInd = 3;
    data.nSectors = length(data.x_elim);
    data.tvec = [0 365*10];
    data.EdInd    = 41;%education sector index
    data.HospInd  = [32,43,44];%hospitality sector indices
    
    contacts = struct;
    contacts.sectorcontacts = readtable('../data/sectorcontacts.csv');
    contacts.sectorcontactfracs = readtable('../data/sec_contact_dist_UK.csv');
    data.contacts = contacts;
    
    data.ageindex = {1,2:4,5:13,14:21};
    
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



