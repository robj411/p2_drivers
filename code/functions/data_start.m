% load in basic data and make general definitions
%
% data: struct of general model parameters

function data = data_start(lbfile)

    %% closure strategies
    data      = struct; %rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    closurefile = '../data/closures.xlsx';
    sheets = sheetnames(closurefile);
    for i = 1:numel(sheets)
        thissheet = sheets{i};
        closurei = table2array(readtable(closurefile,'FileType','spreadsheet','Sheet',sheets{i}));
        data.(thissheet) = closurei;
    end
    
    %% general variables
    data.adInd = 3;
    data.nSectors = length(data.x_elim);
    data.tvec = [0 365*10];
    data.EdInd    = 41;%education sector index
    data.HospInd  = [32,43,44];%hospitality sector indices
    data.ageindex = {1,2:4,5:13,14:21};
    
    %% contacts
    contacts = struct;
    contacts.sectorcontacts = readtable('../data/sectorcontacts.csv');
    contacts.sectorcontactfracs = readtable('../data/sec_contact_dist_UK.csv');
    data.contacts = contacts;
    
    %% indices
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
    
    %% vaccine rollout
    fulltable = readtable(lbfile,'FileType','spreadsheet','Sheet','Vx timeline');
    % fulltable = readtable('../data/20250328 pandemic delivery scenarios - for ICL.xlsx','FileType','spreadsheet','Sheet',1);
    % fulltable = readtable('../data/20240611 LB Daily Vaccine Delivery.xlsx','FileType','spreadsheet','Sheet',1);
    colnames = regexprep(fulltable.Properties.VariableNames(2+[1:3]),'s','');
    vaxtab = table2array(fulltable);
    nScen = (length(fulltable.Properties.VariableNames)-2)/6;
    scenarios = cell(1,nScen);
    for i = 1:nScen
        scenbpsv = array2table(vaxtab(:, 2+(i-1)*6 + [1:3]));
        scenspec = array2table(vaxtab(:, 2+(i-1)*6 + [4:6]));
        scenbpsv.Properties.VariableNames = colnames;
        scenspec.Properties.VariableNames = colnames;
        scenarios{i} = {scenbpsv, scenspec};        
    end
    
    data.scenarios = scenarios;

end



