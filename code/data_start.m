function data = data_start()

    load('../country_mats/Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields    = fieldnames(data);
    ikeep     = [6,7,8,13,14,16,17,18];
    data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    data.adInd = 3;
    data.lx = length(data.B);
    data.tvec = [-75 365];
    data.alp = 1;
    data.EdInd    = 41;%education sector index
    data.HospInd  = [32,43,44];%hospitality sector indices

end
