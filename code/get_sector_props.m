
    addpath('../');
    
    load('Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields    = fieldnames(data);
    ikeep     = [6,7,8,13,14,16,17,18];
    data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    lx        = length(data.B);
    data.tvec = [-75 365*3+1];
    CD        = readtable('country_data.csv');
    nsamples  = 1000;
    inputs    = zeros(nsamples,493);
    sectorprops    = zeros(nsamples,lx);
    popprops    = zeros(nsamples,lx+1);
    R0samples    = zeros(nsamples,1);
    outputs   = zeros(nsamples,4);
    
    dis = get_dis_params('Covid Wildtype');  
    
    
    rng default;
    
    for i = 1:nsamples;
    
        ldata     = data;
        ldata     = p2RandCountry(ldata,CD);

        [ldata,~,~]    = p2Params(ldata,'Covid Wildtype',dis);%to define wnorm and Td_CWT
        [ldata,dis2,p2] = p2Params(ldata,'Covid Wildtype',dis);  
        popprops(i,:) = ldata.NNs([1:45,48])' ./ sum(ldata.NNs([1:45,48]));
        sectorprops(i,:) = ldata.obj' ./ sum(ldata.obj);
        R0samples(i) = dis2.R0sample;

    end
scatter(R0samples,sectorprops(:,1))    
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
nexttile   
hist(sectorprops(:,1))    
nexttile   
hist(sectorprops(:,31))  
nexttile   
hist(popprops(:,1))    
nexttile   
hist(popprops(:,31))   
    
       