function param_struct = sample_disease_parameters(nsamples)

sevenpathogens = readtable('../data/sevenpathogens.csv');

ircolumns = cell2mat(cellfun(@(a) strmatch(a, sevenpathogens.Properties.VariableNames),...
    {'ihr','ifr'},'uniform',false));

disparams = sevenpathogens(:,setdiff(1:size(sevenpathogens,2),ircolumns));

param_struct = struct;

for i = 1:size(disparams,2)
    thisparam = disparams.Properties.VariableNames{i};
    samples = disparams.(thisparam);
    if strmatch('R0',thisparam)
        samples = samples - 1;
    end
    pHat = lognfit(samples);
    newparams = lognrnd(pHat(1),pHat(2),nsamples,1);
    if strmatch('R0',thisparam)
        newparams = newparams + 1;
    end
    param_struct.(thisparam) = newparams;
end


ihrs = table2array(sevenpathogens(:,ircolumns(:,1)));
ifrs = table2array(sevenpathogens(:,ircolumns(:,2)));
hfrs = ifrs./ihrs;

maxihrs = mean(ihrs');
maxhfrs = mean(hfrs');

pHat = betafit(maxihrs);
maxihr = betarnd(pHat(1),pHat(2),nsamples,1);

pHat = betafit(maxhfrs);
%%!! max val = 1
maxhfr = betarnd(pHat(1),pHat(2),nsamples,1);
maxhfr = unifrnd(0,1,nsamples,1);


ihrrr = exp(table2array(readtable('ihrrr.csv')));
hfrrr = exp(table2array(readtable('hfrrr.csv')));

maxihrrr = max(ihrrr');
maxhfrrr = max(hfrrr');

param_struct.ihr = zeros(nsamples, size(ihrrr,2));
param_struct.ifr = zeros(nsamples, size(ihrrr,2));

for i = 1:nsamples
    param_struct.ihr(i,:) = ihrrr(i,:) ./ maxihrrr(i) .* maxihr(i);
    hfr = hfrrr(i,:) ./ maxhfrrr(i) .* maxhfr(i);
    param_struct.ifr(i,:) = param_struct.ihr(i,:) .* hfr;
end

end