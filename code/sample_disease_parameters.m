function param_struct = sample_disease_parameters(nsamples)

sevenpathogens = readtable('../data/sevenpathogens.csv');

%% single parameters
ircolumns = cell2mat(cellfun(@(a) strmatch(a, sevenpathogens.Properties.VariableNames),...
    {'ihr','ifr'},'uniform',false));

disparams = sevenpathogens(:,setdiff(1:size(sevenpathogens,2),ircolumns));

disparams.Ti = repmat(365000000,size(disparams,1),1);

sample_struct = struct;
% replicate non covid
for i = 1:size(disparams,2)
    thisparam = disparams.Properties.VariableNames{i};
    samples = disparams.(thisparam);
    sample_struct.(thisparam) = [samples; samples(1:4); samples(1:4)];
end
    
param_struct = struct;
hyper_param_struct = struct;

for i = 1:size(disparams,2)
    thisparam = disparams.Properties.VariableNames{i};
    samples = sample_struct.(thisparam);
    pHat = gamfit(samples);
    if strmatch('R0',thisparam)
%         samples = samples(samples<5.5) - 1.5;
        [p1 p2] = normfit(samples);
        pHat = [p1 p2];
    end
    hyper_param_struct.(thisparam) = pHat;
end

set1 = {'Tsr','Tsh','R0'};
set2 = {'Threc','Thd'};
sets = {set1; set2};

rho = eye(length(disparams.Properties.VariableNames));
% disparams.R0(disparams.R0>4) = 4;
for seti = 1:length(sets)
    pis = cell2mat(cellfun(@(a) strmatch(a, disparams.Properties.VariableNames),...
        sets{seti},'uniform',false));
    orig_samples = table2array(disparams(:,pis));
    correlation = corr(orig_samples);
    rho(pis,pis) = correlation;
    disp(sets{seti})
    disp(correlation)
end

mu  = zeros(1,length(disparams.Properties.VariableNames));   
Z = mvnrnd(mu, rho, nsamples); %Generate multivariate corralated random number
U = normcdf(Z,0,1);     %Compute the CDF

for i = 1:size(disparams,2)
    thisparam = disparams.Properties.VariableNames{i};
    pHat = hyper_param_struct.(thisparam);
    newparams = gaminv(U(:,i),pHat(1),pHat(2));
    if strmatch('R0',thisparam)
%         newparams = min(newparams + 1.5,4);
        pd = makedist('Normal',pHat(1),pHat(2));
        t = truncate(pd,1.5,4);
        newparams = icdf(t,U(:,i));
    end
    newparams(isnan(newparams)) = mean(sample_struct.(thisparam));
    param_struct.(thisparam) = newparams;
end

% figure('Position', [100 100 400 300]); hist(param_struct.R0)

% use beta for ps
thisparam = 'ps';
samples = disparams.(thisparam);
pHat = betafit(samples);
hyper_param_struct.(thisparam) = pHat;
newparams = betarnd(pHat(1),pHat(2),nsamples,1);
param_struct.(thisparam) = newparams;
disp(hyper_param_struct)

%% ihr and hfr 

ihrs = table2array(sevenpathogens(:,ircolumns(:,1)));
ifrs = table2array(sevenpathogens(:,ircolumns(:,2)));
hfrs = ifrs./ihrs;

maxihrs = mean(ihrs');
maxhfrs = mean(hfrs');

pHat = betafit(maxihrs);
maxihr = betarnd(pHat(1),pHat(2),nsamples,1);
%%!! ihr cannot exceed probability symptomatic
while(sum(maxihr>param_struct.ps)>0)
    resample = find(maxihr>param_struct.ps);
    newihr = betarnd(pHat(1),pHat(2),length(resample),1);
    maxihr(resample) = newihr;
end

pHat = betafit(maxhfrs);
%%!! max val = 1
maxhfr = betarnd(pHat(1),pHat(2),nsamples,1);
% maxhfr = unifrnd(0,1,nsamples,1);


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

%% CEPI ifr

newifr = readtable('../data/IFR Distributions.xlsx','ReadRowNames',true).IFR(1:17);
% ps > ihr > ifr
% hfr > ifr
upperbound = sevenpathogens.ps;
vals = mean(ihrs')';
lowerbound = mean(ifrs')';

transformedvals = (vals-lowerbound)./(upperbound-lowerbound);

pHat = betafit(transformedvals);
relvals = betarnd(pHat(1),pHat(2),nsamples,1);


p1=7; newparams = gamrnd(p1,1/p1,nsamples,1);

for i = 1:nsamples
    param_struct.ifr(i,:) = newifr'*newparams(i);
    param_struct.ihr(i,:) = relvals(i) * (param_struct.ps(i) - param_struct.ifr(i,:)) + param_struct.ifr(i,:); 
end


end