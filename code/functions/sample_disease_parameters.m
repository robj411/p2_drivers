% samples many pathogen values from distributions

% nsamples: the number of samples to take

% param_struct: struct containing all values sampled

function param_struct = sample_disease_parameters(nsamples)

sevenpathogens = readtable('../data/sevenpathogens.csv');

%% single parameters
ircolumns = cell2mat(cellfun(@(a) strmatch(a, sevenpathogens.Properties.VariableNames),...
    {'ihr','ifr'},'uniform',false));

disparams = sevenpathogens(:,setdiff(1:size(sevenpathogens,2),ircolumns));

% no waning
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
        t = truncate(pd,2,4.5);
        newparams = icdf(t,U(:,i));
    end
    newparams(isnan(newparams)) = mean(sample_struct.(thisparam));
    param_struct.(thisparam) = newparams;
end

% figure('Position', [100 100 400 300]); hist(param_struct.R0)

%% use beta for ps and red, with parameters 5,5
% thisparam = 'ps';
% samples = disparams.(thisparam);
% pHat = betafit(samples);
% hyper_param_struct.(thisparam) = pHat;
% newparams = betarnd(pHat(1),pHat(2),nsamples,1);
newparams = betarnd(5,5,nsamples,1);
param_struct.ps = newparams;
newparams = betarnd(5,5,nsamples,1);
param_struct.red = newparams;

param_struct.frac_presymptomatic = unifrnd(0,1,nsamples,1);

disp(hyper_param_struct)

%% ihr and hfr 

ihrs = table2array(sevenpathogens(:,ircolumns(1:17,1)));
ifrs = table2array(sevenpathogens(:,ircolumns(1:17,2)));
hfrs = ifrs./ihrs;

maxihrs = max(ihrs');
maxifrs = max(ifrs');
%%!! is the last val
maxhfrs = hfrs(:,end); %max(hfrs(hfrs<1)');
% repeat non-covid outcomes
% maxihrs = [maxihrs maxihrs];
% maxhfrs = [maxhfrs; maxhfrs];

offdiag = corr(maxihrs',maxhfrs);

Z = mvnrnd([0 0], [1 offdiag; offdiag 1], nsamples); %Generate multivariate corralated random number
U = normcdf(Z,0,1);     %Compute the CDF

pHat = betafit(maxihrs);
pHat = pHat./min(.5,min(pHat));
maxihr = betainv(U(:,1),pHat(1),pHat(2));
%%!! ihr cannot exceed probability symptomatic
while(sum(maxihr>param_struct.ps)>0)
    resample = find(maxihr>param_struct.ps);
    newihr = betarnd(pHat(1),pHat(2),length(resample),1);
    maxihr(resample) = newihr;
end

pHat = betafit(maxhfrs);
pHat = pHat./min(.5,min(pHat));
%%!! max val = 1
maxhfr = betainv(U(:,2),pHat(1),pHat(2));
% maxhfr = unifrnd(0,1,nsamples,1);
% figure; hist(maxhfr.*maxihr)
% mean(maxhfr.*maxihr)

ihrrr = exp(table2array(readtable('../data/ihrrr.csv')));
hfrrr = exp(table2array(readtable('../data/hfrrr.csv')));

maxihrrr = max(ihrrr');
maxhfrrr = max(hfrrr');

param_struct.ihr = zeros(nsamples, size(ihrrr,2));
param_struct.ifr = zeros(nsamples, size(ihrrr,2));
hfr = zeros(nsamples, size(ihrrr,2));

for i = 1:nsamples
    param_struct.ihr(i,:) = ihrrr(i,:) ./ maxihrrr(i) .* maxihr(i);
    hfr(i,:) = hfrrr(i,:) ./ maxhfrrr(i) .* maxhfr(i);
    param_struct.ifr(i,:) = param_struct.ihr(i,:) .* hfr(i,:);
end

plotfun = @(x) log(x'./repmat(x(:,1),1,size(x,2))');
plotfun = @(x) log(x');
h = figure('Position', [100 100 900 300]); 
subplot(1,3,1); plot(plotfun(param_struct.ihr),'color',[.5 .5 .5 .05]); hold on
subplot(1,3,1); plot(plotfun(ihrs)); title('IHR')
subplot(1,3,2); plot(plotfun(hfr),'color',[.5 .5 .5 .05]); hold on
subplot(1,3,2); plot(plotfun(hfrs)); title('HFR')
subplot(1,3,3); plot(plotfun(param_struct.ifr),'color',[.5 .5 .5 .05]); hold on
subplot(1,3,3); plot(plotfun(ifrs)); title('IFR')
saveas(h,'../figures/ratesbyage','jpg');
close gcf

%% CEPI ifr

% newifr = readtable('../data/IFR Distributions.xlsx','ReadRowNames',true).IFR(1:17);
% % ps > ihr > ifr
% % hfr > ifr
% upperbound = sevenpathogens.ps;
% vals = mean(ihrs')';
% lowerbound = mean(ifrs')';
% 
% transformedvals = (vals-lowerbound)./(upperbound-lowerbound);
% 
% pHat = betafit(transformedvals);
% relvals = betarnd(pHat(1),pHat(2),nsamples,1);
% 
% 
% p1=7; newparams = gamrnd(p1,1/p1,nsamples,1);
% 
% for i = 1:nsamples
%     param_struct.ifr(i,:) = newifr'*newparams(i);
%     param_struct.ihr(i,:) = relvals(i) * (param_struct.ps(i) - param_struct.ifr(i,:)) + param_struct.ifr(i,:); 
% end


end