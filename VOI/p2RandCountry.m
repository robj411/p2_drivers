function data = p2RandCountry(data,CD,income_level)



remote_quantile = unifrnd(0,1);
labsh_quantile = unifrnd(0,1);

data.remote_quantile = remote_quantile;

%%!! hard coded
if strcmp(income_level,'LLMIC')
    country_indices = strcmp(CD.igroup,'LIC') | strcmp(CD.igroup,'LMIC');
    internet_coverage = betainv(remote_quantile,1.78,3.11);
    labsh = betainv(labsh_quantile,5.09,4.52);
    Hmax = gamrnd(1.296265366, 1/0.049499412 );
elseif strcmp(income_level,'MIC')
    country_indices = strcmp(CD.igroup,'UMIC') | strcmp(CD.igroup,'LMIC');
    internet_coverage = betainv(remote_quantile,3.77,2.91);    
    labsh = betainv(labsh_quantile,6.29,6.54);
    Hmax = gamrnd(1.360805416, 1/0.027592345);
elseif strcmp(income_level,'HIC')
    country_indices = strcmp(CD.igroup,'HIC');
    internet_coverage = betainv(remote_quantile,9.57,1.39);
    labsh = betainv(labsh_quantile,7.97,6.87);
    Hmax = gamrnd(2.045582620,   1/0.021471378);
else
    country_indices = strcmp(CD.igroup,'LIC') | strcmp(CD.igroup,'LMIC') |...
        strcmp(CD.igroup,'UMIC') | strcmp(CD.igroup,'HIC');
    internet_coverage = 0;
    labsh = 0;
    Hmax = lognrnd(2.83235908, 0.87983858);
end



%Npop
nonempind = find(~isnan(CD.Npop1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,4:24));
defivalue = 50*10^6*randvalue'/sum(randvalue);
data.Npop = defivalue;

%NNs
nonempind = find(~isnan(CD.NNs1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,25:69));%number of workers by sector in real country
defivalue = randvalue/sum(table2array(CD(randindex,3+[5:13])));%proportion of adult population by sector in real country
defivalue = sum(data.Npop(5:13))*defivalue;%number of workers by sector in artificial country
defivalue = [defivalue,data.Npop(1),sum(data.Npop(2:4)),sum(data.Npop(5:13))-sum(defivalue),sum(data.Npop(14:end))]';
data.NNs  = defivalue;

%CM
nonempind = find(~isnan(CD.CMaa) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,70:325));
defivalue = reshape(randvalue,16,16);
data.CM   = defivalue;

%comm
nonempind = find(~isnan(CD.comm) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,326));
defivalue = randvalue;
data.comm = defivalue;

%travelA3
nonempind     = find(~isnan(CD.travelA3) & country_indices);
randindex     = nonempind(randi(numel(nonempind)));
randvalue     = table2array(CD(randindex,327));
defivalue     = randvalue;
data.travelA3 = defivalue;

%schoolA1&schoolA2
nonempind     = find(~isnan(CD.schoolA1) & country_indices);
randindex     = nonempind(randi(numel(nonempind)));
randvalue1    = table2array(CD(randindex,328));
randvalue2    = table2array(CD(randindex,329));
data.schoolA1 = randvalue1;
data.schoolA2 = randvalue2;

%workp
nonempind  = find(~isnan(CD.workp) & country_indices);
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,330));
defivalue  = randvalue;
data.workp = defivalue;

% change sizes of tourism and agricultural sectors
% rescale sectors
% get populations
agInd = 1;
FAAind = 32;
agPop = data.NNs(agInd);
FAAPop = data.NNs(FAAind);
otherindices = setdiff([1:45,48],[agInd,FAAind]);

workingagepop = sum(data.NNs([1:45,48]));
leftoverpop = workingagepop - agPop - FAAPop;

popindices = find(~isnan(CD.NNs1)&~isnan(CD.Npop1) & country_indices);
workers = table2array(CD(popindices,25:69));%number of workers by sector in real country
adults = sum(table2array(CD(popindices,8:16)),2);

minFAAworkers = min(workers(:,FAAind) ./ adults);
minagWorkers = min(workers(:,agInd) ./ adults);
maxFAAworkers = max(workers(:,FAAind) ./ sum(workers,2));
maxagWorkers = max(workers(:,agInd) ./ sum(workers,2));

FAAfrac = unifrnd( minFAAworkers , maxFAAworkers );
agFrac = unifrnd( minagWorkers, min(maxagWorkers, 1-FAAfrac) );
leftover_frac = 1 - FAAfrac - agFrac;

data.NNs(agInd) = agFrac * workingagepop;
data.NNs(FAAind) = FAAfrac * workingagepop;
data.NNs(otherindices) = data.NNs(otherindices)./leftoverpop .* workingagepop .* leftover_frac;

%obj: income per worker
nonempind                   = find(~isnan(CD.obj1)&~isnan(CD.NNs1) & country_indices);
randindex                   = nonempind(randi(numel(nonempind)));
randvalue                   = table2array(CD(randindex,331:375));%gva by sector in real country
defivalue                   = randvalue./table2array(CD(randindex,25:69));%gva per worker by sector in real country
defivalue(isnan(defivalue)) = 0;
defivalue(isinf(defivalue)) = 0;
defivalue                   = data.NNs(1:45).*defivalue';%gva by sector in artificial country
%%!! need to check food and accomm gva as % of gdp
data.obj                    = defivalue;


pointiness = 5.93; % for beta distribution
sec_to_international = 3.66; % scales fraction to fraction
international_const = 0.0959; % constant

remaining_international_tourism = lognrnd(-1.39,0.39);
GDP = sum(data.obj);
FAAfrac = data.obj(FAAind)/GDP;
alpha = pointiness * min(sec_to_international * FAAfrac + international_const,1);
beta = pointiness - alpha;
frac_tourism_international = betarnd(alpha,beta);
FAAmax = 1 - frac_tourism_international + frac_tourism_international*remaining_international_tourism;
data.frac_tourism_international = frac_tourism_international;

data.unmit = ones(size(data.x_elim));
data.unmit(FAAind) = FAAmax;
data.x_elim(FAAind) = min(FAAmax,data.x_elim(FAAind));
data.x_econ(FAAind,:) = min(FAAmax,data.x_econ(FAAind,:));
data.x_schc(FAAind,:) = min(FAAmax,data.x_schc(FAAind,:));

%wfh
nonempind = find(~isnan(CD.wfhl1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,376:465));
defivalue = reshape(randvalue,45,2)';

mins = min(table2array(CD(nonempind,376:420)));
maxs = max(table2array(CD(nonempind,421:465)));
% props = reshape([mins,maxs],45,2)';
% props(props==0) = 1e-6;
% sigmoid = -log(1./props - 1);
% meansig = mean(sigmoid);
% sdsig = (meansig-sigmoid(1,:))./norminv(.95);
% qsig = norminv(quantile,meansig,sdsig);
% newprop = 1./(1+exp(-qsig));
newprop = unifinv(remote_quantile,mins,maxs);

data.wfh  = defivalue;
data.wfh  = [newprop; newprop];

%Tres
nonempind = find(~isnan(CD.Tres) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,472));
defivalue = randvalue;
data.Tres = defivalue;

%t_tit
nonempind  = find(~isnan(CD.t_tit) & country_indices);
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,470));
defivalue  = randvalue;
data.t_tit = defivalue;

%trate
nonempind  = find(~isnan(CD.trate) & country_indices);
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,471));
defivalue  = randvalue;
data.trate = defivalue;

%sdl
nonempind = find(~isnan(CD.sdl) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,473));
defivalue = randvalue;
data.sdl  = defivalue;

%sdb
nonempind = find(~isnan(CD.sdb) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,474));
defivalue = randvalue;
data.sdb  = defivalue;

%Hmax
nonempind = find(~isnan(CD.Hmax) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,469));
defivalue = randvalue;
data.Hmax = Hmax; %defivalue;

%t_vax
nonempind  = find(~isnan(CD.t_vax) & country_indices);
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,466));
defivalue  = randvalue;
data.t_vax = defivalue;

%arate
nonempind  = find(~isnan(CD.arate) & country_indices);
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,467));
defivalue  = randvalue;
data.arate = defivalue;

%puptake
nonempind    = find(~isnan(CD.puptake) & country_indices);
randindex    = nonempind(randi(numel(nonempind)));
randvalue    = table2array(CD(randindex,468));
defivalue    = randvalue;
data.puptake = defivalue;

%la
nonempind = find(~isnan(CD.la1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,476:493));
defivalue = randvalue;
data.la   = defivalue;

%vly
na        = [data.Npop(1:17)',sum(data.Npop(18:end))];%length is 18 to match life table
gdp       = 365*sum(data.obj);
la        = data.la;
lg        = [dot(la(1),na(1))/sum(na(1)),...
             dot(la(2:4),na(2:4))/sum(na(2:4)),...
             dot(la(5:13),na(5:13))/sum(na(5:13)),...
             dot(la(14:end),na(14:end))/sum(na(14:end))];
for k = 1:length(lg); 
    lgh(k) = sum(1./((1+0.03).^[1:lg(k)]));
end 
vsl       = 160*gdp/sum(na);
defivalue = vsl/(dot(lgh,[na(1);sum(na(2:4));sum(na(5:13));sum(na(14:end))])/sum(na));
data.vly  = defivalue;
data.gdp = gdp;
data.gdppc = gdp/sum(data.Npop);

%vsy
working_years = 45;
discount_rate = 0.03;
rate_of_return_one_year = 0.08;
PV = (1-(1+discount_rate)^(-working_years))/discount_rate;

students_affected = 1 - 0.86 * internet_coverage;
mean_annual_income = labsh*gdp/workingagepop;

educationloss_all_students = ...
    PV*mean_annual_income*rate_of_return_one_year*students_affected*sum(data.Npop(2:4));
educationloss_per_student = ...
    PV*mean_annual_income*rate_of_return_one_year*students_affected;



defivalue = 0.5454*gdp/sum(data.Npop(2:4));
data.vsy  = educationloss_per_student;
data.educationloss_all_students  = educationloss_all_students;

end