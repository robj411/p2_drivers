function data = p2RandCountry(data,CD,income_level,dis_ref,R0betafun)

lx = data.lx;

remote_quantile = unifrnd(0,1);
labsh_quantile = unifrnd(0,1);
bmi_quantile = unifrnd(0,1);
international_tourism_quant = unifrnd(0,1);
remaining_tourism_quant = unifrnd(0,1);
Hmax_quantile = unifrnd(0,1);

data.remote_quantile = remote_quantile;

%%!! hard coded
if strcmp(income_level,'LLMIC')
    country_indices = strcmp(CD.igroup,'LIC') | strcmp(CD.igroup,'LMIC');
    internet_coverage = betainv(remote_quantile,1.78,3.11);
    labsh = betainv(labsh_quantile,5.09,4.52);
    bmi = norminv(bmi_quantile,24.29,2.28);
    Hmax = gaminv(Hmax_quantile, 1.296265366, 1/0.049499412 );
elseif strcmp(income_level,'MIC')
    country_indices = strcmp(CD.igroup,'UMIC') | strcmp(CD.igroup,'LMIC');
    internet_coverage = betainv(remote_quantile,3.77,2.91);    
    labsh = betainv(labsh_quantile,6.29,6.54);
    bmi = norminv(bmi_quantile,25.96,2.22);
    Hmax = gaminv(Hmax_quantile, 1.360805416, 1/0.027592345);
elseif strcmp(income_level,'HIC')
    country_indices = strcmp(CD.igroup,'HIC');
    internet_coverage = betainv(remote_quantile,9.57,1.39);
    labsh = betainv(labsh_quantile,7.97,6.87);
    bmi = norminv(bmi_quantile,26.81,1.52);
    Hmax = gaminv(Hmax_quantile, 2.045582620,   1/0.021471378);
else
    country_indices = strcmp(CD.igroup,'LIC') | strcmp(CD.igroup,'LMIC') |...
        strcmp(CD.igroup,'UMIC') | strcmp(CD.igroup,'HIC');
    internet_coverage = 0;
    labsh = 0;
    Hmax = logninv(Hmax_quantile, 2.83235908, 0.87983858);
end

% bmi
bmi = min(max(25.9,bmi), 29.9);
% three outcomes (infection, hospitalisation, death)
% two age groups (20 to 64, 65 plus)
bmi_gradients = [0.0166 0.0518 0.0534; 0.0045 0.026 0.0111];
bmi_intercepts= [0.524 -0.484 -0.531; 0.872 0.254 0.68];
bmi_sigma= [0.00185 0.00501 0.0168; 0.00155 0.0059 0.0121];

bmi_rr_quantile = repmat(unifrnd(0,1,1,3), 2, 1);
bmi_rr = norminv(bmi_rr_quantile, bmi.*bmi_gradients + bmi_intercepts, bmi_sigma);
data.bmi_rr = bmi_rr;
data.bmi = bmi;
data.bmi_rr_quantile = bmi_rr_quantile(1,:);

% self isolation compliance
data.self_isolation_compliance = unifrnd(0,1,1,1);

% contacts
data.B = unifrnd(max(data.B/2-1,0),data.B*2+1);
data.C = unifrnd(max(data.C/2-1,0),data.C*2+1);
data.hospitality_contact_quantile = unifrnd(0,1);
data.hospA2 = logninv(data.hospitality_contact_quantile,log(0.5),0.5);
data.hospA3 = logninv(data.hospitality_contact_quantile,log(1),0.5);
data.hospA4 = logninv(data.hospitality_contact_quantile,log(1.5),0.5);

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
data.NNs(data.NNs==0) = 1;
data.ntot     = size(data.NNs,1);

%CM
nonempind = find(~isnan(CD.CMaa) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,70:325));
defivalue = reshape(randvalue,16,16);
data.CM   = defivalue;

%comm
nonempind = find(~isnan(CD.comm) & country_indices);
mincomm = min(CD.comm(nonempind));
maxcomm = max(CD.comm(nonempind));
data.comm = unifrnd(mincomm,maxcomm,1,1);

%travelA3
nonempind     = find(~isnan(CD.travelA3) & country_indices);
mina3 = min(CD.travelA3(nonempind));
maxa3 = max(CD.travelA3(nonempind));
data.travelA3 =  unifrnd(mina3,maxa3,1,1);

%schoolA1&schoolA2
school_contact_quantile = unifrnd(0,1,1,1);
mina1 = min(CD.schoolA1);
maxa1 = max(CD.schoolA1);
mina2 = min(CD.schoolA2);
maxa2 = max(CD.schoolA2);
data.schoolA1 = unifinv(school_contact_quantile,mina1,maxa1) ; 
data.schoolA2 = unifinv(school_contact_quantile,mina2,maxa2);

% sample_uniform({'schoolA1','schoolA2'},CD,country_indices)

%workp
nonempind  = find(~isnan(CD.workp) & country_indices);
minwp = min(CD.workp(nonempind));
maxwp = max(CD.workp(nonempind));
data.workp = unifrnd(minwp,maxwp,1,1);


%%
%wfh
nonempind = find(~isnan(CD.wfhl1) & country_indices);
mins = min(table2array(CD(nonempind,376:420)));
maxs = max(table2array(CD(nonempind,421:465)));
newprop = unifinv(remote_quantile,mins,maxs);
data.wfh  = [newprop; newprop];

%Tres
nonempind = find(~isnan(CD.Tres) & country_indices & CD.Tres<365);
mint = min(CD.Tres(nonempind));
maxt = max(CD.Tres(nonempind));
data.Tres = unifrnd(mint,maxt,1,1);

%t_tit
nonempind  = find(~isnan(CD.t_tit) & country_indices & CD.t_tit<365);
mint = min(CD.t_tit(nonempind));
maxt = max(CD.t_tit(nonempind));
data.t_tit = unifrnd(mint,maxt,1,1);

%trate
nonempind  = find(~isnan(CD.trate) & country_indices);
mint = min(CD.trate(nonempind));
maxt = max(CD.trate(nonempind));
data.trate = unifrnd(mint,maxt,1,1);

%sdl
nonempind = find(~isnan(CD.sdl) & country_indices);
mins = min(CD.sdl(nonempind));
maxs = max(CD.sdl(nonempind));
data.sdl  = unifrnd(mins,maxs,1,1);

%sdb
nonempind = find(~isnan(CD.sdb) & country_indices);
mins = min(CD.sdb(nonempind));
maxs = max(CD.sdb(nonempind));
data.sdb  = unifrnd(mins,maxs,1,1);

data.Hmax = Hmax; 

%t_vax
nonempind  = find(~isnan(CD.t_vax) & country_indices);
mint = min(CD.t_vax(nonempind));
maxt = max(CD.t_vax(nonempind));
data.t_vax = unifrnd(mint,maxt,1,1);

%arate
nonempind  = find(~isnan(CD.arate) & country_indices);
mint = min(CD.arate(nonempind));
maxt = max(CD.arate(nonempind));
data.arate = unifrnd(mint,maxt,1,1);

%puptake
nonempind    = find(~isnan(CD.puptake) & country_indices);
mint = min(CD.puptake(nonempind));
maxt = max(CD.puptake(nonempind));
data.puptake = unifrnd(mint,maxt,1,1);

%la
nonempind = find(~isnan(CD.la1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,476:493));
defivalue = randvalue;
data.la   = defivalue;


%% change sizes of sectors
% rescale sectors
% get populations
% agInd = 1;
FAAind = 32;
pointiness = 1000;

adultindices = [1:lx,lx+data.adInd];
workingagepop = sum(data.NNs(adultindices));

adult_props = data.NNs(adultindices)./workingagepop;
small_sectors = find(adult_props<1e-3);
resample_sectors = setdiff(adultindices, adultindices(small_sectors));

workingagepop2 = sum(data.NNs(resample_sectors));
adult_props2 = data.NNs(resample_sectors)./workingagepop2;

newvals = gamrnd(adult_props2*pointiness,1);
vals = newvals ./ sum(newvals);
% min(vals)
% max(abs(adult_props2 - vals))
% scatter(adult_props2,vals)

data.NNs(resample_sectors) = workingagepop2 .* vals;


% agPop = data.NNs(agInd);
% FAAPop = data.NNs(FAAind);
% otherindices = setdiff(adultindices,[agInd,FAAind]);
% 
% leftoverpop = workingagepop - agPop - FAAPop;

% popindices = find(~isnan(CD.NNs1)&~isnan(CD.Npop1) & country_indices);
% workers = table2array(CD(popindices,25:69));%number of workers by sector in real country
% adults = sum(table2array(CD(popindices,8:16)),2);

% minFAAworkers = min(workers(:,FAAind) ./ adults);
% minagWorkers = min(workers(:,agInd) ./ adults);
% maxFAAworkers = max(workers(:,FAAind) ./ sum(workers,2));
% maxagWorkers = max(workers(:,agInd) ./ sum(workers,2));
% 
% FAAfrac = unifrnd( minFAAworkers , maxFAAworkers );
% agFrac = unifrnd( minagWorkers, min(maxagWorkers, 1-FAAfrac) );
% leftover_frac = 1 - FAAfrac - agFrac;

% data.NNs(agInd) = agFrac * workingagepop;
% data.NNs(FAAind) = FAAfrac * workingagepop;
% data.NNs(otherindices) = data.NNs(otherindices)./leftoverpop .* workingagepop .* leftover_frac;

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

%% valuations

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

remote_teaching_effectiveness = betarnd(5,5,1,1);
students_affected = 1 - remote_teaching_effectiveness;
mean_annual_income = labsh*gdp/workingagepop;
data.labsh = labsh;

educationloss_all_students = ...
    PV*mean_annual_income*rate_of_return_one_year*students_affected*sum(data.Npop(2:4));
educationloss_per_student = ...
    PV*mean_annual_income*rate_of_return_one_year*students_affected;



defivalue = 0.5454*gdp/sum(data.Npop(2:4));
data.vsy  = educationloss_per_student;
data.educationloss_all_students  = educationloss_all_students;




%% international tourism 

pointiness = 5.93; % for beta distribution
sec_to_international = 3.66; % scales fraction to fraction
international_const = 0.0959; % constant

remaining_international_tourism = logninv(remaining_tourism_quant,-1.39,0.39);
GDP = sum(data.obj);
FAAfrac = data.obj(FAAind)/GDP;
alpha = pointiness * min(sec_to_international * FAAfrac + international_const,1);
beta = pointiness - alpha;
frac_tourism_international = betainv(international_tourism_quant,alpha,beta);
FAAmax = 1 - frac_tourism_international + frac_tourism_international*remaining_international_tourism;
data.frac_tourism_international = frac_tourism_international;

data.unmit = ones(size(data.x_elim));
%%!! set min to one for now until we can allow for increased tourism
data.unmit(FAAind) = min(FAAmax,1);
data.x_elim(FAAind) = min(FAAmax,data.x_elim(FAAind));
data.x_econ(FAAind,:) = min(FAAmax,data.x_econ(FAAind,:));
data.x_schc(FAAind,:) = min(FAAmax,data.x_schc(FAAind,:));



%% population disease parameters

%Contact Matrix
[basic_contact_matrix,data] = p2MakeDs(data,data.NNs,ones(data.lx,1),zeros(1,data.lx));
data.basic_contact_matrix = basic_contact_matrix;

% get covid wt doubling time

dis = population_disease_parameters(data,dis_ref,R0betafun);

data.Td_CWT = dis.Td;



end