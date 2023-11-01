function data = p2RandCountry(data,CD,income_level,country_parameter_distributions)

%% start
lx = data.lx;

contacts = data.contacts;

%% generate quantiles
internet_coverage_quantile = unifrnd(0,1);
labsh_quantile = unifrnd(0,1);
bmi_quantile = unifrnd(0,1);
Hmax_quantile = unifrnd(0,1);
pt_quantile = unifrnd(0,1);
schoolA1_frac_quantile = unifrnd(0,1,1,1);
schoolA2_frac_quantile = unifrnd(0,1,1,1);
remaining_international_tourism_quantile = unifrnd(0,1);

international_tourism_quant = unifrnd(0,1);
contacts.hospitality_frac = betarnd(10,40,1,1);

data.remote_quantile = internet_coverage_quantile;
data.response_time_quantile = unifrnd(0,1);

%% values from distributions
pindices = find(strcmp(country_parameter_distributions.igroup,income_level) | ...
    strcmp(country_parameter_distributions.igroup,'all') & ...
    ~strcmp(country_parameter_distributions.distribution,'NA'));
cpd = country_parameter_distributions(pindices,:);
for i = 1:size(cpd,1)
    varname = cpd.parameter_name{i};
    expression = strcat('=',cpd.distribution{i},'(',...
        cpd.parameter_name{i},'_quantile,',...
        num2str(cpd.Parameter_1(i)),',',...
        num2str(cpd.Parameter_2(i)),');');
    eval([varname, expression]);
end

contacts.pt = pt;
contacts.schoolA1_frac = schoolA1_frac;
contacts.schoolA2_frac = schoolA2_frac;
data.Hmax = Hmax; 

if strcmp(income_level,'LLMIC')
    country_indices = strcmp(CD.igroup,'LIC') | strcmp(CD.igroup,'LMIC');
elseif strcmp(income_level,'UMIC')
    country_indices = strcmp(CD.igroup,'UMIC');
elseif strcmp(income_level,'HIC')
    country_indices = strcmp(CD.igroup,'HIC');
end

%% map to parameters
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
contacts.B = unifrnd(max(contacts.B/2-1,0),contacts.B*2+1);
contacts.C = unifrnd(max(contacts.C/2-1,0),contacts.C*2+1);

%Npop
nonempind = find(~isnan(CD.CMaa) & ~isnan(CD.Npop1) & country_indices);
[~,idx] = sort(CD.average_contacts(nonempind));
nonempind = nonempind(idx);
%weights = unifrnd(0,1,1,1) + cumsum(CD.popsum(nonempind)/sum(CD.popsum(nonempind)));
%demoindex = find(weights>1,1) ; %randsample(nonempind,1,true,weights);
demoindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(demoindex,4:24));
defivalue = 50*10^6*randvalue'/sum(randvalue);
data.Npop = defivalue;

%NNs
nonempind = find(~isnan(CD.NNs1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
%weights = CD.popsum(nonempind)/sum(CD.popsum(nonempind));
%randindex = randsample(nonempind,1,true,weights);
randvalue = table2array(CD(randindex,25:69));%number of workers by sector in real country
defivalue = randvalue/sum(table2array(CD(randindex,3+[5:13])));%proportion of adult population by sector in real country
defivalue = sum(data.Npop(5:13))*defivalue;%number of workers by sector in artificial country
defivalue = [defivalue,data.Npop(1),sum(data.Npop(2:4)),sum(data.Npop(5:13))-sum(defivalue),sum(data.Npop(14:end))]';
data.NNs  = defivalue;
data.NNs(data.NNs==0) = 1;
data.ntot     = size(data.NNs,1);

%CM -- match to pop for now
% nonempind = find(~isnan(CD.CMaa) & country_indices);
% randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(demoindex,70:325));
defivalue = reshape(randvalue,16,16);
contacts.CM   = defivalue;

%comm
% nonempind = find(~isnan(CD.comm) & country_indices);
% mincomm = min(CD.comm(nonempind));
% maxcomm = max(CD.comm(nonempind));
% contacts.comm = unifrnd(mincomm,maxcomm,1,1);

%travelA3
% nonempind     = find(~isnan(CD.travelA3) & country_indices);
% mina3 = min(CD.travelA3(nonempind));
% maxa3 = max(CD.travelA3(nonempind));
% normalise France values: 18% travel is pt, and pt has 2.5% of contacts,
% which is 0.555
% contacts.travelA3 =  pt/0.18*0.555 ; %unifrnd(mina3,maxa3,1,1);

%schoolA1&schoolA2
% mina1 = min(CD.schoolA1);
% maxa1 = max(CD.schoolA1);
% mina2 = min(CD.schoolA2);
% maxa2 = max(CD.schoolA2);
% contacts.schoolA1 = unifinv(school_contact_quantile,mina1,maxa1) ; 
% contacts.schoolA2 = unifinv(school_contact_quantile,mina2,maxa2);

% sample_uniform({'schoolA1','schoolA2'},CD,country_indices)

%workp
nonempind  = find(~isnan(CD.workp) & country_indices);
minwp = min(CD.workp(nonempind));
maxwp = max(CD.workp(nonempind));
contacts.workp = unifrnd(minwp,maxwp,1,1);


%%
%wfh
nonempind = find(~isnan(CD.wfhl1) & country_indices);
mins = min(table2array(CD(nonempind,376:420)));
maxs = max(table2array(CD(nonempind,421:465)));
newprop = unifinv(internet_coverage_quantile,mins,maxs);
data.wfh  = [newprop; newprop];

% date of importation
data.t_import = unifrnd(10,20,1,1);

%Tres
% nonempind = find(~isnan(CD.Tres) & country_indices & CD.Tres<365);
% mint = 30; %min(CD.Tres(nonempind));
% maxt = 90; %max(CD.Tres(nonempind));
% data.Tres = unifinv(data.response_time_quantile, mint,maxt);
data.Hres = unifinv(data.response_time_quantile, 1,20);

%t_tit
% nonempind  = find(~isnan(CD.t_tit) & country_indices & CD.t_tit<365);
% mint = min(CD.t_tit(nonempind));
% maxt = max(CD.t_tit(nonempind));
% data.t_tit = data.Tres; %unifrnd(mint,maxt,1,1);

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

%t_vax
nonempind  = find(~isnan(CD.t_vax) & country_indices);
mint = min(CD.t_vax(nonempind));
maxt = max(CD.t_vax(nonempind));
data.t_vax = unifinv(data.response_time_quantile, mint,maxt);

%arate
% nonempind  = find(~isnan(CD.arate) & country_indices);
% mint = min(CD.arate(nonempind));
% maxt = max(CD.arate(nonempind));
data.arate = 10^5 * 0.005; %unifrnd(mint,maxt,1,1);

%puptake
% nonempind    = find(~isnan(CD.puptake) & country_indices);
% mint = min(CD.puptake(nonempind));
% maxt = max(CD.puptake(nonempind));
% data.puptake = unifrnd(mint,maxt,1,1);
data.puptake = 0.80;

%la
nonempind = find(~isnan(CD.la1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
%weights = CD.popsum(nonempind)/sum(CD.popsum(nonempind));
%randindex = randsample(nonempind,1,true,weights);
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
% randindex                   = nonempind(randi(numel(nonempind)));
weights = CD.popsum(nonempind)/sum(CD.popsum(nonempind));
randindex = randsample(nonempind,1,true,weights);
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
for k = 1:length(lg)
    lgh(k) = sum(1./((1+0.03).^(1:lg(k))));
end 
vsl       = 160*gdp/sum(na);
defivalue = vsl/(dot(lgh,[na(1);sum(na(2:4));sum(na(5:13));sum(na(14:end))])/sum(na));
data.vly  = defivalue;
data.gdp = gdp;
data.gdppc = gdp/sum(data.Npop);

%vsy
agefracs = data.Npop(2:4)/sum(data.Npop(2:4));
midages = [7 12 17];
workstarts = 20 - midages;
working_years = 45;
workends = working_years + workstarts;
discount_rate = 0.03;
rate_of_return_one_year = 0.08;
PV = ((1-(1+discount_rate).^(-workends))/discount_rate - (1-(1+discount_rate).^(-workstarts))/discount_rate)*agefracs;

data.remote_teaching_effectiveness = unifrnd(0,1,1,1);
mean_annual_income = labsh*gdp/workingagepop;
data.labsh = labsh;

educationloss_all_students = ...
    PV*mean_annual_income*rate_of_return_one_year*sum(data.Npop(2:4));
educationloss_per_student = ...
    PV*mean_annual_income*rate_of_return_one_year;



defivalue = 0.5454*gdp/sum(data.Npop(2:4));
data.vsy  = educationloss_per_student;
data.educationloss_all_students  = educationloss_all_students;




%% international tourism 

pindex = find(strcmp(country_parameter_distributions.parameter_name,'tourism_pointiness'));
pindex2 = find(strcmp(country_parameter_distributions.parameter_name,'sec_to_international'));

pointiness = country_parameter_distributions.Parameter_1(pindex); % for beta distribution
sec_to_international = country_parameter_distributions.Parameter_1(pindex2); % scales fraction to fraction
international_const = country_parameter_distributions.Parameter_2(pindex2); % constant

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



%% contacts

%Contact Matrix
data.contacts = contacts;
data = get_basic_contacts(data);
basic_contact_matrix = p2MakeDs(data,data.NNs,ones(data.lx,1),zeros(1,data.lx));
data.contacts.basic_contact_matrix = basic_contact_matrix;

% get covid wt doubling time

% dis = population_disease_parameters(data,dis_ref,R0betafun);

% data.Td_CWT = dis.Td;




end