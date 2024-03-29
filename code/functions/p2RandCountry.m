% simulate a random country by drawing from distributions and data

% data: struct of general model parameters
% CD: table of country data values
% income_level: string indicating income level (e.g. HIC)
% country_parameter_distributions: pre-specified, named distributions and
% parameters
% social_dist_coefs: table of parameters for social distancing function

% data: struct of general model parameters

function data = p2RandCountry(data,CD,income_level,country_parameter_distributions, social_dist_coefs)

%% start
nSectors = data.nSectors;

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
data.remote_teaching_effectiveness = unifrnd(0,1,1,1);
data.self_isolation_compliance = unifrnd(0,1,1,1);

sdtab_ncol = size(social_dist_coefs,1);
randrow = randi([1 sdtab_ncol],1,1);
data.sd_baseline = social_dist_coefs.baseline(randrow);
data.sd_death_coef = exp(social_dist_coefs.deathcoef(randrow));
data.sd_mandate_coef = exp(social_dist_coefs.mandatecoef(randrow));

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
data.labsh = labsh;


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

% contacts
contacts.B = unifrnd(max(contacts.B/2-1,0),contacts.B*2+1);
contacts.C = unifrnd(max(contacts.C/2-1,0),contacts.C*2+1);

% population by age
nonempind = find(~isnan(CD.CMaa) & ~isnan(CD.Npop1) & country_indices);
[~,idx] = sort(CD.average_contacts(nonempind));
nonempind = nonempind(idx);
demoindex = nonempind(randi(numel(nonempind)));
cols = strmatch('Npop', CD.Properties.VariableNames);
randvalue = table2array(CD(demoindex,cols));
defivalue = 50*10^6*randvalue'/sum(randvalue);
data.Npop = defivalue;

% population by stratum
nonempind = find(~isnan(CD.NNs1) & country_indices);
randindex = nonempind(randi(numel(nonempind)));
colNNs = strmatch('NNs', CD.Properties.VariableNames);
randvalue = table2array(CD(randindex,colNNs));%number of workers by sector in real country
defivalue = randvalue/sum(table2array(CD(randindex,3+[5:13])));%proportion of adult population by sector in real country
defivalue = sum(data.Npop(5:13))*defivalue;%number of workers by sector in artificial country
NNs = [defivalue,data.Npop(1),sum(data.Npop(2:4)),sum(data.Npop(5:13))-sum(defivalue),sum(data.Npop(14:end))]';
data.NNs  = NNs;
data.NNs(data.NNs==0) = 1;
data.nStrata     = size(data.NNs,1);

% contact matrix -- match to pop for now
% nonempind = find(~isnan(CD.CMaa) & country_indices);
% randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(demoindex,70:325));
defivalue = reshape(randvalue,16,16);
contacts.CM   = defivalue;

%workp = number of contacts in workplace
contacts.workp = sample_uniform("workp",CD,country_indices);


%%
%wfh = work from home

nonempind = find(~isnan(CD.wfhl1) & country_indices);
mins = min(table2array(CD(nonempind,strmatch('wfhl', CD.Properties.VariableNames))));
maxs = max(table2array(CD(nonempind,strmatch('wfhu', CD.Properties.VariableNames))));
newprop = unifinv(internet_coverage_quantile,mins,maxs);
data.wfh  = [newprop; newprop];

% date of importation
data.t_import = unifrnd(0,20,1,1);

%Hres = hospital occupancy at response time in origin country
data.Hres = unifinv(data.response_time_quantile, 1,20);

%trate = testing rate
data.trate = sample_uniform("trate",CD,country_indices);

%sdl
data.sdl = sample_uniform("sdl",CD,country_indices);

%sdb
data.sdb = sample_uniform("sdb",CD,country_indices);

%t_vax = time to start vaccine administration
data.t_vax = 1000; 

%arate = vaccine administration rate
data.vaccination_rate_pc = 0.01;%unifrnd(0.5,1.5,1,1)/100;

%puptake = population uptake
data.vaccine_uptake = 0.8; %unifrnd(.4,.8,1,1);

%la = life expectancy
nonempind = find(~isnan(CD.la1) & country_indices);
cols = strmatch('la', CD.Properties.VariableNames);
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,cols));
data.la   = randvalue;


%% change sizes of sectors
% rescale sectors
% get populations
% agInd = 1;
FAAind = 32;
pointiness = 1000;

adultindices = [1:nSectors,nSectors+data.adInd];
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

%% obj: income per worker
nonempind                   = find(~isnan(CD.obj1)&~isnan(CD.NNs1) & country_indices);
randindex                   = nonempind(randi(numel(nonempind)));
% weights = CD.popsum(nonempind)/sum(CD.popsum(nonempind));
% randindex = randsample(nonempind,1,true,weights);
cols1 = strmatch('obj', CD.Properties.VariableNames);
randvalue                   = table2array(CD(randindex,cols1));%gva by sector in real country
defivalue                   = randvalue./table2array(CD(randindex,colNNs));%gva per worker by sector in real country
defivalue(isnan(defivalue)) = 0;
defivalue(isinf(defivalue)) = 0;
defivalue                   = data.NNs(1:45).*defivalue';%gva by sector in artificial country
%%!! need to check food and accomm gva as % of gdp
data.obj                    = defivalue;

%% valuations

discount_rate = 0.03;
gdp = 365*sum(data.obj);
data.gdp = gdp;
data.gdppc = gdp/sum(data.Npop);

%vly = value of a life year
pop_sizes        = [data.Npop(1:17)',sum(data.Npop(18:end))];%length is 18 to match life table
pop_sizes_4 = [pop_sizes(1);sum(pop_sizes(2:4));sum(pop_sizes(5:13));sum(pop_sizes(14:end))];

% map life expectancy to our age groups
life_expectancy  = data.la;
age_map = {1, 2:4, 5:13, 14:length(pop_sizes)};
life_expectancy_4 = zeros(size(age_map));
for k = 1:length(life_expectancy_4)
    index = age_map{k};
    life_expectancy_4(k) = dot(life_expectancy(index),pop_sizes(index))/sum(pop_sizes(index));
end
% get discounted values
discounted_life_expectancy = zeros(size(life_expectancy_4));
for k = 1:length(life_expectancy_4)
    discounted_life_expectancy(k) = sum(1./((1+discount_rate).^(1:life_expectancy_4(k))));
end 

vsl       = 160*gdp/sum(pop_sizes);
discounted_value_of_a_life_year = vsl/(dot(discounted_life_expectancy,pop_sizes_4)/sum(pop_sizes_4));
data.vly  = discounted_value_of_a_life_year;

%% vsy = value of a school year
rate_of_return_one_year = 0.08;
agefracs = data.Npop(2:4)/sum(data.Npop(2:4));
midages = [7 12 17];
working_years = 45;
workstarts = 20 - midages;
workends = working_years + workstarts;

% present value
PV = ((1-(1+discount_rate).^(-workends))/discount_rate - (1-(1+discount_rate).^(-workstarts))/discount_rate)*agefracs;

mean_annual_income = labsh*gdp/workingagepop;
educationloss_all_students = ...
    PV*mean_annual_income*rate_of_return_one_year*sum(data.Npop(2:4));
educationloss_per_student = ...
    PV*mean_annual_income*rate_of_return_one_year;

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
data.contacts = get_basic_contacts(data, contacts);
basic_contact_matrix = p2MakeDs(data,data.NNs,ones(data.nSectors,1),zeros(1,data.nSectors));
data.contacts.basic_contact_matrix = basic_contact_matrix;


end