
%% global variables

income_levels = {'LLMIC','UMIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
vaccination_levels = {'BAU','100 days','BPSV'};
countrytype_levels = {'Origin','Secondary'};

nsamples  = 100;
n_income = numel(income_levels);

synthetic_countries_base = cell(nsamples,length(income_levels));
synthetic_countries = cell(nsamples,length(income_levels));
synthetic_countries_dis = cell(nsamples,length(income_levels));
synthetic_countries_dis_basis = cell(nsamples,1);
synthetic_countries_p2 = cell(nsamples,length(income_levels),length(vaccination_levels));

%% country variables

[CD, country_parameter_distributions] = load_country_data();
data = data_start();
lx = data.lx;

%% disease variables

rng(0);
alldissamples = sample_disease_parameters(nsamples);

R0_to_beta = @(dis) [dis.R0, dis.R0/dis.CI];

% get basic disease profiles
names = fieldnames(alldissamples);
dis = struct;
for i = 1:nsamples
    rng(i);
    for fn = 1:numel(names)
        thisfield = names{fn};
        samples = alldissamples.(thisfield);
        dis.(thisfield) = samples(i,:);
    end
    synthetic_countries_dis_basis{i} = dis;
end

%% countries by disease

% 2: get all countries
% get covid doubling time
% 3: get beta and R0, where beta depends on R0
betas = zeros(nsamples,n_income);
R0s = zeros(nsamples,n_income);
for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    for il = 1:n_income
        rng(i);
        income_level = income_levels{il};
        
        % country data. random samples
        ldata1     = p2RandCountry(data,CD,income_level,country_parameter_distributions);
        synthetic_countries_base{i,il} = ldata1;
        
        for vl = 1:length(vaccination_levels)
            [ldata,dis2,p2] = p2Params(ldata1,dis,R0_to_beta,vl);
            synthetic_countries_p2{i,il,vl} = p2;
        end
        synthetic_countries_dis{i,il} = dis2;
        synthetic_countries{i,il}     = ldata;
    end
end

%% set up simulation
columnnames = {'VLY','VSY',...
    'School_contacts','School_age','Working_age','Elders',...
    'Unemployment_rate','GDP','Labour_share','Work_contacts',...
    'Hospitality_contacts','Hospital_capacity','Test_rate',...
    'Hospital_response','Response_time',...
    'Social_distancing_max','Social_distancing_rate','Self_isolation_compliance','Candidate_infectees',...
    'R0','beta','Mean_IHR',...
    'Mean_HFR','Mean_IFR','Probability_symptomatic','Latent_period',...
    'Asymptomatic_period','Symptomatic_period','Time_to_hospitalisation','Time_to_discharge',...
    'Time_to_death','Agriculture','Food_sector','International_tourism',...
    'Remote_teaching_effectiveness','Importation_time','Hospital_occupancy_at_response',...
    'Internet','Start_vaccination','Vaccination_rate','Vaccine_uptake',...
    'Cost','Deaths','School','GDP_loss'};
inputs    = zeros(nsamples,length(columnnames)-4);
outputs   = zeros(nsamples,4);

%% simulate

hts = zeros(nsamples,n_income);

for il = 1:n_income
    income_level = income_levels{il};
    for ms = 1:length(strategies)
        strategy = strategies{ms};
        for vl = 1:length(vaccination_levels)
            for ct = 1:length(countrytype_levels)
                ht = 730*ones(1,nsamples);
                countrytype = countrytype_levels(ct);
                parfor i = 1:nsamples
                    dis2 = synthetic_countries_dis{i,il};
                    p2 = synthetic_countries_p2{i,il,vl};
                    ldata= synthetic_countries{i,il};

                    [rdata, xoptim] = get_strategy_design(ldata,strategy,countrytype,p2);

                    sec         = nan(1,4);
                    int = 5;
                    try
%                         disp([il ms vl ct i]);
                        [~,f,g] = p2Run(rdata,dis2,strategy,int,xoptim,p2);
                        ht(i) = f(find(f(:,1) > p2.Tres,1),3);
                        [cost,~]    = p2Cost(ldata,dis2,p2,g);
                        sec(1)      = sum(cost([3,6,7:10],:),'all');
                        sec(2)      = sum(cost([3],:),'all');
                        sec(3)      = sum(cost([6],:),'all');
                        sec(4)      = sum(cost([7:10],:),'all');
                    catch
                        disp(strcat('VOI_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'.NA'))
                        disp([il ms vl ct i]);
                    end
                    if any(sec<0)
                        disp(strcat('VOI_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'.0'))
                        disp(i)
                    end

                    outputs(i,:) = sec;

                    gdp = sum(ldata.obj);
                    popsize = sum(ldata.Npop);
                    working_age = popsize - sum(ldata.NNs([46,47,49]));
                    unemployment_rate = ldata.NNs(48)/working_age;
                    contacts = ldata.contacts;
                    inputs(i,:)  = [ldata.vly ldata.vsy ...
                        contacts.schoolA2 ldata.NNs(47)/popsize working_age/popsize ldata.NNs(49)/popsize...
                        unemployment_rate ldata.gdp ldata.labsh contacts.workrel ...
                        contacts.hospitality_frac ldata.Hmax ldata.trate  ...
                        ldata.Hres p2.Tres ldata.sdl ldata.sdb ldata.self_isolation_compliance dis2.CI ...
                        dis2.R0 dis2.beta mean(dis2.ihr) ...
                        mean(dis2.hfr) mean(dis2.ifr) dis2.ps dis2.Tlat ...
                        dis2.Tay dis2.Tsr dis2.Tsh dis2.Threc ...
                        dis2.Thd ldata.obj([1 32])'/gdp ldata.frac_tourism_international ...
                        ldata.remote_teaching_effectiveness ldata.t_import ht(i)...
                        ldata.remote_quantile ldata.t_vax ldata.arate ldata.puptake];
                end
                
                if vl==1 & ms==1 & ct==1
                    hts(:,il) = ht;
                end
                T                          = array2table([inputs,outputs]);
                T.Properties.VariableNames = columnnames;

                writetable(T,strcat('results/VOI_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'_',string(countrytype),'.csv'));
            end
        end
    end
end

