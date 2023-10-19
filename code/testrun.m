il=3;ms=4;vl=1;

[il ms vl]


%% global variables

income_levels = {'LLMIC','UMIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
vaccination_levels = {'BAU','100 days'};
countrytype_levels = {'Spillover','Secondary'};

nsamples  = 20;
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
        rng(i);
        income_level = income_levels{il};
        
        % country data. defines wnorm and Td_CWT
        ldata1     = p2RandCountry(data,CD,income_level);
        synthetic_countries_base{i,il} = ldata1;
        
        % get beta and R0
        dis2 = population_disease_parameters(ldata1,dis,R0_to_beta);
        betas(i,il) = dis2.beta;
        R0s(i,il) = dis2.R0;
end

%% generate new R0 using beta

% 4: get relationship between R0 and beta.
% Resample R0.
% use covid doubling time

beta_to_R0 = @(dis) [dis.beta.*dis.CI, dis.beta];

betameans = betas(:,il);%prod(betas').^(1/length(income_levels));

betas2 = zeros(nsamples,n_income);
R0s2 = zeros(nsamples,n_income);
rts = zeros(nsamples,n_income);
for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    dis.beta = betameans(i);
        ldata0 = synthetic_countries_base{i,il};
            [ldata,dis2,p2] = p2Params(ldata0,dis,beta_to_R0,vl);
            synthetic_countries_p2{i,il,vl} = p2;
        synthetic_countries_dis{i,il} = dis2;
        synthetic_countries{i,il}     = ldata;
        rts(i,il) = get_response_time(ldata,dis2,20);
        betas2(i,il) = dis2.beta;
        R0s2(i,il) = dis2.R0;
        if dis2.Td==Inf
            disp([i il]);
            break;
        end
end
% scatter(R0s(:),R0s2(:))

%% set up simulation
columnnames = {'VLY','VSY',...
    'BMI','BMI_infection','BMI_hospitalisation','BMI_death',...
    'School_contacts','School_age','Working_age','Elders',...
    'Unemployment_rate','GDP','Labour_share','Work_contacts',...
    'Hospitality_contacts','Community_contacts','Hospital_capacity','Test_rate',...
    'Response_time_quantile','Actual_response_time',...
    'Social_distancing_max','Social_distancing_rate','Self_isolation_compliance','Candidate_infectees',...
    'R0','beta','Mean_IHR',...
    'Mean_HFR','Mean_IFR','Probability_symptomatic','Latent_period',...
    'Asymptomatic_period','Symptomatic_period','Time_to_hospitalisation','Time_to_discharge',...
    'Time_to_death','Agriculture','Food_sector','International_tourism',...
    'Remote_teaching_effectiveness',...
    'Internet','Start_vaccination','Vaccination_rate','Vaccine_uptake',...
    'Cost','Deaths','School','GDP_loss'};
inputs    = zeros(nsamples,length(columnnames)-4);
outputs   = zeros(nsamples,4);

%% simulate
iseqtab = [];
hts = zeros(nsamples,n_income);

    income_level = income_levels{il};
        strategy = strategies{ms};
        for i = 1:nsamples
            for ct = 1:length(countrytype_levels)
                ht = 730*ones(1,nsamples);
                countrytype = countrytype_levels(ct);
            
                dis2 = synthetic_countries_dis{i,il};
                p2 = synthetic_countries_p2{i,il,vl};
                ldata= synthetic_countries{i,il};
                
                [rdata, xoptim] = get_strategy_design(ldata,strategy,countrytype,p2);
                disp([rdata.t_import p2.Tres])

                sec         = nan(1,4);
                int = 5;
                try
                    [~,f,g,isequence] = p2Run(rdata,dis2,strategy,int,xoptim,p2);
                    iseqtab = [iseqtab; [isequence [rdata.t_import ct i].*ones(size(isequence,1),3)]];
                    if any(f(:,3) > ldata.Hres)
                        ht(i) = f(find(f(:,3) > ldata.Hres,1),1);
                    end
                    [cost,~]    = p2Cost(ldata,dis2,p2,g);
                    sec(1)      = sum(cost([3,6,7:10],:),'all');
                    sec(2)      = sum(cost([3],:),'all');
                    sec(3)      = sum(cost([6],:),'all');
                    sec(4)      = sum(cost([7:10],:),'all');
                catch
                    disp(strcat('VOI_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'.NA'))
                    disp([il ms vl i]);
            %             dis2
                end
                if any(sec<0)
                    disp(strcat('VOI_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'.0'))
                    disp(i)
                end
                
                
            end
            if vl==1 & ms==1
                hts(:,il) = ht;
            end
        end
writematrix(iseqtab,strcat('results/iseq.csv'));
