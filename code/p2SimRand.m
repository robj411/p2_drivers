
%% global variables

income_levels = {'LLMIC','MIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};

nsamples  = 4096;
n_income = numel(income_levels);

synthetic_countries_base = cell(nsamples,length(income_levels));
synthetic_countries = cell(nsamples,length(income_levels));
synthetic_countries_dis = cell(nsamples,length(income_levels));
synthetic_countries_dis_basis = cell(nsamples,1);
synthetic_countries_p2 = cell(nsamples,length(income_levels));

%% country variables

CD        = readtable('../data/country_data.csv');
data = data_start();
lx = data.lx;

%% disease variables

dis_ref = get_dis_params('Covid Wildtype'); 

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
CIs = zeros(nsamples,n_income);
other_parameters = zeros(nsamples*n_income,7);
for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    for il = 1:n_income
        rng(i);
        income_level = income_levels{il};
        
        % country data. defines wnorm and Td_CWT
        ldata1     = p2RandCountry(data,CD,income_level,dis_ref,R0_to_beta);
        synthetic_countries_base{i,il} = ldata1;
        
        % get beta and R0
        dis2 = population_disease_parameters(ldata1,dis,R0_to_beta);
        betas(i,il) = dis2.beta;
        R0s(i,il) = dis2.R0;
        CIs(i,il) = dis2.CI;
        other_parameters(i+(il-1)*nsamples,:) = [dis2.ps, dis2.sig1, dis2.sig2, dis2.g1, dis2.Tsr, dis2.Tsh, mean(dis2.ihr)];
    end
end

%% generate new R0 using beta

% 4: get relationship between R0 and beta.
% Resample R0.
% use covid doubling time

sigmoid = @(x) 1./(1+exp(-x));
logit = @(x) log(x./(1-x));

beta_to_R0 = @(dis) [dis.beta.*dis.CI, dis.beta];

betameans = mean(betas,2);

betas2 = zeros(nsamples,n_income);
R0s2 = zeros(nsamples,n_income);
CIs2 = zeros(nsamples,n_income);
for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    dis.beta = betameans(i);
    for il = 1:n_income
        ldata0 = synthetic_countries_base{i,il};
        
        [ldata,dis2,p2] = p2Params(ldata0,dis,beta_to_R0);

        synthetic_countries_dis{i,il} = dis2;
        synthetic_countries_p2{i,il} = p2;
        synthetic_countries{i,il}     = ldata;
        betas2(i,il) = dis2.beta;
        R0s2(i,il) = dis2.R0;
        CIs2(i,il) = dis2.CI;
    end
end
scatter(R0s(:),R0s2(:))

%% which samples cannot be contained

schoolRs = zeros(nsamples,n_income);
econRs = zeros(nsamples,n_income);
p3s = zeros(nsamples,n_income);
for il = 1:n_income
     for i = 1:nsamples
        dis2 = synthetic_countries_dis{i,il};
        p2 = synthetic_countries_p2{i,il};
        ldata = synthetic_countries{i,il};
        Rs = test_lockdown(p2,ldata,dis2);
        econRs(i,il) = Rs(1);
        schoolRs(i,il) = Rs(2);
        p3s(i,il) = Rs(3);
     end
end

[econrows,econcols] = find(econRs>1);
[schoolrows,schoolcols] = find(schoolRs>1);
disp(length(schoolrows))

%% set up simulation
columnnames = {'Econ_LD_R','School_LD_R',...
    'BMI','BMI_infection','BMI_hospitalisation','BMI_death',...
    'School_contacts','School_age','Working_age','Elders',...
    'Unemployment_rate','GDP','Labour_share','Work_contacts',...
    'Hospitality_contacts','Community_contacts','Hospital_capacity','Test_rate','Test_start',...
    'Response_time','Social_distancing_max','Social_distancing_rate','Self_isolation_compliance','Candidate_infectees',...
    'R0','beta','Mean_IHR',...
    'Mean_HFR','Mean_IFR','Probability_symptomatic','Latent_period',...
    'Asymptomatic_period','Symptomatic_period','Time_to_hospitalisation','Time_to_discharge',...
    'Time_to_death','Agriculture','Food_sector','International_tourism',...
    'Internet',...
    'Cost','Deaths','School','GDP_loss'};
inputs    = zeros(nsamples,length(columnnames)-6);
outputs   = zeros(nsamples,4);

%% simulate

for il = 1:n_income
    income_level = income_levels{il};
    for ms = 1:length(strategies)
        strategy = strategies{ms};

        parfor i = 1:nsamples
            dis2 = synthetic_countries_dis{i,il};
            p2 = synthetic_countries_p2{i,il};
            ldata = synthetic_countries{i,il};

            int = 5;
            xoptim = 0;
            if strcmp(strategy,'Elimination')
                xoptim      = [ones(1*lx,1);ldata.x_econ(:,2);ldata.x_elim(:,1);ones(2*lx,1)];
                ldata.hw    = [zeros(1,lx);ldata.wfh(2,:);ldata.wfh(1,:);zeros(2,lx)];
                ldata.imand = [2];
                ldata.inext = [2,2,3,2,5];
            elseif strcmp(strategy,'Economic Closures')
                xoptim      = [ones(1*lx,1);ldata.unmit;ldata.x_econ(:,2);ldata.x_econ(:,1);ones(lx,1)];
                ldata.hw    = [zeros(1,lx);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,lx)];
                ldata.imand = [3];
                ldata.inext = [2,3,3,4,5];
            elseif strcmp(strategy,'School Closures')
                xoptim      = [ones(1*lx,1);ldata.unmit;ldata.x_schc(:,2);ldata.x_schc(:,1);ones(lx,1)];
                ldata.hw    = [zeros(1,lx);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,lx)];
                ldata.imand = [3];
                ldata.inext = [2,3,3,4,5];
            elseif strcmp(strategy,'No Closures')
                xoptim      = [ones(1*lx,1); repmat(ldata.unmit,4,1)];
                ldata.hw    = [zeros(5,lx)];
                ldata.imand = [10];
                ldata.inext = [2,2,5];
            else
                error('Unknown Mitigation Strategy!');
            end

            sec         = nan(1,4);

            try
                [ldata,~,g] = p2Run(ldata,dis2,strategy,int,xoptim,p2);
                [cost,~]    = p2Cost(ldata,dis2,p2,g);
                sec(1)      = sum(cost([3,6,7:10],:),'all');
                sec(2)      = sum(cost([3],:),'all');
                sec(3)      = sum(cost([6],:),'all');
                sec(4)      = sum(cost([7:10],:),'all');
            catch
                disp(strcat('VOI_',string(strategy),'_',string(income_level),'.NA'))
                disp(i);
        %             dis2
            end
            if any(sec<0)
                disp(strcat('VOI_',string(strategy),'_',string(income_level),'.0'))
                disp(i)
            end

            outputs(i,:) = sec;

            gdp = sum(ldata.obj);
            popsize = sum(ldata.Npop);
            working_age = popsize - sum(ldata.NNs([46,47,49]));
            unemployment_rate = ldata.NNs(48)/working_age;
            inputs(i,:)  = [ldata.bmi ldata.bmi_rr_quantile...
                ldata.schoolA2 ldata.NNs(47)/popsize working_age/popsize ldata.NNs(49)/popsize...
                unemployment_rate ldata.gdp ldata.labsh ldata.workp ...
                ldata.hospitality_contact_quantile ldata.comm ldata.Hmax ldata.trate ldata.t_tit ...
                ldata.Tres ldata.sdl ldata.sdb ldata.self_isolation_compliance dis2.CI ...
                dis2.R0 dis2.beta mean(dis2.ihr) ...
                mean(dis2.hfr) mean(dis2.ifr) dis2.ps dis2.Tlat ...
                dis2.Tay dis2.Tsr dis2.Tsh dis2.Threc ...
                dis2.Thd ldata.obj([1 32])'/gdp ldata.frac_tourism_international ...
                ldata.remote_quantile];

        end

        T                          = array2table([econRs(:,il),schoolRs(:,il),inputs,outputs]);
        T.Properties.VariableNames = columnnames;

        writetable(T,strcat('results/VOI_',string(strategy),'_',string(income_level),'.csv'));

    end
end

