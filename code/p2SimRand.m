
income_levels = {'LLMIC','MIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};

load('../country_mats/Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
fields    = fieldnames(data);
ikeep     = [6,7,8,13,14,16,17,18];
data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
data.adInd = 3;
data.lx = length(data.B);
data.tvec = [-75 365];
data.alp = 1;
data.EdInd    = 41;%education sector index
data.HospInd  = [32,43,44];%hospitality sector indices


dis_ref = get_dis_params('Covid Wildtype'); 

CD        = readtable('../data/country_data.csv');
nsamples  = 1000;

alldissamples = sample_disease_parameters(nsamples);

synthetic_countries_base = cell(nsamples,length(income_levels));
synthetic_countries = cell(nsamples,length(income_levels));
synthetic_countries_dis = cell(nsamples,length(income_levels));
synthetic_countries_dis_basis = cell(nsamples,1);
synthetic_countries_p2 = cell(nsamples,length(income_levels));

n_income = numel(income_levels);

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

% 2: get all countries
% get covid doubling time
% 3: get beta and R0, where beta depends on R0
betas = [];
R0s = [];
CIs = [];
for i = 1:nsamples
    rng(i);
    dis = synthetic_countries_dis_basis{i};
    for il = 1:n_income
        income_level = income_levels{il};
        
        % country data. defines wnorm and Td_CWT
        ldata1     = p2RandCountry(data,CD,income_level,dis_ref,R0_to_beta);
        synthetic_countries_base{i,il} = ldata1;
        
        % get beta and R0
        dis2 = population_disease_parameters(ldata1,dis,R0_to_beta);
        betas(i+(il-1)*nsamples) = dis2.beta;
        R0s(i+(il-1)*nsamples) = dis2.R0;
        CIs(i+(il-1)*nsamples) = dis2.CI;
        
    end
end
figure; tiledlayout(1,3)
nexttile
scatter(R0s(1:nsamples),CIs(1:nsamples),30,'filled')
nexttile
scatter(CIs(1:nsamples),betas(1:nsamples),30,'filled')
nexttile
scatter(R0s(1:nsamples),betas(1:nsamples),30,'filled')


% 4: get relationship between R0 and beta.
% Resample R0.
% use covid doubling time
betamodel = glmfit(R0s,betas);
bmodelintercept = betamodel(1);
bmodelgrad = betamodel(2);
beta_to_R0 = @(dis) [(bmodelgrad*dis.R0 + bmodelintercept)*dis.CI , bmodelgrad*dis.R0 + bmodelintercept];

for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    for il = 1:n_income
        ldata0 = synthetic_countries_base{i,il};
        
        [ldata,dis2,p2] = p2Params(ldata0,dis,beta_to_R0);

        synthetic_countries_dis{i,il} = dis2;
        synthetic_countries_p2{i,il} = p2;
        synthetic_countries{i,il}     = ldata;
        betas(i+(il-1)*nsamples) = dis2.beta;
        R0s(i+(il-1)*nsamples) = dis2.R0;
        CIs(i+(il-1)*nsamples) = dis2.CI;
    end
end
figure; tiledlayout(1,3)
nexttile
scatter(R0s,CIs,30,'filled')
nexttile
scatter(CIs,betas,30,'filled')
nexttile
scatter(R0s,betas,30,'filled')



columnnames = {'School_contacts','School_age','Working_age','Elders','Population_size','Unemployment_rate',...
    'GDP','Labour_share','workp',...
    'Hospital_capacity','Test_rate','Test_start','Response_time',...
    'Social_distancing_min','Social_distancing_rate','Candidate_infectees','R0','beta','Max_IHR','Max_IFR',...
    'Agriculture','Food_sector','International_tourism','Internet',...
    'Cost','Deaths','School','GDP_loss'};
inputs    = zeros(nsamples,length(columnnames)-4);
outputs   = zeros(nsamples,4);



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

        %         inputs(i,:)  = [i,NaN,NaN,ldata.Npop',ldata.NNs(1:45)',...
        %                         ldata.CM(:)',ldata.comm,ldata.travelA3,ldata.schoolA1,ldata.schoolA2,ldata.workp,...
        %                         ldata.obj',ldata.wfh(1,:),ldata.wfh(2,:),...
        %                         ldata.t_vax,ldata.arate,ldata.puptake,ldata.Hmax,ldata.t_tit,ldata.trate,ldata.Tres,ldata.sdl,ldata.sdb,...
        %                         NaN,ldata.la];
            outputs(i,:) = sec;

            gdp = sum(ldata.obj);
            popsize = sum(ldata.Npop);
            working_age = popsize - sum(ldata.NNs([46,47,49]));
            unemployment_rate = ldata.NNs(48)/working_age;
            inputs(i,:)  = [ldata.schoolA2 ldata.NNs(47)/popsize working_age/popsize ldata.NNs(49)/popsize popsize...
                unemployment_rate ldata.gdp ldata.labsh ldata.workp ldata.Hmax ldata.trate ldata.t_tit ldata.Tres ...
                ldata.sdl ldata.sdb dis2.CI dis2.R0 dis2.beta max(dis2.ihr) max(dis2.ifr) ...
                ldata.obj([1 32])'/gdp ldata.frac_tourism_international ldata.remote_quantile];
        %         sectorprops(i,:) = ldata.NNs([1:45,48])' ./ sum(ldata.NNs([1:45,48]));
        %         R0samples(i) = ;

        end

        T                          = array2table([inputs,outputs]);
        T.Properties.VariableNames = columnnames;

        %     'SEC','VLYL','VSYL','GDPL'};
        writetable(T,strcat('results/VOI_',string(strategy),'_',string(income_level),'.csv'));
        %plots = p2Plot(data,f,p2,g,cost,ccost_t,sec(1),inp1,inp2,inp3);


    end
end

