
addpath('functions');

%% global variables

income_levels = {'LLMIC','UMIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
vaccination_levels = [365, 100];
bpsv_levels = [0, 1];

nsamples  = 2048;
n_income = numel(income_levels);

synthetic_countries = cell(nsamples,length(income_levels));
synthetic_countries_dis = cell(nsamples,length(income_levels));
synthetic_countries_dis_basis = cell(nsamples,1);
synthetic_countries_p2 = cell(nsamples,length(income_levels),length(vaccination_levels),length(bpsv_levels));

%% country variables

[CD, country_parameter_distributions, social_dist_coefs] = load_country_data();
data = data_start();

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

for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    for il = 1:n_income
        % shuffle countries
        rng(il+n_income*(i-1));
        income_level = income_levels{il};
        % country data. random samples
        ldata1     = p2RandCountry(data,CD,income_level,country_parameter_distributions,social_dist_coefs);
        % get combined country and disease parameters
        dis1 = population_disease_parameters(ldata1,dis,R0_to_beta);
        % sample a response time
        ldata1.response_time = get_response_time(ldata1,dis1,ldata1.Hres);
        % save temporarily:
        synthetic_countries_dis{i,il} = dis1;
        synthetic_countries{i,il}     = ldata1;
    end
end
clear synthetic_countries_dis_basis

%% shuffle response times - does not depend on income level
origin_countries = randi(length(income_levels),1,nsamples);
for i = 1:nsamples
    for il = 1:n_income
        synthetic_countries{i,il}.response_time = synthetic_countries{i,origin_countries(i)}.response_time;
    end
end

%% complete p2, dis and data structs

for i = 1:nsamples
    for il = 1:n_income
        dis1 = synthetic_countries_dis{i,il};
        ldata1 = synthetic_countries{i,il};
        % get p2 parameters: depend on vaccine scenario
        for vl = 1:length(vaccination_levels)
            for bl = 1:length(bpsv_levels)
                [ldata,dis2,p2] = p2Params(ldata1,dis1,vaccination_levels(vl),bpsv_levels(bl));
                synthetic_countries_p2{i,il,vl,bl} = p2;
            end
        end
        synthetic_countries_dis{i,il} = dis2;
        synthetic_countries{i,il}     = ldata;
    end
end



%% save inputs
%%!! assuming there are no changes between scenarios and strategies, only
%%income levels
inputcolumnnames = tabulate_inputs(ldata,p2,dis2);
inputs    = zeros(nsamples,length(inputcolumnnames));
for il = 1:n_income
    income_level = income_levels{il};
    for i = 1:nsamples
        %% load stored objects
        dis2 = synthetic_countries_dis{i,il};
        ldata = synthetic_countries{i,il};                   
        p2 = synthetic_countries_p2{i,il,1,1};

        [~, vals] = tabulate_inputs(ldata,p2,dis2);
        inputs(i,:)  = vals;
    end
    T = array2table(inputs);
    T.Properties.VariableNames = inputcolumnnames;
    writetable(T,strcat('results/inputs_',string(income_level),'.csv'));
end


%% set up simulation

outputcolumnnames = {'Mitigated_deaths','End_mitigation','End_simulation', 'Remaining_susceptible','End_hosp','Exit_wave',...
    'Deaths','Cost','dYLLs','School','GDP_loss'};
columnnames = [outputcolumnnames ];
outputs   = zeros(nsamples,length(outputcolumnnames));

%% simulate

for il = 1:n_income
    income_level = income_levels{il};
    for ms = 1:length(strategies)
        strategy = strategies{ms};
        for vl = 1:length(vaccination_levels)
            for bl = 1:length(bpsv_levels)
                parfor i = 1:nsamples
                    %% load stored objects
                    dis2 = synthetic_countries_dis{i,il};
                    ldata = synthetic_countries{i,il};                   
                    p2 = synthetic_countries_p2{i,il,vl,bl};
                    try
                        %% run model
                        [~,returned] = p2Run(ldata,dis2,strategy,p2);
        %                         figure('Position', [100 100 400 300]); plot(returned.Tout,returned.Htot)
                        
                        %% outputs: costs
                        costs    = p2Cost(ldata,dis2,p2,returned);
                        sec         = nan(1,4);
                        sec(2)      = sum(costs.value_dYLL); % dylls
                        sec(3)      = sum(costs.value_SYL); % school
                        sec(4)      = sum(costs.GDP_lost);  % gdp
                        sec(1)      = sum(sec(2:4)); % cost
                        total_deaths = returned.deathtot(end);
                        
                        %% store some intermediate values
                        % store final time point
                        endsim = max(returned.Tout);
                        % store time mitigation ends
                        endmit = returned.isequence(end,1);
                        % get index of endmit time
                        [~,exitwave] = min(abs(returned.Tout-endmit));
                        % get fraction of deaths that happen after
                        % mitigation ends
                        mitdeaths = returned.deathtot(exitwave);
                        exitwavefrac = 1-mitdeaths/returned.deathtot(end);
                        % number in hospital at end
                        endhosp = returned.Htot(end);
                        % total still susceptible at end of simulation
                        endsusc = returned.Stotal(end)/returned.Stotal(1);
                        % hospital occupancy at response time
                        ht = returned.Htot(find(returned.Tout > p2.Tres,1));
                        
                        %% store outputs
                        outputs(i,:) = [mitdeaths endmit endsim endsusc endhosp exitwavefrac total_deaths sec];
                        
                        if any(sec<0)
                            disp(strcat(string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'_',string(bpsv_levels(bl)),'_',string(i),' 0'))
                            disp(i)
                        end
                    catch
                        disp(strcat(string(income_level),'_',string(strategy),'_',string(vaccination_levels(vl)),'_',string(bpsv_levels(bl)),'_',string(i),' NA'))
                        disp([il ms vl bl i]);
                    end
                end   
                % write results
                T = array2table(outputs);
                T.Properties.VariableNames = columnnames;
                writetable(T,strcat('results/outputs_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'_',string(bpsv_levels(bl)),'.csv'));
            end
        end
    end
end

