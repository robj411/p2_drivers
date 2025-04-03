
addpath('functions');

%% global variables

income_levels = {'LLMIC','UMIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
nsamples  = 2048;
n_income = numel(income_levels);

%% country variables

[CD, country_parameter_distributions, utr_coefs] = load_country_data();
data = data_start();
nScen = length(data.scenarios);

synthetic_countries = cell(nsamples,length(income_levels));
synthetic_countries_dis = cell(nsamples,length(income_levels));
synthetic_countries_dis_basis = cell(nsamples,1);
synthetic_countries_p2 = cell(nsamples,length(income_levels),nScen);

%% disease variables

% prep self isolation compliance to depend on R0
R0_quant = zeros(nsamples,1);

rng(0);
[alldissamples, R0_dist] = sample_disease_parameters(nsamples);

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
    % get R0 quantile
    R0_quant(i) = cdf(R0_dist,dis.R0);
end

% generate correlated self-isolation variables
si_quant = correlate_random_var(R0_quant, 0.7);

%% countries by disease

for i = 1:nsamples
    dis = synthetic_countries_dis_basis{i};
    for il = 1:n_income
        % shuffle countries
        rng(il+n_income*(i-1));
        income_level = income_levels{il};
        % country data. random samples
        ldata1     = p2RandCountry(data,CD,income_level,country_parameter_distributions,utr_coefs);
        % get combined country and disease parameters
        [dis1, ldata1] = population_disease_parameters(ldata1,dis,R0_to_beta,R0_dist);
        % convert to beta random variable
        ldata1.self_isolation_compliance = betainv(si_quant(i), 5,5);
        % save temporarily:
        synthetic_countries_dis{i,il} = dis1;
        synthetic_countries{i,il}     = ldata1;
    end
end
% clear synthetic_countries_dis_basis

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
        for sl = 1:nScen
            [ldata,dis2,p2] = p2Params(ldata1,dis1,sl);
            synthetic_countries_p2{i,il,sl} = p2;
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

outputcolumnnames = {'State_changes','Breach_before','Breach_after','Mitigated_deaths','End_mitigation','End_simulation', 'Remaining_susceptible','End_hosp','Exit_wave',...
    'Deaths1','Deaths2','Deaths3','Deaths4','Deaths','Cost','YLL','School','GDP_loss'};
columnnames = [outputcolumnnames ];
outputs   = zeros(nsamples,length(outputcolumnnames));

%% simulate

for il = 1:n_income
    income_level = income_levels{il};
    for sl = 1:nScen
        allcosts = zeros(nsamples,length(strategies));
        for ms = 1:length(strategies)        
            strategy = strategies{ms};
            parfor i = 1:nsamples
                %% load stored objects
                dis2 = synthetic_countries_dis{i,il};
                ldata = synthetic_countries{i,il};                   
                p2 = synthetic_countries_p2{i,il,sl};
                try
                    %% run model
                    [dataout,returned] = p2Run(ldata,dis2,strategy,p2);
    %                         figure('Position', [100 100 400 300]); plot(returned.Tout,returned.Htot)

                    %% outputs: costs
                    costs    = p2Cost(ldata,dis2,p2,returned);
                                        
                    sec         = nan(1,4);
                    sec(2)      = sum(costs.value_YLL); % ylls
                    sec(3)      = sum(costs.value_SYL); % school
                    sec(4)      = sum(costs.GDP_lost);  % gdp
                    sec(1)      = sum(sec(2:4)); % cost
                    
                    total_deaths = returned.deathtot(end);
                    deaths1 = returned.death1(end);
                    deaths2 = returned.death2(end);
                    deaths3 = returned.death3(end);
                    deaths4 = returned.death4(end);

                    %% store some intermediate values
                    % store final time point
                    endsim = max(returned.Tout);
                    % store time mitigation ends
                    endmit = returned.isequence(find(returned.isequence(:,2)>4,1),1);
                    % get index of endmit time
                    [~,exitwave] = min(abs(returned.Tout-endmit));
                    % mitigation ends
                    mitdeaths = returned.deathtot(exitwave);
                    % get fraction of deaths that happen after
                    exitwavefrac = 1-mitdeaths/returned.deathtot(end);
                    % number in hospital at end
                    endhosp = returned.Htot(end);
                    % total still susceptible at end of simulation
                    endsusc = returned.Stotal(end)/returned.Stotal(1);
                    % hospital occupancy at response time
                    ht = returned.Htot(find(returned.Tout > p2.Tres,1));
                    % was hospital capacity breached in the response period
                    breach_before = max(returned.Htot(1:exitwave)) - p2.Hmax;
                    % was hospital capacity breached in the exit wave
                    breach_after = max(returned.Htot(exitwave:end)) - p2.Hmax;

                    %% store outputs
                    outputs(i,:) = [size(returned.isequence,1) breach_before breach_after mitdeaths endmit endsim endsusc endhosp exitwavefrac deaths1 deaths2 deaths3 deaths4 total_deaths sec];

                    if any(sec<0)
                        disp(strcat(string(strategy),'_',string(income_level),'_scen',string(sl),'_',string(i),' 0'))
                        disp(i)
                    end
                catch
                    disp(strcat(string(income_level),'_',string(strategy),'_scen',string(sl),'_',string(i),' NA'))
                    disp([il sl ms i]);
                end
            end   
            % write results
            T = array2table(outputs);
            T.Properties.VariableNames = columnnames;
            writetable(T,strcat('results/outputs_',string(strategy),'_',string(income_level),'_scen',string(sl),'.csv'));
            allcosts(:,ms) = T.Cost;
        end
        disp([il sl]);
        [~,id] = min(allcosts');
        disp(tabulate(id));
    end
end

% !\Progra~1\R\R-4.4.1\bin\x64\Rscript cepi_voi.R
!\Progra~1\R\R-4.4.1\bin\x64\Rscript exceedance_probabilities.R

