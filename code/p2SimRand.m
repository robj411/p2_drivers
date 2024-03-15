
addpath('functions');

%% global variables

income_levels = {'LLMIC','UMIC','HIC'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
vaccination_levels = [365, 100];
bpsv_levels = [0, 1];

nsamples  = 128;
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
        rng(i);
        income_level = income_levels{il};
        % country data. random samples
        ldata1     = p2RandCountry(data,CD,income_level,country_parameter_distributions,social_dist_coefs);
        % get combined country and disease parameters
        dis1 = population_disease_parameters(ldata1,dis,R0_to_beta);
        % sample a response time
        ldata1.rts = get_response_time(ldata1,dis1,ldata1.Hres);
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
clear synthetic_countries_dis_basis

%% set up simulation

inputcolumnnames = {'VLY','VSY',...
    'School_contacts','School_age','Working_age','Elders',...
    'Unemployment_rate','GDP','Labour_share','Work_contacts',...
    'Hospitality_contacts',...
    'Hospital_capacity','Test_rate','Hospital_response','Response_time',...
    'Vaccination_rate','Vaccine_uptake',...
    'Social_distancing_baseline','Social_distancing_death','Social_distancing_mandate','Fraction_infectiousness_averted','Candidate_infectees',...
    'Agriculture','Food_sector','International_tourism',...
    'Remote_teaching_effectiveness','Importation_time',...
    'Remote_quantile','Hospital_occupancy_at_response',... 
    'Doubling_time','Generation_time',...
    'R0','beta','Frac_presymptomatic',...
    'Mean_IHR','Mean_HFR','Mean_IFR',...
    'Max_IHR','Max_HFR','Max_IFR',...
    'Probability_symptomatic','Latent_period',...
    'Asymptomatic_period','Symptomatic_period','Time_to_hospitalisation','Time_to_discharge',...
    'Time_to_death',...
    'End_mitigation','End_simulation', 'Remaining_susceptible','Exit_wave'};%,...
    %'Elimination_R0','School_R0_1','School_R0_2','Econ_R0_1','Econ_R0_2'};
outputcolumnnames = {'Cost','dYLLs','School','GDP_loss','Deaths'};
columnnames = [inputcolumnnames outputcolumnnames ];
inputs    = zeros(nsamples,length(inputcolumnnames));
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
                        endmit = iseq(end,1);
                        % get index of endmit time
                        [~,exitwave] = min(abs(returned.Tout-endmit));
                        % get fraction of deaths that happen after
                        % mitigation ends
                        exitwavefrac = 1-returned.deathtot(exitwave)/returned.deathtot(end);
                        % total still susceptible at end of simulation
                        endsusc = returned.Stotal(end)/returned.Stotal(1);
                        % hospital occupancy at response time
                        ht = returned.Htot(find(returned.Tout > p2.Tres,1));
                        
                        %% store inputs
                        gdp = sum(ldata.obj);
                        popsize = sum(ldata.Npop);
                        Npop = [ldata.Npop(1:(length(dis2.ihr)-1)) ; sum(ldata.Npop(length(dis2.ihr):length(ldata.Npop)))];
                        meanihr = dis2.ihr * Npop / popsize;
                        meanifr = dis2.ifr * Npop / popsize;
                        meanhfr = dis2.hfr * Npop / popsize;
                        maxihr = max(dis2.ihr);
                        maxifr = max(dis2.ifr);
                        maxhfr = max(dis2.hfr);
                        working_age = popsize - sum(ldata.NNs([46,47,49]));
                        unemployment_rate = ldata.NNs(48)/working_age;
                        contacts = ldata.contacts;
                        inputs(i,:)  = [ldata.vly ldata.vsy ...
                            contacts.schoolA2 ldata.NNs(47)/popsize working_age/popsize ldata.NNs(49)/popsize...
                            unemployment_rate ldata.gdp ldata.labsh contacts.workrel ...
                            contacts.hospitality_frac...
                            p2.Hmax ldata.trate ldata.Hres p2.Tres...
                            ldata.vaccination_rate_pc ldata.vaccine_uptake ...
                            ldata.sd_baseline ldata.sd_death_coef ldata.sd_mandate_coef p2.frac_sym_infectiousness_averted dis2.CI ...
                            ldata.obj([1 32])'/gdp ldata.frac_tourism_international ...
                            ldata.remote_teaching_effectiveness ldata.t_import...
                            ldata.remote_quantile ht...
                            dis2.Td dis2.generation_time ...
                            dis2.R0 dis2.beta dis2.frac_presymptomatic ...
                            meanihr meanhfr meanifr ...
                            maxihr maxhfr maxifr ...
                            dis2.ps dis2.Tlat ...
                            dis2.Tay dis2.Tsr dis2.Tsh dis2.Threc ...
                            dis2.Thd ...
                            endmit endsim endsusc exitwavefrac];% ...
%                             dis2.R0s];
                        outputs(i,:) = [sec total_deaths];

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
                T = array2table([inputs outputs]);
                T.Properties.VariableNames = columnnames;
                writetable(T,strcat('results/VOI_',string(strategy),'_',string(income_level),'_',string(vaccination_levels(vl)),'_',string(bpsv_levels(bl)),'.csv'));
            end
        end
    end
end

