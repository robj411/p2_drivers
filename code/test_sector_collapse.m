
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
CI_cell = cell(length(income_levels),1);
contacts_cell = cell(length(income_levels),1);
adinds = [1:45,48];
for il = 1:n_income
    CI_cell{il} = zeros(nsamples,7);
    contacts_cell{il} = zeros(nsamples,length(adinds));
end
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
        synthetic_countries{i,il}     = ldata1;
        
        % get collapsed CI
        NNs = ldata1.NNs;
        zs = zeros(size(NNs));
        [CI4, strata_weights, contacts46] = get_candidate_infectees4(length(NNs), dis1, NNs,zs, zs, 0, 0, NNs, ldata1.contacts.basic_contact_matrix);
        
        
        meancontacts = dot(strata_weights, contacts46);
        varcontacts = sum((meancontacts - contacts46).^2 .* strata_weights);
        work_frac = ldata1.contacts.work_frac;
        gini = sum(sum(abs(contacts46 - contacts46').*strata_weights.*strata_weights'))/2/meancontacts;
        
        employmentrate = ldata1.employmentrate;
        
        contacts_cell{il}(i,:) = contacts46;
        CI_cell{il}(i,:) = [CI4, dis1.CI, ...
            work_frac, varcontacts, ...
            employmentrate, ldata1.pupil_teacher_ratio,...
            gini];
    end
end

    scatter(CI_cell{il}(:,1)./CI_cell{il}(:,2),CI_cell{il}(:,7),'.')
    scatter((CI_cell{il}(:,3)),CI_cell{il}(:,5),'.')
colcol = 7;
y = (CI_cell{il}(:,1))./CI_cell{il}(:,2);
for il = 1:n_income
    figure;
    tiledlayout(3,3,"Padding","tight")
    nexttile
    histogram(CI_cell{il}(:,5))
    xlabel('Employment')
    nexttile
    scatter(CI_cell{il}(:,4),CI_cell{il}(:,5),'.','CData',CI_cell{il}(:,colcol))
    nexttile
    scatter(y,CI_cell{il}(:,5),'.','CData',CI_cell{il}(:,colcol))
    nexttile(5)
    histogram(CI_cell{il}(:,4))
    xlabel('Variance')
    nexttile(6)
    scatter(y,CI_cell{il}(:,4),'.','CData',CI_cell{il}(:,colcol))
    ylim([0,100])
    nexttile(9)
    histogram(y)
    xlabel('Ratio')
    saveas(gcf,sprintf('../figures/sector_collapse_%s',string(income_levels(il))),'jpg');

end


close all


find(CI_cell{il}(:,4)<.5)
find(CI_cell{il}(:,7)<.05)
find(CI_cell{il}(:,4)>2000)
synthetic_countries{2020,il}.contacts.sectorcontacts
synthetic_countries{1006,il}.contacts.sectorcontacts
synthetic_countries{2020,il}.pupil_teacher_ratio
synthetic_countries{1006,il}.pupil_teacher_ratio
synthetic_countries{1006,il}.NNs
synthetic_countries{2020,il}.NNs


logit = @(x) log(x./(1-x));
il=3;
y = (CI_cell{il}(:,2)-CI_cell{il}(:,1))./CI_cell{il}(:,2);
X = array2table([y,CI_cell{il}(:,[3,5:7]),log(CI_cell{il}(:,4))]);

md = fitrgam(X,'Var1')
scatter(y,predict(md,X))


for i = [94, 1186]
    NNs = synthetic_countries{i,il}.NNs(adinds);
    [contacts iy] = sort(contacts_cell{il}(i,:));
    p = NNs(iy)/sum(NNs);

    x = cumsum([0 p']); % start of bar
    y = zeros(length(x),1);
    dx = diff(x); % width of bar
    pcontacts = cumsum(contacts);
    dy = cumsum(contacts);


    figure, hold on
    for ii=1:length(dx)
        rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)])
    end
    axis([0 1 0 max(dy)])
end

