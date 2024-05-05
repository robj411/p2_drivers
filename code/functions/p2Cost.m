% computes socio-economic costs from epi model outputs

% p2: struct of p2 intervention parameters
% data: struct of general model parameters
% dis: struct of pathogen parameters
% returnobject: output from epi model

% cost: 

function costs = p2Cost(data,dis,p2,returnobject)

t  = returnobject.Tout;
nSectors = length(data.obj);
nStrata = nSectors+4;

workers = returnobject.workers;
homeworkers = returnobject.homeworkers;
Iamat = returnobject.Iamat;
Ismat = returnobject.Ismat;
hospmat = returnobject.hospmat;
deathmat = returnobject.deathmat;
selfisolation = returnobject.selfisolation;

costs = struct;

days_per_infectious_day = 2;
days_per_infectious_day_asym = days_per_infectious_day./dis.Tay.*dis.Tsr;

p3 = selfisolation.p3;
p4 = selfisolation.p4; % does not matter: symptomatic isolate for a fixed number of days based on symptom onset. testing therefore affects their infectiousness but not their time spent at home.

self_isolation_compliance = data.self_isolation_compliance;
frac_cases_found = p3 / p2.frac_asym_infectiousness_averted;
frac_isolating = frac_cases_found * self_isolation_compliance; % maybe should not scale by compliance because why do you test if you do not comply?

%% VLYL

deaths    = deathmat;
costs.deaths = deaths(end,:);

lyl       = deaths(end,:).*data.lgh;
costs.life_years = lyl;

vlyl      = lyl*data.vly;
costs.value_dYLL = vlyl;

% ccost_t(:,1:ln) = deaths.*data.lgh.*data.vly;

%% VSYL

Stu              = nSectors+2;
students         = data.NNs(Stu);
% cost(4,nSectors+[1,2]) = students;

%Student Supply
days_per_infectious_day_sym = days_per_infectious_day./dis.Ts(Stu).*dis.Tsr;
prob_hosp = dis.ph(Stu).*[1 1-dis.hv1 1-dis.hv2];
isoasym       = sum(reshape(Iamat(:,Stu,:),[],3),2).* days_per_infectious_day_asym .* frac_isolating; 

symstudents = reshape(Ismat(:,Stu,:),[],3);
sym_no_hosp       = sum(symstudents.*(1-prob_hosp),2).* self_isolation_compliance .* days_per_infectious_day_sym; 
sym_to_hosp       = sum(symstudents.*prob_hosp,2).* self_isolation_compliance + hospmat(:,Stu) ; 
% deaths       = deathmat(:,Stu) ; 
isosym          = sym_no_hosp + sym_to_hosp;% + deaths;%numbers of students

%Student Demand
closure = (1-returnobject.workers(:,data.EdInd)) .* (1 - data.remote_teaching_effectiveness);
not_learning        = closure.*students + (1-closure).*isosym + (1-2*closure).*isoasym;
not_learning_int     = trapz(t,not_learning)/365;%= (diff(t)'*presl)/365;

costs.value_SYL = not_learning_int * data.vsy;


%% SGDPL

notEd = [1:(data.EdInd-1),(data.EdInd+1):nSectors];
worker_numbers = data.NNs(notEd);

% labour supply
isoasym       = sum(Iamat(:,notEd,:),3).* days_per_infectious_day_asym .* frac_isolating; 

days_per_infectious_day_sym = days_per_infectious_day./dis.Ts(notEd).*dis.Tsr;
prob_hosp = dis.ph(notEd).*[1 1-dis.hv1 1-dis.hv2];

sym_no_hosp_total = zeros(size(Ismat(:,notEd,:)));
sym_to_hosp_total = zeros(size(Ismat(:,notEd,:)));
for i = 1:size(sym_no_hosp_total,1)
    for j = 1:size(prob_hosp,2)
        sym_no_hosp_total(i,:,j) = Ismat(i,notEd,j) .* (1-prob_hosp(:,j)');
        sym_to_hosp_total(i,:,j) = Ismat(i,notEd,j) .* prob_hosp(:,j)';
    end
end
sym_no_hosp       = sum(sym_no_hosp_total,3).* self_isolation_compliance .* days_per_infectious_day_sym'; 
sym_to_hosp       = sum(sym_to_hosp_total,3).* self_isolation_compliance + hospmat(:,notEd) + deathmat(:,notEd) ; 
% deaths       = deathmat(:,Stu) ; 
isosym          = sym_no_hosp + sym_to_hosp;% + deaths;

% Labour demand
hw            = homeworkers(:,notEd);
x             = workers(:,notEd);

workers_absent = isosym + (1-hw).*isoasym;
worker_presence        = x.*(1 - workers_absent./worker_numbers');
worker_presence_int     = trapz(t,worker_presence);

GDP_in = sum(worker_presence_int .* data.obj(notEd)');
max_GDP = (max(t)-min(t))*sum(data.obj(notEd));
costs.GDP_lost = max_GDP - GDP_in;


end