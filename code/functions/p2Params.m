% get parameters that characterise p2 (pandemic preparedness)
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% vaccine_day: 100 or 365, depending on SARS-X vaccine scenario
% bpsv: 0 or 1, depending on BPSV scenario (present or absent)
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% p2: struct of p2 intervention parameters

function [data,dis,p2] = p2Params(data,dis,vaccine_day,bpsv)

%% PREPAREDNESS PARAMETERS

p2 = struct;

% social distancing
% p2.sdl   = data.sdl;                        %Social Distancing Asymptote
% p2.sdb   = data.sdb;                        %Social Distancing Steepness

% response time and testing start time set by global alert
p2.Tres = data.response_time;
p2.t_tit = data.response_time;

% testing impacts
self_isolation_compliance = data.self_isolation_compliance;
p2.trate = data.trate;                      %Test-Isolate-Trace Rate
data = rmfield(data,'trate');
% amount averted by isolating needs to account for discovery rule
% time taken to test
time_to_test = 2;
% time taken to isolate after symptoms
time_to_isolate = 1;
% total amount avertible if isolating after symptom onset
frac_by_symptoms = max(1 - dis.frac_presymptomatic - time_to_isolate/dis.Tsr, 0);
% remaining amount can be averted by testing
frac_by_testing = 1 - frac_by_symptoms;
% sym: assume isolation on symptoms, so infectiousness averted is frac that
% is presymptomatic, unless test taken in presymptomatic period
% if people isolate when symptomatic, a fraction of infectiousness is
% averted always
frac_averted = frac_by_symptoms;
p2.frac_sym_infectiousness_averted = self_isolation_compliance * frac_averted * 1;
% if people test when presymptomatic, a fraction of presymptomatic
% infectiousness is averted when testing
frac_averted = max(frac_by_testing - time_to_test/dis.Tsr, 0);
p2.frac_presym_infectiousness_averted = self_isolation_compliance * frac_averted;
% assume they test after one day.
p2.frac_asym_infectiousness_averted = self_isolation_compliance * min(1,1-time_to_test./dis.Tay);
p2.time_to_test = time_to_test;

% Hospital Capacity
Hmax  = data.Hmax*sum(data.Npop)/10^5;   %Hospital Capacity
data = rmfield(data,'Hmax');
p2.hosp_release_trigger   = max(1,0.25*Hmax);%lower threshold can't be less than 1 occupant
p2.Hmax  = max(4*p2.hosp_release_trigger,Hmax);
% p2.SHmax = 2*p2.Hmax;

% stopping criteria
p2.hosp_final_threshold = 100;
p2.final_doubling_time_threshold = 30;

%% Vaccine Uptake

t_vax    = data.t_vax;                      %Vaccine Administration Time

vaccination_rate_pc    = data.vaccination_rate_pc;  %Vaccine Administration Rate
vaccine_uptake  = data.vaccine_uptake;                    %Vaccine Uptake

Npop    = data.Npop;
Npop4 = data.Npop4;
over14frac = Npop(4)/sum(Npop4(2));

p2.t_vax2 = p2.Tres + 7 + vaccine_day;  

if bpsv == 1
    t_vax = p2.Tres;
end

vaccination_rate = vaccination_rate_pc*sum(Npop);

t_vax2 = p2.t_vax2;
t_bspv = t_vax2 - t_vax;

uptake  = vaccine_uptake*[0 over14frac 1 1]; 


%Vaccine Administration Rate
t_ages     = min((uptake.*Npop4)/vaccination_rate,Inf);%vaccination_rate may be 0

if bpsv==1
    % if primer starts before booster, there are t_bspv days of primer.
    % these people must then be boosted.
    if t_bspv < t_ages(4)
        t_ages = [t_ages(4)-t_bspv, t_ages(4), t_ages(3), t_ages(2)];
        p2.group_order = [4,4,3,2];
    else
        t_ages = [t_ages(4), t_bspv-t_ages(4), t_ages(4), t_ages(3), t_ages(2)];
        p2.group_order = [4,0,4,3,2];
    end
else
    % if booster starts before primer, just work down the ages
    t_ages     = [t_ages(4),t_ages(3),t_ages(2)];
    p2.group_order = [4,3,2];
end    
tpoints    = cumsum([min(t_vax, p2.t_vax2), t_ages]);
p2.tpoints = tpoints;
p2.arate = vaccination_rate;


%% COST PARAMETERS

% update age groups and indices for number of groups
lihr = length(dis.ihr);
Npop     = data.Npop';
Npop     = [Npop(1:(lihr-1)),sum(Npop(lihr:end))];%last age range for disease is 80+
ageindex = data.ageindex;
ageindex{4} = min(ageindex{4}):lihr;

% get life expectancy for age groups
life_expectancy         = data.la(1:lihr);
life_expectancy(lihr) = dot(data.la(lihr:end),[data.Npop(lihr),sum(data.Npop((lihr+1):end))]) / sum(data.Npop(lihr:end)+1e-10);

% get population-weighted average ifr
weighted_ifr       = Npop.*dis.ifr;
% use weight to compute life expectancy lost per death
life_years_lost_per_death         = arrayfun(@(x) dot(life_expectancy(x{1}),weighted_ifr(x{1}))/sum(weighted_ifr(x{1})), ageindex);
% apply discounting
% discount_rate = 0.03;
% discounted_life_years_per_death = zeros(size(life_years_lost_per_death));
% for k = 1:length(life_years_lost_per_death)
%     discounted_life_years_per_death(k) = sum(1./((1+discount_rate).^[1:life_years_lost_per_death(k)]));
% end  
data.lgh   = [repmat(life_years_lost_per_death(data.adInd),1,45),life_years_lost_per_death];

end