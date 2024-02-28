% get parameters that characterise p2 (pandemic preparedness)

% dis: struct of pathogen parameters
% p2: struct of p2 intervention parameters
% data: struct of general model parameters
% vaccine_day: 100 or 326, depending on SARS-X vaccine scenario
% bpsv: 0 or 1, depending on BPSV scenario (present or absent)

% dis: struct of pathogen parameters
% p2: struct of p2 intervention parameters
% data: struct of general model parameters

function [data,dis,p2] = p2Params(data,dis,vaccine_day,bpsv)

%% PREPAREDNESS PARAMETERS

p2 = struct;

% social distancing
p2.sdl   = data.sdl;                        %Social Distancing Asymptote
p2.sdb   = data.sdb;                        %Social Distancing Steepness

% response time and testing start time set by global alert
p2.Tres = data.rts;
p2.t_tit = data.rts;

% testing impacts
p2.self_isolation_compliance = data.self_isolation_compliance;
p2.trate = data.trate;                      %Test-Isolate-Trace Rate
p2.dur    = 1;
p2.frac_sym_infectiousness_remaining = min(1,p2.dur./(dis.Tsr+dis.Tsh)*2);
p2.frac_asym_infectiousness_remaining = min(1,p2.dur./dis.Tay);

% Hospital Capacity
p2.Hmax  = data.Hmax*sum(data.Npop)/10^5;   %Hospital Capacity
p2.thl   = max(1,0.25*p2.Hmax);%lower threshold can't be less than 1 occupant
p2.Hmax  = max(4*p2.thl,p2.Hmax);
p2.SHmax = 2*p2.Hmax;

% stopping criteria
p2.hosp_final_threshold = 1000;
p2.final_doubling_time_threshold = 30;

%% Vaccine Uptake

t_vax    = data.t_vax;                      %Vaccine Administration Time

if bpsv==0 && vaccine_day==326
    vaccination_rate_pc = 0.5/100;
    vaccine_uptake  = .4;
else
    vaccination_rate_pc    = data.vaccination_rate_pc;  %Vaccine Administration Rate
    vaccine_uptake  = data.vaccine_uptake;                    %Vaccine Uptake
end

Npop    = data.Npop;
NNage   = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
over14frac = Npop(4)/sum(Npop(2:4));

p2.t_vax2 = p2.Tres + 7 + vaccine_day;  

if bpsv == 1
    t_vax = p2.Tres;
end

vaccination_rate = vaccination_rate_pc*sum(data.Npop);

t_vax2 = p2.t_vax2;
t_bspv = t_vax2 - t_vax;

uptake  = vaccine_uptake*[0 over14frac 1 1]; 


%Vaccine Administration Rate
t_ages     = min((uptake.*NNage)/vaccination_rate,Inf);%vaccination_rate may be 0

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

pop_size         = [data.Npop(1:16)',sum(data.Npop(17:end))];%length is 17 to match ifr
life_expectancy         = [data.la(1:16),...
              dot(data.la(17:end),[data.Npop(17),sum(data.Npop(18:end))])/sum(data.Npop(17:end)+1e-10)];
weighted_ifr       = pop_size.*dis.ifr;
life_years_lost_per_death         = [dot(life_expectancy(1),weighted_ifr(1))/sum(weighted_ifr(1)),...
              dot(life_expectancy(2:4),weighted_ifr(2:4))/sum(weighted_ifr(2:4)),...
              dot(life_expectancy(5:13),weighted_ifr(5:13))/sum(weighted_ifr(5:13)),...
              dot(life_expectancy(14:end),weighted_ifr(14:end))/sum(weighted_ifr(14:end))];
discounted_life_years_per_death = zeros(size(life_years_lost_per_death));
for k = 1:length(life_years_lost_per_death)
    discounted_life_years_per_death(k) = sum(1./((1+0.03).^[1:life_years_lost_per_death(k)]));
end  
data.lgh   = [repmat(discounted_life_years_per_death(3),1,45),discounted_life_years_per_death];

end