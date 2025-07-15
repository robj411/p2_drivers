% get parameters that characterise p2 (pandemic preparedness)
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% scenario: integer from 1 to 5 indicating which scenario
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% p2: struct of p2 intervention parameters

function [data,dis,p2] = p2Params(data,dis,scenario)

    %% PREPAREDNESS PARAMETERS

    p2 = struct;

    % response time and testing start time set by global alert
    p2.Tres = data.response_time;
    p2.t_tit = data.response_time;

    % testing impacts
    self_isolation_compliance = data.self_isolation_compliance;
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

    % stopping criteria
    p2.hosp_final_threshold = 1000;
    p2.final_doubling_time_threshold = 30;

    %% Vaccine Uptake

    scentab = data.scenarios{scenario};
    bpsv = sum(scentab(:,1))>0;

    Npop    = data.Npop;
    Npop4 = data.Npop4;
    over14frac = Npop(4)/sum(Npop4(2));

    vaccine_uptake  = data.vaccine_uptake;                    %Vaccine Uptake
    uptake  = vaccine_uptake*[0 over14frac 1 1]; 
    doses4 = Npop4 .* uptake;
    doses = sum(doses4);

    bpsv_sched = scentab(:,1);
    sarsx_sched = scentab(:,2);
    cumdoses = cumsum(sarsx_sched*50000000);
    endrollout = find(cumdoses > doses,1);
    first_bpsv_nonzero = find(bpsv_sched>0,1);

    vaccine_day = find(sarsx_sched>0,1);
    p2.t_vax2 = p2.Tres + 7 + vaccine_day;  

    if bpsv == 1
        cum_bpsv = cumsum(bpsv_sched*50000000);
        keep_indices = find(diff([0; cum_bpsv])>0);
        interptime = interp1([0; cum_bpsv(keep_indices)],[0; keep_indices],doses4(4));
        % find(cum_bpsv>doses4(4),1)
        t_bpsv = min(interptime, find(bpsv_sched~=0,1,'last')) - first_bpsv_nonzero + 1;
        t_vax = p2.Tres + first_bpsv_nonzero;
    else
        t_vax = 10*365;
    end

    t_vax2 = p2.t_vax2;
    start_gap = t_vax2 - t_vax;

    cum_sarsx = cumsum(sarsx_sched*50000000);
    keep_indices = find(diff([0; cum_sarsx])>0);
    interptime = @(x) interp1([0; cum_sarsx(keep_indices)],[0; keep_indices],x);

    t_ages = [];
    t_ages(4) = interptime(doses4(4)) - vaccine_day;
    if max(cumdoses) <= doses
        if max(cumdoses) <= sum(doses4(3:4))
            t_ages(3) = length(cumdoses) - sum(t_ages(4)) - vaccine_day;
            t_ages(2) = 0;
        else
            t_ages(3) = interptime(sum(doses4(3:4))) - t_ages(4) - vaccine_day;
            t_ages(2) = length(cumdoses) - sum(t_ages(3:4)) - vaccine_day;
        end
    else
        t_ages(3) = interptime(sum(doses4(3:4))) - t_ages(4) - vaccine_day;
        t_ages(2) = interptime(sum(doses4(2:4))) - sum(t_ages(3:4)) - vaccine_day;
    end

    % vaccination_rate = vaccination_rate_pc*sum(Npop);
    % vaccination_rate_pc    = data.vaccination_rate_pc;  %Vaccine Administration Rate


    %Vaccine Administration Rate
    % t_ages     = min((uptake.*Npop4)/vaccination_rate,Inf);%vaccination_rate may be 0

    if bpsv==1
        % if primer starts before booster, there are t_bspv days of primer.
        % these people must then be boosted.
        if start_gap < t_bpsv
            t_ages = [start_gap, t_ages(4), t_ages(3), t_ages(2)];
            p2.group_order = [4,4,3,2];
        else
            t_ages = [t_bpsv, start_gap-t_bpsv, t_ages(4), t_ages(3), t_ages(2)];
            p2.group_order = [4,0,4,3,2];
        end
    else
        % if booster starts before primer, just work down the ages
        t_ages     = [t_ages(4),t_ages(3),t_ages(2)];
        p2.group_order = [4,3,2];
    end    
    tpoints    = cumsum([min(t_vax, p2.t_vax2), t_ages]);
    p2.tpoints = tpoints;
    p2.sarsx_per_day = sarsx_sched(vaccine_day:end)*sum(Npop);
    p2.bpsv_per_day = bpsv_sched(first_bpsv_nonzero:end)*sum(Npop);

    

    %Vaccination Rollout by Sector
    NNbar = data.NNs;
    nSectors = data.nSectors;
    sumWorkingAge = sum(NNbar([1:nSectors,nSectors+3]));
    NNnext              = NNbar;
    NNnext(nSectors+[1,2])    = 1;
    NNnext([1:nSectors,nSectors+3]) = NNnext([1:nSectors,nSectors+3])/sumWorkingAge;
    NNnext(end)         = 1;
    p2.NNnext = NNnext;

    %% copy over parameters

    p2.trate = data.trate;                      %Test-Isolate-Trace Rate
    data = rmfield(data,'trate');

    % Hospital Capacity
    Hmax  = data.Hmax*sum(Npop)/10^5;   %Hospital Capacity
    data = rmfield(data,'Hmax');
    p2.hosp_release_trigger   = max(1,0.25*Hmax);%lower threshold can't be less than 1 occupant
    p2.Hmax  = max(4*p2.hosp_release_trigger,Hmax);
    % p2.SHmax = 2*p2.Hmax;


    %% COST PARAMETERS

    % update age groups and indices for number of groups
    lihr = length(dis.ihr);
    Npop     = data.Npop';
    Npop     = [Npop(1:(lihr-1)),sum(Npop(lihr:end))];%last age range for disease is 80+
    ageindex = data.ageindex;
    ageindex{4} = min(ageindex{4}):lihr;

    % get life expectancy for age groups
    life_expectancy         = data.life_expectancy(1:lihr);
    life_expectancy(lihr) = dot(data.life_expectancy(lihr:end),[data.Npop(lihr),sum(data.Npop((lihr+1):end))]) / sum(data.Npop(lihr:end)+1e-10);

    % get population-weighted average ifr
    weighted_ifr       = Npop.*dis.ifr;
    % use weight to compute life expectancy lost per death
    life_years_lost_per_death         = arrayfun(@(x) dot(life_expectancy(x{1}),weighted_ifr(x{1}))/sum(weighted_ifr(x{1})), ageindex); 
    data.yll   = [repmat(life_years_lost_per_death(data.adInd),1,45),life_years_lost_per_death];

end




