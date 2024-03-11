% compute current rates of transfer from S and R compartments to V
% compartments
%
% p2: struct of p2 intervention parameters
% t: current time, day 
% nStrata: number of strata
% R: numbers in recovered compartment
% S: numbers in susceptible compartment
% DE: numbers in died compartment
% Rv1: numbers in vaccinated recovered compartment
% Sv1: numbers in vaccinated susceptible compartment
%
% v1rates: rate of susptible unvaccinated to BPSV 
% v1rater: rate of recovered unvaccinated to BPSV  
% v2rates: rate of susptible unvaccinated to SARS-X  
% v2rater: rate of recovered unvaccinated to SARS-X  
% v12rates: rate of susptible BPSV to SARS-X   
% v12rater: rate of recovered BPSV to SARS-X  

%S (and R) ./nonVax accounts for the (inefficient) administration of vaccines to exposed, infectious and hospitalised people
%nonVax approximates S+E+I+H+R (unvaccinated) and D (partially)
%S (or R) ./nonVax is approximately 1 (or 0) when prevalence is low but is closer to 0 (or 1) when prevalence is high
%nonVax is non-zero as long as uptake is less than 100%
    
function [v1rates, v1rater, v2rates, v2rater, v12rates, v12rater] = ...
        get_vax_rates(p2, t, nStrata, R,S,DE, Rv1,Sv1)
    
    
    v1rates = zeros(nStrata,1);
    v1rater = zeros(nStrata,1);
    v2rates = zeros(nStrata,1);
    v2rater = zeros(nStrata,1);
    v12rates = zeros(nStrata,1);
    v12rater = zeros(nStrata,1);
        
    tpoints = p2.tpoints;
    if t >= tpoints(1) && t <= max(tpoints) 
        current_time = find(t>tpoints,1,'last');
        current_group = p2.group_order(current_time);
        arate = p2.arate;
        NNnext = p2.NNnext;
        targets = zeros(size(NNnext));
        if current_group == 3
            adInd = [1:(nStrata-4),(nStrata-4) + 3];
            targets(adInd) = 1;
        elseif current_group ~= 0 % current_group=0 in the gap between BPSV and specific/booster vaccine
            targets((nStrata-4) + current_group) = 1;
        end
        
        total_to_vax = targets.*NNnext.*arate;
        if t > p2.t_vax2
            if current_group == 4
                % populate v2rate and v12rate from S, Sv1, R, Rv1
                denom2 = R+S+DE+1e-15;
                denom12 = Rv1+Sv1+1e-15;
                vrate =  total_to_vax ./ (denom12 + denom2);
                v12rates = vrate.*Sv1;
                v12rater = vrate.*Rv1;
                v2rates = vrate.*S;
                v2rater = vrate.*R;
            else
                % populate v2rate from S, R
                denom2 = R+S+DE+1e-15;
                v2rate =  total_to_vax ./ denom2;
                v2rates = v2rate.*S;
                v2rater = v2rate.*R;
            end
            
        else
            % populate v1 from S, R
            denom = R+S+DE+1e-15;
            v1rate = total_to_vax ./ denom;
            v1rates = v1rate.*S;
            v1rater = v1rate.*R;
            
        end
    end
end

