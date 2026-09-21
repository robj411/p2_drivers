% compute the fraction of people identified by testing and then the
% fraction of infectiousness averted as a consequence

% p2: struct of p2 intervention parameters
% Ip: number of cases per 10^5 population

% p3: fraction of infectiousness averted

function frac_cases_found = get_case_ID_rate(test_rate, Ip)
    
    a0    = 2.197;
    a1    = 0.1838;
    a2    = -1.024;
    
    frac_cases_found = 1./(1+exp(a0+a1*Ip+a2*log10(test_rate)));
    
    frac_cases_found(Ip >= test_rate) = min(frac_cases_found(Ip >= test_rate),test_rate(Ip >= test_rate)./10^5);
    
    frac_cases_found = max(frac_cases_found, test_rate/10^5 );
    
    
    
end
    
            
            