% compute the fraction of people identified by testing and then the
% fraction of infectiousness averted as a consequence

% p2: struct of p2 intervention parameters
% Ip: number of cases per 10^5 population

% p3: fraction of infectiousness averted

function frac_cases_found = get_case_ID_rate(p2, Ip)

    trate = p2.trate;
    
    b0    = 2.197;
    b1    = 0.1838;
    b2    = -1.024;
    
    frac_cases_found = 1./(1+exp(b0+b1*Ip+b2*log10(trate)));
    
    frac_cases_found(Ip >= trate) = min(frac_cases_found(Ip >= trate),trate/10^5);
    
    frac_cases_found = max(frac_cases_found, trate/10^5 );
    
    
    
end
    
            
            