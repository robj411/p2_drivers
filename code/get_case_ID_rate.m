function p3 = get_case_ID_rate(p2, Ip)

    trate = p2.trate;
    
    b0    = 2.197;
    b1    = 0.1838;
    b2    = -1.024;
    
    frac_cases_found = 1./(1+exp(b0+b1*Ip+b2*log10(trate)));
    
    frac_cases_found(Ip >= trate) = min(frac_cases_found(Ip >= trate),trate/10^5);
    
    frac_cases_found = max(frac_cases_found, trate/10^5 );
    
    frac_infectiousness_averted = 1 - p2.frac_asym_infectiousness_remaining;
    p3 = p2.self_isolation_compliance .* frac_cases_found * frac_infectiousness_averted;
    
end
    
            
            