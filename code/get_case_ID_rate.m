function p3 = get_case_ID_rate(p2, Ip)

    trate = p2.trate;
    
    b0    = 2.197;
    b1    = 0.1838;
    b2    = -1.024;
    
    coef = 1./(1+exp(b0+b1*Ip+b2*log10(trate)));
    
    coef(Ip >= trate) = min(coef(Ip >= trate),trate/10^5);
    
    coef = max(coef, trate/10^5 .* Ip/10^5);
    
    p3 = p2.self_isolation_compliance .* coef / p2.dur;
    
end
    
            
            