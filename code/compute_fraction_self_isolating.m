function [p3, p4] = compute_fraction_self_isolating(sumI, sumN, p2, t, mandate)


    test_start_time = p2.t_tit;
    test_end_time = max(p2.tpoints);
    valid_timepoints = (t>=test_start_time).*(t<=test_end_time);
        
    if any(valid_timepoints==1) & mandate~=5
        Ip = 10^5*sumI/sumN;
        p3 = get_case_ID_rate(p2, Ip);
        p3 = p3.*(valid_timepoints); 
        p4 = 1 - p2.frac_sym_infectiousness_remaining;
    else
        p3 = zeros(size(t));
        p4 = p3;
    end
    
    p4 = p3;
end