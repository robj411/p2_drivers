% computes fractions of infectious people self isolating. Depends on
% mandate, testing, and symptomaticity.
%
% sumI: total number of infectious people
% sumN: total number of people
% p2: struct of p2 intervention parameters
% t: current time (day)
% mandate: integer corresponding to states of government mandate
%
% p3: fraction of asymptomatic infectious people's infectiousness averted
% p4: fraction of symptomatic infectious people's infectiousness averted


function [p3, p4] = fraction_averted_self_isolating(sumI, sumN, p2, t, mandate)


    test_start_time = p2.t_tit;
    test_end_time = max(p2.tpoints);
    valid_timepoints = (t>=test_start_time);%.*(t<=test_end_time);
        
    if any(valid_timepoints==1) & mandate~=5
        Ip = 10^5*sumI/sumN;
        test_rate = p2.trate .* max(min((t-test_start_time)./30, 1),0);
        frac_cases_found = get_case_ID_rate(test_rate, Ip);
        
        p3 = frac_cases_found * p2.frac_asym_infectiousness_averted ;
        p3 = p3.*(valid_timepoints); 
        
        presym_averted = frac_cases_found * p2.frac_presym_infectiousness_averted ;
        presym_averted = presym_averted.*(valid_timepoints); 
        sym_averted = p2.frac_sym_infectiousness_averted + presym_averted;
        p4 = sym_averted;
    else
        p3 = zeros(size(t));
        p4 = p3;
        p4(t>=test_start_time) = p2.frac_sym_infectiousness_averted;
    end
    
%     p4 = p3;
end




