% function to compute the impact of social distancing on transmission,
% which is a function of current closures and current deaths
%
% baseline: minimum value for sd, to which function tends for high values
% of deaths and mandate
% death_coef: coefficient for number of deaths per million
% mandate_coef: coefficient for current mandate
% deaths_per_mil: current deaths per million
% rel_stringency: current mandate
%
% sd: amount by which to scale transmission. Value between 0 and 1.


% function sd_out = social_distancing(lower_bound,rate,deaths_per_10k,current_drop ) 
function sd = social_distancing(baseline, death_coef, mandate_coef,...
                                    deaths_per_mil, rel_stringency) 
    
    % same function as in R code where parameters fitted
    sd = 1./(1+death_coef.*deaths_per_mil + mandate_coef.*rel_stringency) .* (1-baseline) + baseline;
    
    % if sd computed is higher than mandate, we should return 1
    sd = min(1, max(sd, baseline));     

end



