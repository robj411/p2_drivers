% function to compute the total reduction in transmission,
% which is a function of current closures, current mandate and current deaths
%
% baseline: minimum value for reduction, to which function tends for high values
% of deaths and mandate
% death_coef: coefficient for number of deaths per million
% mandate_coef: coefficient for current mandate
% deaths_per_mil: current deaths per million
% rel_stringency: current mandate
%
% reduction: amount by which to scale transmission. Value between 0 and 1.


function reduction = transmission_reduction(baseline, death_coef, mandate_coef,...
                                    deaths_per_mil, rel_stringency) 
    
    % same function as in R code where parameters fitted
    reduction = 1./(1+death_coef.*deaths_per_mil + mandate_coef.*rel_stringency) .* (1-baseline) + baseline;
    
    % if utr computed is higher than mandate, we should return 1
    reduction = min(1, max(reduction, baseline));     

end



