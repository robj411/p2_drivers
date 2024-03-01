% function to compute the impact of social distancing on transmission,
% which is a function of current closures and current deaths

% lower_bound
% rate
% deaths_per_10k
% current_drop

% sd_out: amount by which to scale transmission. Value between 0 and 1.


% function sd_out = social_distancing(lower_bound,rate,deaths_per_10k,current_drop ) 
function sd_out = social_distancing(baseline, death_coef, mandate_coef,...
                                    deaths_per_mil,rel_mobility, rel_stringency) 
    
    % same function as in R code where parameters fitted
    sd = 1./(1+death_coef.*deaths_per_mil + mandate_coef.*rel_stringency) .* (1-baseline) + baseline;
    
    % if sd computed is higher than mandate, we should return 1
    sd = min(rel_mobility, max(sd, baseline)); 
    
    % divide through as we already account for the mandate impact on
    % transmission
    sd_out = sd./rel_mobility;
    

end
%        sd_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
