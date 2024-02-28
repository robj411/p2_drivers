% function to compute the impact of social distancing on transmission,
% which is a function of current closures and current deaths

% lower_bound
% rate
% deaths_per_10k
% current_drop

% sd_out: amount by which to scale transmission. Value between 0 and 1.


% function sd_out = social_distancing(lower_bound,rate,deaths_per_10k,current_drop ) 
function sd_out = social_distancing(baseline, death_coef, mandate_coef,deaths_per_10k,current_drop ) 

    rel_mobility = 1-current_drop ;
    
    
%     sd = (lower_bound-rate)+(1-lower_bound+rate)*(1+((lower_bound-1)/(1-lower_bound+rate))).^(deaths_per_10k./10);
    sd = 1./(1+death_coef.*deaths_per_10k + mandate_coef.*current_drop) .* (1-baseline) + baseline;
    sd = min(rel_mobility, max(sd, baseline)); % if sd computed is higher than mandate, we should return 1
    
    sd_out = sd./rel_mobility;
    
%     max_final_reduction = 1 - baseline;
%     maxdrop = min(maxdrop, max_final_reduction-current_drop);
%     
%     final_reduction = current_drop + (1-sd).*maxdrop;
%     final_mobility = 1-final_reduction;
%     sd_out = final_mobility./rel_mobility;

end
%        sd_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
