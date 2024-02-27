% function to compute the impact of social distancing on transmission,
% which is a function of current closures and current deaths

% lower_bound
% rate
% deaths_per_10k
% current_drop

% sd_out: amount by which to scale transmission. Value between 0 and 1.

function sd_out = social_distancing(lower_bound,rate,deaths_per_10k,current_drop ) 

    sd = (lower_bound-rate)+(1-lower_bound+rate)*(1+((lower_bound-1)/(1-lower_bound+rate))).^(deaths_per_10k./10);
    sd = min(1, max(sd, lower_bound)); %ones(size(deaths_per_10k)); %
    
    cutpoint1 = 1-0.55;
    cutpoint2 = 1-0.21;
    maxdrop1 = 0.38;
    maxdrop2 = 0.61;
    
    rel_mobility = 1-current_drop ;
    
    if current_drop < cutpoint1
        maxdrop = maxdrop1;
    elseif current_drop > cutpoint2
        maxdrop = maxdrop2;
    else
        maxdrop = (maxdrop2*(current_drop-cutpoint1) + maxdrop1*(cutpoint2-current_drop))./(cutpoint2-cutpoint1);
    end
    
    max_final_reduction = 0.91;
    maxdrop = min(maxdrop, max_final_reduction-current_drop);
    
    final_reduction = current_drop + (1-sd).*maxdrop;
    final_mobility = 1-final_reduction;
    sd_out = final_mobility./rel_mobility;

end
%        sd_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
