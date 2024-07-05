% wrapper to compute scalar that multiplies beta, the transmission rate, as
% in some circumstances it should be 1, and to impose a lower bound
%
% ddk: deaths per 10^6 population at this time point
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate
% rel_stringency: fraction, approximation of impact of mandate on R
% t: current time
%
% betamod: scalar between 0 and 1

function sd = betamod_wrapped(ddk, data, mandate, rel_stringency, t)
    
    if mandate==1   % means no mandate
        sd = ones(size(ddk));
    else
        baseline = 1 - (1-data.sd_baseline)./(1+data.sd_decay_rate).^t;
        death_coef = data.sd_death_coef;
        mandate_coef = data.sd_mandate_coef;
        
        sd_with_mandate = social_distancing(baseline, death_coef, mandate_coef,ddk, rel_stringency);
        if any(mandate==data.imand)
            sd_with_mandate = min(sd_with_mandate, social_distancing(baseline, death_coef, mandate_coef, 20, rel_stringency));
        end
        sd = min(1, sd_with_mandate ./ (1-rel_stringency).^.5) ;
    end
end



