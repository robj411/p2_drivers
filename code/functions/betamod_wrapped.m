% wrapper to compute scalar that multiplies beta, the transmission rate, as
% in some circumstances it should be 1, and to impose a lower bound
%
% ddk: deaths per 10^6 population at this time point
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate
% rel_stringency: fraction, approximation of impact of mandate on R
%
% betamod: scalar between 0 and 1

function rel_betamod = betamod_wrapped(ddk, data, mandate, rel_stringency)
    
    if mandate==1   % means no mandate
        rel_betamod = ones(size(ddk));
    else
        baseline = data.sd_baseline;
        death_coef = data.sd_death_coef;
        mandate_coef = data.sd_mandate_coef;
        
        betamod = min(1-rel_stringency,social_distancing(baseline, death_coef, mandate_coef,ddk, rel_stringency));
%         betamod = social_distancing(baseline, death_coef, mandate_coef,ddk, rel_stringency);
        if any(mandate==data.imand)
            betamod = min(betamod, social_distancing(baseline, death_coef, mandate_coef, 20, rel_stringency));
        end
        rel_betamod = min(1, betamod ./ (1-rel_stringency)) ;
    end
%     rel_betamod = ones(size(ddk));
end



