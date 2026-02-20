% wrapper to compute scalar that multiplies beta, the transmission rate, as
% in some circumstances it should be 1, and to impose a lower bound
%
% deaths_per_mill: deaths per 10^6 population at this time point
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate
% rel_stringency: fraction, approximation of impact of mandate on R
% t: current time
%
% betamod: scalar between 0 and 1

function utr = betamod_wrapped(deaths_per_mill, data, mandate, rel_stringency, t)
    
    if mandate==1   % means no mandate
        utr = ones(size(deaths_per_mill));
    else
        % baseline = data.utr_baseline; %1 - (1-data.utr_baseline).*(data.utr_decay_rate);
        baseline = 1 - (1-data.utr_baseline) ./ (1 + data.utr_decay_rate).^t;
        death_coef = data.utr_death_coef;
        mandate_coef = data.utr_mandate_coef;
        
        reduction_with_mandate = transmission_reduction(baseline, death_coef, mandate_coef,deaths_per_mill, rel_stringency);
        if any(mandate==data.imand)
            reduction_with_mandate = min(reduction_with_mandate, transmission_reduction(baseline, death_coef, mandate_coef, 20, rel_stringency));
        end
        utr = min(1, reduction_with_mandate ./ (1-rel_stringency)) ;
    end
%     utr = ones(size(deaths_per_mill));
end



