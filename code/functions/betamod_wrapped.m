% wrapper to compute scalar that multiplies beta, the transmission rate, as
% in some circumstances it should be 1, and to impose a lower bound

% ddk: deaths per 10^5 population at this time point
% p2: struct of p2 intervention parameters
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate

% betamod: scalar between 0 and 1

function betamod = betamod_wrapped(ddk, p2, data, mandate)

    if mandate==1 % means no mandate
        betamod = ones(size(ddk));
    else
        betamod = social_distancing(p2.sdl,p2.sdb,ddk,data.rel_mobility(mandate));
        if any(mandate==data.imand)
            betamod = min(betamod, social_distancing(p2.sdl,p2.sdb,2,data.rel_mobility(mandate)));
        end
    end

end