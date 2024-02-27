function betamod = betamod_wrapped(ddk, p2, data, mandate)


    if mandate==1
        betamod = ones(size(ddk));
    else
        betamod = social_distancing(p2.sdl,p2.sdb,ddk,data.rel_mobility(mandate));
        if any(mandate==data.imand)
            betamod = min(betamod, social_distancing(p2.sdl,p2.sdb,2,data.rel_mobility(mandate)));
        end
    end

end