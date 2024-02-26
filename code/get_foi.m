function foi = get_foi(dis, hospital_occupancy, p2, data, mandate,...
        Ina,Ins,Inav1,Insv1,Inav2,Insv2,D,NN0)

    phi = 1 .* dis.rr_infection;  %+data.amp*cos((t-32-data.phi)/(365/2*pi));

    betamod = betamod_wrapped(10^5*sum(dis.mu.*hospital_occupancy)/sum(NN0), ...
        p2, data, mandate);
    
    red = dis.red;
    trv1 = dis.trv1;
    trv2 = dis.trv2;
    beta = dis.beta;
    
    
    I       = red*Ina+Ins +(1-trv1).*(red*Inav1+Insv1) + (1-trv2).*(red*Inav2+Insv2) ;   
    foi     = phi.*beta.*betamod.*(D*(I./NN0));
    
end