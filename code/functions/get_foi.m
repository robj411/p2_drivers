% computes force of infection (foi)
%
% dis: struct of pathogen parameters
% hospital_occupancy: number of people in hospital 
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate       
% Ina: number of unvaccinated infectious not self isolating asymptomatic people
% Ins: number of unvaccinated infectious not self isolating symptomatic people
% Inav1: number of BPSV-vaccinated infectious not self isolating asymptomatic people
% Insv1: number of BPSV-vaccinated infectious not self isolating symptomatic people
% Inav2: number of SARS-X--vaccinated infectious not self isolating asymptomatic people
% Insv2: number of SARS-X--vaccinated infectious not self isolating symptomatic people
% D: contact matrix
%
% foi: force of infection, vector

function foi = get_foi(dis, hospital_occupancy, data, mandate,...
        Ina,Ins,Inav1,Insv1,Inav2,Insv2,D)
    
    NN0 = data.NNs;
    phi = 1 .* dis.rr_infection;  %+data.amp*cos((t-32-data.phi)/(365/2*pi));

    sd = betamod_wrapped(10^6*sum(dis.mu.*hospital_occupancy)/sum(NN0), ...
        data, mandate);
    
    red = dis.red; % relative reduction in infectiousness of asymptomatic
    trv1 = dis.trv1; % relative reduction in infectiousness of BPSV-vaccinated
    trv2 = dis.trv2; % relative reduction in infectiousness of SARS-X--vaccinated
    beta = dis.beta;
    
    I       = red*Ina+Ins +(1-trv1).*(red*Inav1+Insv1) + (1-trv2).*(red*Inav2+Insv2) ;   
%     foi     = phi.*beta.*betamod.*(D*(I./NN0));

    D0 = data.Dvec(:,:,1);
    Ifrac = I./NN0;
    foi0 = D0*Ifrac;
    foi1 = D*Ifrac;
    sd_so_far = median((foi1+1e-10)./(foi0+1e-10));
%     sd = betamod*data.rel_mobility(mandate);
    new_betamod = sd./sd_so_far;
    foi     = phi.*beta.*(new_betamod.*D)*Ifrac;
    

    
end