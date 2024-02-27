% computes force of infection (foi)

% dis: struct of pathogen parameters
% hospital_occupancy: number of people in hospital 
% p2: struct of p2 intervention parameters
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate       
% Ina: number of unvaccinated infectious not self isolating asymptomatic people
% Ins: number of unvaccinated infectious not self isolating symptomatic people
% Inav1: number of BPSV-vaccinated infectious not self isolating asymptomatic people
% Insv1: number of BPSV-vaccinated infectious not self isolating symptomatic people
% Inav2: number of SARS-X--vaccinated infectious not self isolating asymptomatic people
% Insv2: number of SARS-X--vaccinated infectious not self isolating symptomatic people
% D: contact matrix
% NN0: number of people per stratum

% foi: force of infection, vector

function foi = get_foi(dis, hospital_occupancy, p2, data, mandate,...
        Ina,Ins,Inav1,Insv1,Inav2,Insv2,D,NN0)
    
    phi = 1 .* dis.rr_infection;  %+data.amp*cos((t-32-data.phi)/(365/2*pi));

    betamod = betamod_wrapped(10^5*sum(dis.mu.*hospital_occupancy)/sum(NN0), ...
        p2, data, mandate);
    
    red = dis.red; % relative reduction in infectiousness of asymptomatic
    trv1 = dis.trv1; % relative reduction in infectiousness of BPSV-vaccinated
    trv2 = dis.trv2; % relative reduction in infectiousness of SARS-X--vaccinated
    beta = dis.beta;
    
    I       = red*Ina+Ins +(1-trv1).*(red*Inav1+Insv1) + (1-trv2).*(red*Inav2+Insv2) ;   
    foi     = phi.*beta.*betamod.*(D*(I./NN0));
    
end