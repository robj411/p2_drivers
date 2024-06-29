% computes force of infection (foi). Accounts for social distancing and
% corrects for sd already explained by mandate.
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
% contact_matrix: contact matrix
%
% foi: force of infection, vector

function foi = get_foi(dis, hospital_occupancy, data, mandate,...
        Ina,Ins,Inav1,Insv1,Inav2,Insv2,contact_matrix)
    
    NN0 = data.NNs;
    phi = 1 .* dis.rr_infection;  %+data.amp*cos((t-32-data.phi)/(365/2*pi));
    
    red = dis.red; % relative reduction in infectiousness of asymptomatic
    trv1 = dis.trv1; % relative reduction in infectiousness of BPSV-vaccinated
    trv2 = dis.trv2; % relative reduction in infectiousness of SARS-X--vaccinated
    beta = dis.beta;
    
    %% number of infectious people
    
    I       = red*Ina+Ins +(1-trv1).*(red*Inav1+Insv1) + (1-trv2).*(red*Inav2+Insv2) ;   
    
    %% seed
    seed = 10^(-data.seedsize);
    if (data.imand==2 && ... % elimination strategy
        mandate < 5 && mandate > 1 || ... % using closures
        data.imand==3 && ... % reactive closure strategy
        mandate == 3) || ... % currently in closures
        sum(I) == 0 % importation has not yet occurred
        seed = 0;
    end
    
    %% correct for mandate
    
    contact_matrix_open = data.Dvec(:,:,1);
    Ifrac = I./NN0;
    foi0 = contact_matrix_open*Ifrac;
    foi1 = contact_matrix*Ifrac;
    sd_so_far = ((foi1'*NN0+1e-10)./(foi0'*NN0+1e-10));
    
    %% social distancing
    weighted_deaths = sum(dis.mu.*hospital_occupancy); 
    sd = betamod_wrapped(10^6*weighted_deaths/sum(NN0), ...
        data, mandate, 1-sd_so_far);
%     new_betamod = 1;%sd./sd_so_far;
    
    %% foi
    foi     = phi.*beta.*(sd.*contact_matrix)*Ifrac + seed;
    
end



