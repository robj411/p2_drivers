% get effective reproduction number

% nStrata: number of strata
% dis: struct of pathogen parameters
% S: susceptible unvaccinated
% Sv1: susceptible BPSV-vaccinated
% Sv2: susceptible SARS-X--vaccinated
% N: population by stratum
% contact_matrix: contact matrix
% beta: rate of infection
% betamod: scalar that multiplies beta
% p3: fraction of asymptomatic infectious people's infectiousness averted
% p4: fraction of symptomatic infectious people's infectiousness averted

% R: effective reproduction number

function R = get_R(nStrata, dis, S,Sv1,Sv2, beta, p3, p4, ddk, data, mandate)

    N = data.NNs;
    
    contact_matrix0 = data.Dvec(:,:,1);
    contact_matrix = data.Dvec(:,:,mandate);
    
    foi0 = contact_matrix0*(1./N);
    foi1 = contact_matrix*(1./N);
    sd_so_far = (foi1+1e-10)./(foi0+1e-10);

    sd = betamod_wrapped(ddk, data, mandate);
    betamod = sd./sd_so_far;
    
    CI = get_candidate_infectees(nStrata, dis, S,Sv1,Sv2, p3, p4, N, betamod.*contact_matrix);
    R  = beta .* CI;
    

end
