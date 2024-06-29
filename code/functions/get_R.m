% get effective reproduction number
%
% nStrata: number of strata
% dis: struct of pathogen parameters
% S: susceptible unvaccinated
% Sv1: susceptible BPSV-vaccinated
% Sv2: susceptible SARS-X--vaccinated
% beta: rate of infection
% p3: fraction of asymptomatic infectious people's infectiousness averted
% p4: fraction of symptomatic infectious people's infectiousness averted
% ddk: deaths per 10^6 population at this time point
% data: struct of general model parameters
% mandate: integer corresponding to states of government mandate  
%
% R: effective reproduction number

function R = get_R(nStrata, dis, S,Sv1,Sv2, beta, p3, p4, ddk, data, mandate)

    N = data.NNs;
    
    contact_matrix = data.Dvec(:,:,mandate);
    
    foi0 = data.basic_foi;
    foi1 = contact_matrix*(1./N);
    sd_so_far = ((foi1'*N+1e-10)./(foi0'*N+1e-10));

    sd = betamod_wrapped(ddk, data, mandate, 1-sd_so_far);
%     betamod = 1;%sd./sd_so_far;
    
    CI = get_candidate_infectees(nStrata, dis, S,Sv1,Sv2, p3, p4, N, sd.*contact_matrix);
    R  = beta .* CI;
    

end



