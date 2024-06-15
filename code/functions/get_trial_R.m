% get effective reproduction number for scaled configuration
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
% openness: vector, economic configuration expressed as a fraction
% home_working: vector, fraction of workforce working from home
%
% R: effective reproduction number

function R = get_trial_R(nStrata, dis, S,Sv1,Sv2, beta, p3, p4, ddk, data, openness, home_working, mandate)

    N = data.NNs;
    
    contact_matrix = p2MakeDs(data,N,openness,home_working);
    
    % current_CI = get_candidate_infectees(nStrata, dis, N,Sv1,Sv2, p3, p4, N, contact_matrix);
    % ref_CI = get_candidate_infectees(nStrata, dis, N,Sv1,Sv2, p3, p4, N, data.Dvec(:,:,1));
    % 
    % foi0 = data.basic_foi;
    % foi1 = contact_matrix*(1./N);
    % sd_so_far = current_CI/ref_CI;%(sum(foi1+1e-10)./sum(foi0+1e-10));

    foi0 = data.basic_foi;
    foi1 = contact_matrix*(1./N);
    sd_so_far = ((foi1'*N+1e-10)./(foi0'*N+1e-10));


    sd = betamod_wrapped(ddk, data, mandate, 1-sd_so_far);
    betamod = sd./sd_so_far;
    % [sd sd_so_far betamod]
    
    CI = get_candidate_infectees(nStrata, dis, S,Sv1,Sv2, p3, p4, N, betamod.* contact_matrix);
    R  = beta .* CI;

end



