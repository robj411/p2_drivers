% get effective reproduction number
%
% nStrata: number of strata
% dis: struct of pathogen parameters
% S: susceptible unvaccinated
% Sv1: susceptible BPSV-vaccinated
% Sv2: susceptible SARS-X--vaccinated
% p3: fraction of asymptomatic infectious people's infectiousness averted
% p4: fraction of symptomatic infectious people's infectiousness averted
% N: population by stratum
% contact_matrix: contact matrix
%
% R: effective reproduction number

function CI = get_candidate_infectees(nStrata, dis, S,Sv1,Sv2, p3, p4, N, contact_matrix)
    
    h = dis.h;
    g2 = dis.g2;
    
    S_sum = S + Sv1 + Sv2;
    S_frac = S ./ S_sum;
    V_frac = Sv1 ./ S_sum;
    B_frac = Sv2 ./ S_sum;
    
    FOI  = contact_matrix .* repmat(S_sum.*dis.rr_infection,1,nStrata)./repmat(N',nStrata,1);
    onesn = ones(nStrata,1);
    
    sig1 = dis.sig1;
    sig2 = dis.sig2;
        
    red = (1-p3).*dis.red;
    trv1 = dis.trv1;
    trv2 = dis.trv2;

    S_frac_mat = repmat(S_frac,1,nStrata);
    V_frac_mat = repmat(V_frac,1,nStrata);
    B_frac_mat = repmat(B_frac,1,nStrata);

    FOIu = S_frac_mat.*FOI;
    FOIv = V_frac_mat.*FOI.*(1-dis.scv1);
    FOIb = B_frac_mat.*FOI.*(1-dis.scv2);

    FOIweighted = FOIu + FOIv + FOIb;

    FOIin = (S_frac_mat + V_frac_mat.*(1-trv1) + B_frac_mat.*(1-trv2)).*FOIweighted;


    F     = zeros(3*nStrata,3*nStrata);
    F(1:nStrata, nStrata+1:end) = [red * FOIin, (1-p4).*FOIin];

    g2hweighted = S_frac.*(g2+h) + V_frac.*(dis.g2_v1+dis.h_v1) + B_frac.*(dis.g2_v2+dis.h_v2);

    vvec = [(sig1+sig2).*onesn;      (dis.g1).*onesn;       g2hweighted.*onesn];

    n = length(vvec);
    V    = zeros(n); 
    V(1:n+1:end) = vvec;
    nmat = eye(nStrata);
    V(nStrata+1:2*nStrata,1:nStrata)   = -sig1.*nmat;
    V(2*nStrata+1:3*nStrata,1:nStrata) = -sig2.*nmat;


    NGM = F/V;
    d = eigs(NGM,1);%largest in magnitude (+/-) 
    CI = max(d); 
    
end

function Rt = get_R_PD(ntot, dis, h, g2, S,Sv2, N, contact_matrix, beta, betamod, p3, p4)


    FOIu = repmat(S.*dis.rr_infection,1,ntot).*beta.*betamod.*contact_matrix./repmat(N',ntot,1);
    FOIv = repmat(Sv2.*dis.rr_infection,   1,ntot).*beta.*betamod.*contact_matrix./repmat(N',ntot,1).*(1-dis.scv2);

    F                                = zeros(6*ntot,6*ntot);
    F(1:ntot,       2*ntot+1:6*ntot) = [dis.red*FOIu,  FOIu,  dis.red*(1-dis.trv2)*FOIu,  (1-dis.trv2)*FOIu];
    F(ntot+1:2*ntot,2*ntot+1:6*ntot) = [dis.red*FOIv,  FOIv,  dis.red*(1-dis.trv2)*FOIv,  (1-dis.trv2)*FOIv];

    onesn                            = ones(ntot,1);
    vvec                             = [(dis.sig1+dis.sig2).*onesn;  (dis.sig1+dis.sig2).*onesn;  (dis.g1+p3).*onesn;...
                                        (g2+h+p4).*onesn;            (dis.g1+p3).*onesn;          (dis.g2_v2+dis.h_v2+p4).*onesn];
    V                                = diag(vvec);
    V(2*ntot+1:3*ntot,1:ntot)        = diag(-dis.sig1.*onesn);
    V(3*ntot+1:4*ntot,1:ntot)        = diag(-dis.sig2.*onesn);
    V(4*ntot+1:5*ntot,ntot+1:2*ntot) = diag(-dis.sig1.*onesn);
    V(5*ntot+1:6*ntot,ntot+1:2*ntot) = diag(-dis.sig2.*onesn);

    NGM = F/V;
    Rt  = eigs(NGM,1);

end
