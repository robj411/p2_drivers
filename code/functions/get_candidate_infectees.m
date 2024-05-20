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
    
    sig1 = dis.sig1;
    sig2 = dis.sig2;
        
    red = (1-p3).*dis.red;
    trv1 = dis.trv1;
    trv2 = dis.trv2;
    
    S_sum = S + Sv1 + Sv2;
    S_frac = S ./ S_sum;
    V_frac = Sv1 ./ S_sum;
    B_frac = Sv2 ./ S_sum;
    
    FOI  = contact_matrix .* repmat(S_sum.*dis.rr_infection,1,nStrata)./repmat(N',nStrata,1);
    onesn = ones(nStrata,1);

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
    
    %% collapse to 4
%     h4 = h(46:49);
%     g24 = g2(46:49);
%     
%     adinds = [1:45,48];
%     collapse_ad = @(x) [x(45+1); x(45+2); sum(x(adinds)); x(end)];
%     
%     S4 = collapse_ad(S);
%     Sv14 = collapse_ad(Sv1);
%     Sv24 = collapse_ad(Sv2);
%     N4 = collapse_ad(N);
%     nStrata4 = 4;
%     
%     S_sum4 = S4 + Sv14 + Sv24;
%     S_frac4 = S4 ./ S_sum4;
%     V_frac4 = Sv14 ./ S_sum4;
%     B_frac4 = Sv24 ./ S_sum4;
%     
%     contact_matrix494 = [contact_matrix(:,45+(1:2)) sum(contact_matrix(:,adinds),2) contact_matrix(:,end)];
%     contact_matrix4 = zeros(nStrata4);
%     for i=1:nStrata4
%         contact_matrix494(adinds,i) = contact_matrix494(adinds,i).*N(adinds)/sum(N(adinds));
%         contact_matrix4(:,i) = collapse_ad(contact_matrix494(:,i));
%     end
%     
% 
%     S_frac_mat4 = repmat(S_frac4,1,nStrata4);
%     V_frac_mat4 = repmat(V_frac4,1,nStrata4);
%     B_frac_mat4 = repmat(B_frac4,1,nStrata4);
%     
%     FOI4  = contact_matrix4 .* repmat(S_sum4.*dis.rr_infection,1,nStrata4)./repmat(N4',nStrata4,1);
%     onesn4 = ones(nStrata4,1);
%     
%     scv1 = dis.scv1;
%     scv2 = dis.scv2;
%     if length(scv1)>4
%         scv1 = scv1(46:49);
%         scv2 = scv2(46:49);
%         trv1 = trv1(46:49);
%         trv2 = trv2(46:49);
%     end
%     FOIu4 = S_frac_mat4.*FOI4;
%     FOIv4 = V_frac_mat4.*FOI4.*(1-scv1);
%     FOIb4 = B_frac_mat4.*FOI4.*(1-scv2);
%     
%     FOIweighted4 = FOIu4 + FOIv4 + FOIb4;
% 
%     FOIin4 = (S_frac_mat4 + V_frac_mat4.*(1-trv1) + B_frac_mat4.*(1-trv2)).*FOIweighted4;
% 
%     
%     F4     = zeros(3*nStrata4,3*nStrata4);
%     F4(1:nStrata4, nStrata4+1:end) = [red * FOIin4, (1-p4).*FOIin4];
% 
%     g2hweighted4 = S_frac4.*(g24+h4) + V_frac4.*(dis.g2_v1(46:49)+dis.h_v1(46:49)) + B_frac4.*(dis.g2_v2(46:49)+dis.h_v2(46:49));
% 
%     vvec4 = [(sig1+sig2).*onesn4;      (dis.g1).*onesn4;       g2hweighted4.*onesn4];
% 
%     n4 = length(vvec4);
%     V4    = zeros(n4); 
%     V4(1:n4+1:end) = vvec4;
%     nmat4 = eye(nStrata4);
%     V4(nStrata4+1:2*nStrata4,1:nStrata4)   = -sig1.*nmat4;
%     V4(2*nStrata4+1:3*nStrata4,1:nStrata4) = -sig2.*nmat4;
% 
% 
%     NGM4 = F4/V4;
%     d4 = eigs(NGM4,1);%largest in magnitude (+/-) 
%     CI = max(d4); 
    
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
