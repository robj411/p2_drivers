function R = get_R(ntot, dis, S,Sv1,Sv2, N, contact_matrix, beta, betamod, p3, p4)

    h = dis.h;
    g2 = dis.g2;
    
    S_sum = S + Sv1 + Sv2;
    S_frac = S ./ S_sum;
    V_frac = Sv1 ./ S_sum;
    B_frac = Sv2 ./ S_sum;
    
    FOI  = contact_matrix .* beta .* betamod .* repmat(S_sum.*dis.rr_infection,1,ntot)./repmat(N',ntot,1);
    onesn = ones(ntot,1);
    
    sig1 = dis.sig1;
    sig2 = dis.sig2;
    
%     if sum(Sv1+Sv2)==0
%             
%         F     = zeros(3*ntot,3*ntot);
%         F(1:ntot, ntot+1:end) = [dis.red*FOI, FOI];
% 
%         vvec = [(sig1+sig2).*onesn;      (dis.g1).*onesn;       (g2 + h).*onesn];%g2 and h are vectors
% 
%         V    = diag(vvec);
%         V(ntot+1:2*ntot,1:ntot)   = diag(-sig1.*onesn);
%         V(2*ntot+1:3*ntot,1:ntot) = diag(-sig2.*onesn);
% 
%     else
    % code should be the same but R0 is computed before trv1 etc are defined 
        
        red = (1-p3).*dis.red;
        trv1 = dis.trv1;
        trv2 = dis.trv2;
        
        S_frac_mat = repmat(S_frac,1,ntot);
        V_frac_mat = repmat(V_frac,1,ntot);
        B_frac_mat = repmat(B_frac,1,ntot);

        FOIu = S_frac_mat.*FOI;
        FOIv = V_frac_mat.*FOI.*(1-dis.scv1);
        FOIb = B_frac_mat.*FOI.*(1-dis.scv2);
        
        FOIweighted = FOIu + FOIv + FOIb;
        
        FOIin = (S_frac_mat + V_frac_mat.*(1-trv1) + B_frac_mat.*(1-trv2)).*FOIweighted;
        
        
        F     = zeros(3*ntot,3*ntot);
        F(1:ntot, ntot+1:end) = [red * FOIin, (1-p4).*FOIin];
                                           
        g2hweighted = S_frac.*(g2+h) + V_frac.*(dis.g2_v1+dis.h_v1) + B_frac.*(dis.g2_v2+dis.h_v2);
        
        vvec = [(sig1+sig2).*onesn;      (dis.g1).*onesn;       g2hweighted.*onesn];

        V    = diag(vvec);
        V(ntot+1:2*ntot,1:ntot)   = diag(-sig1.*onesn);
        V(2*ntot+1:3*ntot,1:ntot) = diag(-sig2.*onesn);

%     end

    NGM = F/V;
    d = eigs(NGM,1);%largest in magnitude (+/-) 
    R = max(d); 
    
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
