function R = get_R(ntot, dis, h, g2, S,Sv1,Sv2, N, contact_matrix, beta, betamod, p3, p4)

    FOI  = contact_matrix .* beta .* betamod .* repmat(S.*dis.rr_infection,1,ntot)./repmat(N',ntot,1);
    
    onesn = ones(ntot,1);
    
    if sum(Sv1+Sv2)==0
    
        F     = zeros(3*ntot,3*ntot);
        F(1:ntot, ntot+1:end) = [dis.red*FOI, FOI];

        vvec = [(dis.sig1+dis.sig2).*onesn;      (dis.g1 + p3).*onesn;       (g2 + h + p4).*onesn];%g2 and h are vectors

        V    = diag(vvec);
        V(ntot+1:2*ntot,1:ntot)   = diag(-dis.sig1.*onesn);
        V(2*ntot+1:3*ntot,1:ntot) = diag(-dis.sig2.*onesn);

    else
        
        sig1 = dis.sig1;
        sig2 = dis.sig2;
        red = dis.red;
        trv1 = dis.trv1;
        trv2 = dis.trv2;
    %     beta = dis.beta;
        g1 = dis.g1;
        red1 = red*(1-trv1);
        red2 = red*(1-trv2);

    %     FOI = beta.*betamod.*contact_matrix./repmat(N',ntot,1);

        FOIu = repmat(S,1,ntot).*FOI;
        FOIv = repmat(Sv1,   1,ntot).*FOI.*(1-dis.scv1);
        FOIb = repmat(Sv2,   1,ntot).*FOI.*(1-dis.scv2);

        F                                = zeros(9*ntot,9*ntot);
        F(1:ntot,          3*ntot+1:9*ntot) = [red*FOIu, FOIu, red1*FOIu, (1-trv1)*FOIu, red2*FOIu,  (1-trv2)*FOIu];
        F(  ntot+1:2*ntot, 3*ntot+1:9*ntot) = [red*FOIv, FOIv, red1*FOIv, (1-trv1)*FOIv, red2*FOIv,  (1-trv2)*FOIv];
        F(2*ntot+1:3*ntot, 3*ntot+1:9*ntot) = [red*FOIb, FOIb, red1*FOIb, (1-trv1)*FOIb, red2*FOIb,  (1-trv2)*FOIb];

        vvec                             = [(sig1+sig2).*onesn;  (sig1+sig2).*onesn;  (sig1+sig2).*onesn;  ...
                                            (g1+p3).*onesn;      (g2+h+p4).*onesn;    ...
                                            (g1+p3).*onesn;      (dis.g2_v1+dis.h_v1+p4).*onesn;    ...
                                            (g1+p3).*onesn;      (dis.g2_v2+dis.h_v2+p4).*onesn];

        V                                = diag(vvec);
        V(2*ntot+1:3*ntot,1:ntot)        = diag(-sig1.*onesn);
        V(3*ntot+1:4*ntot,1:ntot)        = diag(-sig2.*onesn);
        V(4*ntot+1:5*ntot,ntot+1:2*ntot) = diag(-sig1.*onesn);
        V(5*ntot+1:6*ntot,ntot+1:2*ntot) = diag(-sig2.*onesn);
        V(6*ntot+1:7*ntot,2*ntot+1:3*ntot) = diag(-sig1.*onesn);
        V(7*ntot+1:8*ntot,2*ntot+1:3*ntot) = diag(-sig2.*onesn);
    end

    NGM = F/V;
%     Rt  = eigs(NGM,1);
    d = eigs(NGM,1);%largest in magnitude (+/-) 
    R = max(d); 
    
end
