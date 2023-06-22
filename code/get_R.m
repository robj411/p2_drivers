function R = get_R(ntot, dis, h, g2, S, N, contact_matrix, beta, betamod, p3, p4)

    FOI  = contact_matrix .* beta .* betamod .* repmat(S.*dis.rr_infection,1,ntot)./repmat(N',ntot,1);
    
    onesn = ones(ntot,1);
    
    F     = zeros(3*ntot,3*ntot);
    F(1:ntot, ntot+1:end) = [dis.red*FOI, FOI];

    vvec = [(dis.sig1+dis.sig2).*onesn;      (dis.g1 + p3).*onesn;       (g2 + h + p4).*onesn];%g2 and h are vectors
    
    V    = diag(vvec);
    V(ntot+1:2*ntot,1:ntot)   = diag(-dis.sig1.*onesn);
    V(2*ntot+1:3*ntot,1:ntot) = diag(-dis.sig2.*onesn);

    NGM = F/V;
    d = eigs(NGM,1);%largest in magnitude (+/-) 
    R = max(d); 

end

% get_R(data.ntot, dis, dis.h, dis.g2, data.NNs, data.NNs, data.basic_contact_matrix, 1, 1, 0, 0)
% 
% function Rt = rep_num(ntot,dis,h,g2,S,N,D,betamod,p3,p4)
%         
%     FOIu = repmat(S,1,ntot).*dis.beta.*dis.rr_infection.*betamod.*D./repmat(N',ntot,1);%+Shv1
%     
%     F                                = zeros(3*ntot,3*ntot);
%     F(1:ntot, 1*ntot+1:3*ntot) = [dis.red*FOIu,  FOIu]; %,  dis.red*(1-dis.trv1)*FOIu,  (1-dis.trv1)*FOIu
%     
%     onesn                            = ones(ntot,1);
%     vvec                             = [(dis.sig1+dis.sig2).*onesn;  (dis.g1+p3).*onesn;...
%                                         (g2+h+p4).*onesn;  ];  %    (dis.sig1+dis.sig2).*onesn;      (dis.g1+p3).*onesn;         (dis.g2_v1+dis.h_v1+p4).*onesn
%     V                                = diag(vvec);
%     V(1*ntot+1:2*ntot,1:ntot)        = diag(-dis.sig1.*onesn);
%     V(2*ntot+1:3*ntot,1:ntot)        = diag(-dis.sig2.*onesn);
%     
%     NGM = F/V;
%     Rt  = eigs(NGM,1);
%     
% end
% 
