function Td = get_doubling_time(data,dis)

    ntot = data.ntot;
    onesn = ones(data.ntot,1);

    % doubling Time
    J                                  = zeros(7*ntot,7*ntot);
    J(1:ntot,2*ntot+1:3*ntot)          = -dis.beta*dis.red*data.basic_contact_matrix.*repmat(dis.rr_infection,1,ntot);
    J(1:ntot,3*ntot+1:4*ntot)          = -dis.beta*data.basic_contact_matrix.*repmat(dis.rr_infection,1,ntot);
    J(1:ntot,5*ntot+1:6*ntot)          = diag(onesn.*dis.nu);
    J(ntot+1:2*ntot,1*ntot+1:2*ntot)   = diag(onesn.*(-dis.sig1-dis.sig2));
    J(ntot+1:2*ntot,2*ntot+1:3*ntot)   = dis.beta*dis.red*data.basic_contact_matrix.*repmat(dis.rr_infection,1,ntot);
    J(ntot+1:2*ntot,3*ntot+1:4*ntot)   = dis.beta*data.basic_contact_matrix.*repmat(dis.rr_infection,1,ntot);
    J(2*ntot+1:3*ntot,1*ntot+1:2*ntot) = diag(onesn.*dis.sig1);
    J(2*ntot+1:3*ntot,2*ntot+1:3*ntot) = diag(onesn.*-dis.g1);
    J(3*ntot+1:4*ntot,1*ntot+1:2*ntot) = diag(onesn.*dis.sig2);
    J(3*ntot+1:4*ntot,3*ntot+1:4*ntot) = diag(onesn.*(-dis.g2-dis.h));
    J(4*ntot+1:5*ntot,3*ntot+1:4*ntot) = diag(onesn.*dis.h);
    J(4*ntot+1:5*ntot,4*ntot+1:5*ntot) = diag(onesn.*(-dis.g3-dis.mu));
    J(5*ntot+1:6*ntot,2*ntot+1:3*ntot) = diag(onesn.*dis.g1);
    J(5*ntot+1:6*ntot,3*ntot+1:4*ntot) = diag(onesn.*dis.g2);
    J(5*ntot+1:6*ntot,4*ntot+1:5*ntot) = diag(onesn.*dis.g3);
    J(5*ntot+1:6*ntot,5*ntot+1:6*ntot) = diag(onesn.*-dis.nu);
    J(6*ntot+1:7*ntot,4*ntot+1:5*ntot) = diag(onesn.*dis.mu);

    r       = max(real(eig(J)));
    Td      = log(2)/r;
end

