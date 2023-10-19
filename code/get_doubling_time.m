function Td = get_doubling_time(data,dis)

    ntot = data.ntot;
    onesn = ones(data.ntot,1);
    basic_contact_matrix = data.contacts.basic_contact_matrix;
    
    rr_mat = repmat(dis.rr_infection,1,ntot);

    % doubling Time
    % 1: S
    % 2: E
    % 3: Ia
    % 4: Is
    % 5: H
    % 6: R
    % 7: D
    J                                  = zeros(7*ntot,7*ntot);
    J(1:ntot, 2*ntot+1:3*ntot)          = -dis.beta .* dis.red .* basic_contact_matrix .* rr_mat;
    J(1:ntot, 3*ntot+1:4*ntot)          = -dis.beta .* basic_contact_matrix .* rr_mat;
    J(1:ntot, 5*ntot+1:6*ntot)          = diag(onesn .* dis.nu);
    J(ntot+1:2*ntot, 1*ntot+1:2*ntot)   = diag(onesn .* (-dis.sig1-dis.sig2));
    J(ntot+1:2*ntot, 2*ntot+1:3*ntot)   = dis.beta .* dis.red .* basic_contact_matrix .* rr_mat;
    J(ntot+1:2*ntot, 3*ntot+1:4*ntot)   = dis.beta .* basic_contact_matrix .* rr_mat;
    J(2*ntot+1:3*ntot, 1*ntot+1:2*ntot) = diag(onesn .* dis.sig1);
    J(2*ntot+1:3*ntot, 2*ntot+1:3*ntot) = diag(onesn .* -dis.g1);
    J(3*ntot+1:4*ntot, 1*ntot+1:2*ntot) = diag(onesn .* dis.sig2);
    J(3*ntot+1:4*ntot, 3*ntot+1:4*ntot) = diag(onesn .* (-dis.g2-dis.h));
    J(4*ntot+1:5*ntot, 3*ntot+1:4*ntot) = diag(onesn .* dis.h);
    J(4*ntot+1:5*ntot, 4*ntot+1:5*ntot) = diag(onesn .* (-dis.g3-dis.mu));
    J(5*ntot+1:6*ntot, 2*ntot+1:3*ntot) = diag(onesn .* dis.g1);
    J(5*ntot+1:6*ntot, 3*ntot+1:4*ntot) = diag(onesn .* dis.g2);
    J(5*ntot+1:6*ntot, 4*ntot+1:5*ntot) = diag(onesn .* dis.g3);
    J(5*ntot+1:6*ntot, 5*ntot+1:6*ntot) = diag(onesn .* -dis.nu);
    J(6*ntot+1:7*ntot, 4*ntot+1:5*ntot) = diag(onesn .* dis.mu);

    r       = max(real(eig(J)));
    Td      = log(2)/r;
end

