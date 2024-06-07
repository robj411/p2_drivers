% computes time for number of infections to double

% dis: struct of pathogen parameters
% data: struct of general model parameters

% Td: doubling time

function Td = get_doubling_time(data,dis)

    nStrata = data.nStrata;
    onesn = ones(nStrata,1);
    basic_contact_matrix = data.contacts.basic_contact_matrix;
    
    rr_mat = repmat(dis.rr_infection,1,nStrata);

    % doubling Time
    % 1: S
    % 2: E
    % 3: Ia
    % 4: Is
    % 5: H
    % 6: R
    % 7: D
    J                                  = zeros(7*nStrata,7*nStrata);
    J(1:nStrata, 2*nStrata+1:3*nStrata)          = -dis.beta .* dis.red .* basic_contact_matrix .* rr_mat;
    J(1:nStrata, 3*nStrata+1:4*nStrata)          = -dis.beta .* basic_contact_matrix .* rr_mat;
    J(1:nStrata, 5*nStrata+1:6*nStrata)          = diag(onesn .* dis.nu);
    J(nStrata+1:2*nStrata, 1*nStrata+1:2*nStrata)   = diag(onesn .* (-dis.sig1-dis.sig2));
    J(nStrata+1:2*nStrata, 2*nStrata+1:3*nStrata)   = dis.beta .* dis.red .* basic_contact_matrix .* rr_mat;
    J(nStrata+1:2*nStrata, 3*nStrata+1:4*nStrata)   = dis.beta .* basic_contact_matrix .* rr_mat;
    J(2*nStrata+1:3*nStrata, 1*nStrata+1:2*nStrata) = diag(onesn .* dis.sig1);
    J(2*nStrata+1:3*nStrata, 2*nStrata+1:3*nStrata) = diag(onesn .* -dis.g1);
    J(3*nStrata+1:4*nStrata, 1*nStrata+1:2*nStrata) = diag(onesn .* dis.sig2);
    J(3*nStrata+1:4*nStrata, 3*nStrata+1:4*nStrata) = diag(onesn .* (-dis.g2-dis.h));
    J(4*nStrata+1:5*nStrata, 3*nStrata+1:4*nStrata) = diag(onesn .* dis.h);
    J(4*nStrata+1:5*nStrata, 4*nStrata+1:5*nStrata) = diag(onesn .* (-dis.g3-dis.mu));
    J(5*nStrata+1:6*nStrata, 2*nStrata+1:3*nStrata) = diag(onesn .* dis.g1);
    J(5*nStrata+1:6*nStrata, 3*nStrata+1:4*nStrata) = diag(onesn .* dis.g2);
    J(5*nStrata+1:6*nStrata, 4*nStrata+1:5*nStrata) = diag(onesn .* dis.g3);
    J(5*nStrata+1:6*nStrata, 5*nStrata+1:6*nStrata) = diag(onesn .* -dis.nu);
    J(6*nStrata+1:7*nStrata, 4*nStrata+1:5*nStrata) = diag(onesn .* dis.mu);

    r       = max(real(eig(J)));
    Td      = log(2)/r;
end




