% correlate_random_var: function to generate random quantiles correlated to
% a pre-existing vector via a multivariate normal with specified
% correlation
% 
% input_quant: column vector of input quantiles
% correlation: value between 0 and 1
%
% output_quant: column vector of correlated output quantiles

function output_quant = correlate_random_var(input_quant, correlation)


    n     = length(input_quant);                    % length of vector
    rho   = correlation;                   % desired correlation = cos(angle)
    theta = acos(rho);             % corresponding angle
    % map input quant to normal var
    x1    = norminv(input_quant, 0, 1);        % fixed given data
    x2    = normrnd(0, 1, n, 1);      % new random data
    Xctr     = [x1, x2];         % matrix

    Id   = eye(n);                               % identity matrix
    [Q, ~] = qr(Xctr(: , 1));
    P    = Q * Q';       % projection onto space defined by x1
    x2o  = (Id-P) * Xctr(: , 2);                 % x2ctr made orthogonal to x1ctr
    Xc2  = [Xctr(: , 1), x2o];               % bind to matrix
    Y    = Xc2 * diag(1./sqrt(sum(Xc2.^2)));  % scale columns to length 1

    x = Y(: , 2) + (1 / tan(theta)) * Y(: , 1);     % final new vector
    corr(x1, x);                                    % check correlation = rho
    output_quant = normcdf(x,0,1);

end

