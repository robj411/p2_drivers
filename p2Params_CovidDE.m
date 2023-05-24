function dis = p2Params_CovidDE

dis = struct;

%Seeding
% dis.t0 = -80;

%Probabilities
dis.ps  = 0.595;
dis.ihr = 1.85*[0.000016 0.000016 0.000408 0.000408 ...	
          0.010400 0.010400 0.034300 0.034300 ...	
          0.042500 0.042500 0.081600 0.081600 ...	
          0.118000 0.118000 0.166000 0.166000 0.184000];
dis.ifr = 1.85*[0.000016 0.000016 0.000070 0.000070 ...
          0.000309 0.000309 0.000844 0.000844 ...
          0.001610 0.001610 0.005950 0.005950 ...
          0.019300 0.019300 0.042800 0.042800 0.078000];

%Durations
dis.Tlat  = 4.0;
dis.Tay   = 2.1;
dis.Tsr   = 4.0;
dis.Tsh   = 4.0;
dis.Threc = 7.6;
dis.Thd   = 7.6;
dis.Ti    = 365;

%Transmission
dis.red = 0.58;
dis.R0  = 5.0800;


load(sprintf('%sR0.mat','Covid Delta'));
dis.R0values = R0;
dis.R0quantiles = qR0;

end