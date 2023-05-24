function dis = p2Params_SARS

dis = struct;

%Seeding
% dis.t0 = -80;

%Probabilities
dis.ps  = 0.867;
dis.ihr = [0.0578 0.0578 0.0578 0.0578 ...	
           0.0816 0.0816 0.0816 0.0816 ...	
           0.3026 0.3026 0.3026 0.3026 ...	
           0.8670 0.8670 0.8670 0.8670 0.6018];
dis.ifr = dis.ps*[0.017 0.017 0.017 0.017 ...
                  0.024 0.024 0.024 0.024 ...
                  0.089 0.089 0.089 0.089 ...
                  0.255 0.255 0.255 0.255 0.177];

%Durations
dis.Tlat  = 4.6;
dis.Tay   = 2.1;
dis.Tsr   = 4.0;
% dis.Tsh   = 4.0;
% dis.Threc = 26.5-4.0;
% dis.Thd   = 23.7-4.0;
dis.Tsh   = 3.75;
dis.Threc = 26.5-3.75;
dis.Thd   = 23.7-3.75;
dis.Ti    = 365;

%Transmission
dis.red = 0.58;
dis.R0  = 1.7500;%3.0000;


load(sprintf('%sR0.mat','SARS'));
dis.R0values = R0;
dis.R0quantiles = qR0;

end