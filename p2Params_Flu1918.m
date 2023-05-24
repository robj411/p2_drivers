function dis = p2Params_Flu1918

dis = struct;

%Seeding
% dis.t0 = -80;

%Probabilities
dis.ps  = 0.669;
dis.ihr = dis.ps*8*[0.02284 0.00398 0.00478 0.00983 ...
                    0.01700 0.02922 0.02470 0.02205 ...
                    0.01647 0.01195 0.01647 0.01169 ...
                    0.03081 0.04144 0.04941 0.04941 0.04941];
dis.ifr = dis.ps*[0.02284 0.00398 0.00478 0.00983 ...
                  0.01700 0.02922 0.02470 0.02205 ...
                  0.01647 0.01195 0.01647 0.01169 ...
                  0.03081 0.04144 0.04941 0.04941 0.04941];

%Durations
dis.Tlat  = 1.1;
% dis.Tay   = 2.1;
% dis.Tsr   = 4.0;
% dis.Tsh   = 4.0;
dis.Tay   = 2.5;
dis.Tsr   = 2.5;
dis.Tsh   = 2.5;
dis.Threc = 5.0;
dis.Thd   = 5.0;
dis.Ti    = 365;

%Transmission
dis.red = 0.58;
dis.R0  = 2.5000;

load(sprintf('%sR0.mat','Influenza 1918'));
dis.R0values = R0;
dis.R0quantiles = qR0;

end