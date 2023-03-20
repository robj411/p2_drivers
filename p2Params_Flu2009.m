function dis = p2Params_Flu2009

dis = struct;

%Seeding
% dis.t0 = -80;

%Probabilities
dis.ps  = 0.669;
dis.ihr = dis.ps*[0.00697 0.00274 0.00274 0.00274 ...
                  0.00274 0.00561 0.00561 0.00561 ...
                  0.00561 0.00561 0.01060 0.01060 ...
                  0.01060 0.01546 0.01546 0.01546 0.01546];
dis.ifr = dis.ps*[0.000276 0.00011 0.00011 0.00012 ...
                  0.00012 0.00030 0.00030 0.00030 ...
                  0.00030 0.00065 0.00065 0.00065 ...
                  0.00065 0.00980 0.00980 0.00980 0.00980];

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
dis.R0  = 1.5800;

end