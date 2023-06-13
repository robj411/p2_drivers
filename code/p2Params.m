function [data,dis,p2] = p2Params(data,inp2,dis)

%% COUNTRY PARAMETERS:

%Population by Age
nn     = data.Npop';
nn     = [nn(1:16),sum(nn(17:end))];%last age range for disease is 80+
nntot  = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
ranges = [1,3,9,4];
nntot  = repelem(nntot,ranges);
nnprop = nn./nntot;
subs   = 1:4;
subs   = repelem(subs,ranges);

%Population by Sector
adInd    = 3;
lx       = length(data.obj);
ntot     = size(data.NNs,1);
data.NNs(data.NNs==0) = 1;
data.alp = 1;

%Contact Matrix
[Dout,data] = p2MakeDs(data,data.NNs,ones(lx,1),zeros(1,lx));

%% INITIAL DISEASE PARAMETERS:
    
% dis.R0  = inp5*dis.R0;
% dis.ps  = min(inp6*dis.ps ,1);
% dis.ihr = min(inp6*dis.ihr,dis.ps);
% dis.ifr = min(inp6*dis.ifr,dis.ihr);

%Probabilities
phgs    = dis.ihr./dis.ps;
pdgh    = dis.ifr./dis.ihr;
phgs    = accumarray(subs',phgs.*nnprop);
dis.ph  = [repmat(phgs(adInd),lx,1);phgs];
nnh     = nn.*dis.ihr;
nnhtot  = [nnh(1),sum(nnh(2:4)),sum(nnh(5:13)),sum(nnh(14:end))];
nnhtot  = repelem(nnhtot,ranges);
nnhprop = nnh./nnhtot;
pdgh    = accumarray(subs',pdgh.*nnhprop);
dis.pd  = [repmat(pdgh(adInd),lx,1);pdgh];

%Durations
dis.Ts = ((1-dis.ph).*dis.Tsr)   + (dis.ph.*dis.Tsh);
dis.Th = ((1-dis.pd).*dis.Threc) + (dis.pd.*dis.Thd);

%Rates
dis.sig1 = (1-dis.ps)/dis.Tlat;
dis.sig2 = dis.ps/dis.Tlat;
dis.g1   = 1/dis.Tay;
dis.g2   = (1-dis.ph)./dis.Ts;
dis.g3   = (1-dis.pd)./dis.Th;
dis.h    = dis.ph./dis.Ts;
dis.mu   = dis.pd./dis.Th;
dis.nu   = 1/dis.Ti;

%Transmission
Deff  = Dout.*repmat(data.NNs,1,ntot)./repmat(data.NNs',ntot,1);
onesn = ones(ntot,1);
F     = zeros(3*ntot,3*ntot);
F(1:ntot,ntot+1:end)=[dis.red*Deff,Deff];

vvec = [(dis.sig1+dis.sig2).*onesn;      dis.g1.*onesn;       (dis.g2+dis.h).*onesn];%g2 and h are vectors
V    = diag(vvec);
V(ntot+1:2*ntot,1:ntot)   = diag(-dis.sig1.*onesn);
V(2*ntot+1:3*ntot,1:ntot) = diag(-dis.sig2.*onesn);

GD=F/V;
d=eigs(GD,1);%largest in magnitude (+/-) 
R0a=max(d); 

% get R0a as quantile from raw null distribution: quantiles qR0, values R0
if R0a > max(dis.R0values)
    R0quantile = max(dis.R0quantiles(dis.R0quantiles<1));
    disp('Too many candidate infectees');
else
    R0quantile = interp1(dis.R0values,dis.R0quantiles,R0a);
end

dis.R0sample = norminv(R0quantile,dis.R0,.1*dis.R0);

dis.beta = dis.R0sample/R0a;%beta scales actual R0 to fitted R0

%Vaccination
% dis.hrv1 = 1/28;                       %time to develop v-acquired immunity
% dis.scv1 = 0.60;                       %infection-blocking efficacy
dis.heff = 0.87;                       %severe-disease-blocking efficacy
% dis.hv1  = 1-((1-dis.heff)/(1-dis.scv1)); 
% dis.trv1 = 0.52;                       %transmission-blocking efficacy
% dis.nuv1 = 1/365;                      %duration of v-acquired immunity
% 
% dis.Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr)  +((1-dis.hv1)*dis.ph.*dis.Tsh);
% dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./dis.Ts_v1;
% dis.h_v1  = (1-dis.hv1)*dis.ph./dis.Ts_v1;

%% PREPAREDNESS PARAMETERS:

p2 = struct;

p2.Tres  = data.Tres;                       %Response Time
p2.t_tit = data.t_tit;                      %Test-Isolate-Trace Time
p2.trate = data.trate;                      %Test-Isolate-Trace Rate
p2.sdl   = data.sdl;                        %Social Distancing Asymptote
p2.sdb   = data.sdb;                        %Social Distancing Steepness
p2.Hmax  = data.Hmax*sum(data.Npop)/10^5;   %Hospital Capacity
% t_vax    = data.t_vax;                      %Vaccine Administration Time
% arate    = data.arate*sum(data.Npop/10^5);  %Vaccine Administration Rate
% puptake  = data.puptake;                    %Vaccine Uptake

%Response Time
J                                  = zeros(7*ntot,7*ntot);
J(1:ntot,2*ntot+1:3*ntot)          = -dis.beta*dis.red*Dout;
J(1:ntot,3*ntot+1:4*ntot)          = -dis.beta*Dout;
J(1:ntot,5*ntot+1:6*ntot)          = diag(onesn.*dis.nu);
J(ntot+1:2*ntot,1*ntot+1:2*ntot)   = diag(onesn.*(-dis.sig1-dis.sig2));
J(ntot+1:2*ntot,2*ntot+1:3*ntot)   = dis.beta*dis.red*Dout;
J(ntot+1:2*ntot,3*ntot+1:4*ntot)   = dis.beta*Dout;
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
if ~isfield(data,'Td_CWT');
    data.Td_CWT = Td;
end
p2.Tres = data.tvec(1) + (-data.tvec(1) + p2.Tres)*Td/data.Td_CWT;

%Test-Isolate-Trace
p2.t_tit  = data.tvec(1) + (-data.tvec(1) + p2.t_tit)*Td/data.Td_CWT;
p2.dur    = 1;
p2.qg1    = 1/(dis.Tay-p2.dur);
p2.qg2    = (1-dis.ph)./(dis.Ts-p2.dur);
% p2.qg2_v1 = (1-(1-dis.hv1)*dis.ph)./(dis.Ts_v1-p2.dur);
p2.qh     = dis.ph./(dis.Ts-p2.dur);
% p2.qh_v1  = (1-dis.hv1)*dis.ph./(dis.Ts_v1-p2.dur);

%Hospital Capacity
p2.thl   = max(1,0.25*p2.Hmax);%lower threshold can't be less than 1 occupant
p2.Hmax  = max(4*p2.thl,p2.Hmax);
p2.SHmax = 2*p2.Hmax;

%Vaccine Uptake
Npop    = data.Npop;
NNage   = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
% puptake = min(0.99*(1-NNage(1)/sum(NNage)),puptake);%population uptake cannot be greater than full coverage in non-pre-school age groups
% up3fun  = @(u3) puptake*sum(NNage) - u3*(NNage(2)/2 + NNage(3)) - min(1.5*u3,1)*NNage(4);
% if up3fun(0)*up3fun(1)<=0;
%     u3  = fzero(up3fun,[0 1]);
% else
%     u3  = fminbnd(up3fun,0,1);
% end
% u4      = min(1.5*u3,1);
% u1      = 0;
% up2fun  = @(u2) u2*NNage(2) + u3*NNage(3) + u4*NNage(4) - puptake*sum(NNage);
% u2      = fzero(up2fun,[0 1]);
% uptake  = [u1,u2,u3,u4];
% if ((uptake*NNage'/sum(NNage))-puptake)~=0;
%     error('Vaccine uptake error!');
% end

%Vaccine Administration Rate
% t_ages     = min((uptake.*NNage)/arate,Inf);%arate may be 0
% if strcmp(inp2,'Influenza 1918');
%     t_ages     = [t_ages(3),t_ages(4),t_ages(2),t_ages(1)];
%     p2.aratep1 = [0;0;arate;0];%Period 1 - working-age%to be split across all economic sectors in heSimCovid19vax.m
%     p2.aratep2 = [0;0;0;arate];%Period 2 - retired-age
%     p2.aratep3 = [0;arate;0;0];%Period 3 - school-age
%     p2.aratep4 = [0;0;0;0];    %Period 4 - pre-school-age
% else
%     t_ages     = [t_ages(4),t_ages(3),t_ages(2),t_ages(1)];
%     p2.aratep1 = [0;0;0;arate];%Period 1 - retired-age
%     p2.aratep2 = [0;0;arate;0];%Period 2 - working-age%to be split across all economic sectors in heSimCovid19vax.m
%     p2.aratep3 = [0;arate;0;0];%Period 3 - school-age
%     p2.aratep4 = [0;0;0;0];    %Period 4 - pre-school-age
% end    
% tpoints    = cumsum([t_vax,t_ages]);
% p2.startp1 = tpoints(1);
% p2.startp2 = tpoints(2);
% p2.startp3 = tpoints(3);
% p2.startp4 = tpoints(4);
p2.end     = 5000; %tpoints(5);%End of Rollout

%% COST PARAMETERS:

na         = [data.Npop(1:16)',sum(data.Npop(17:end))];%length is 17 to match ifr
la         = [data.la(1:16),...
              dot(data.la(17:end),[data.Npop(17),sum(data.Npop(18:end))])/sum(data.Npop(17:end)+1e-10)];
napd       = na.*dis.ifr;
lg         = [dot(la(1),napd(1))/sum(napd(1)),...
              dot(la(2:4),napd(2:4))/sum(napd(2:4)),...
              dot(la(5:13),napd(5:13))/sum(napd(5:13)),...
              dot(la(14:end),napd(14:end))/sum(napd(14:end))];
for k = 1:length(lg); 
    lgh(k) = sum(1./((1+0.03).^[1:lg(k)]));
end  
data.lgh   = [repmat(lgh(3),1,45),lgh];

end