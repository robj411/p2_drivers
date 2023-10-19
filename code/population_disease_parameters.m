function dis = population_disease_parameters(data,dis,R0betafun)

%% COUNTRY PARAMETERS:

%% INITIAL DISEASE PARAMETERS:
%Population by Age
nn     = data.Npop';
nn     = [nn(1:16),sum(nn(17:end))];%last age range for disease is 80+
nntot  = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
ranges = [1,3,9,4];
nntot  = repelem(nntot,ranges);
nnprop = nn./nntot;
subs   = 1:4;
subs   = repelem(subs,ranges);
lx = data.lx;
adInd = data.adInd;


%Probabilities
phgs    = dis.ihr./dis.ps;
pdgh    = dis.ifr./dis.ihr;
dis.hfr = pdgh;
phgs    = accumarray(subs',phgs.*nnprop);
dis.ph  = [repmat(phgs(adInd),data.lx,1);phgs];
nnh     = nn.*dis.ihr;
nnhtot  = [nnh(1),sum(nnh(2:4)),sum(nnh(5:13)),sum(nnh(14:end))];
nnhtot  = repelem(nnhtot,ranges);
nnhprop = nnh./nnhtot;
pdgh    = accumarray(subs',pdgh.*nnhprop);
dis.pd  = [repmat(pdgh(adInd),data.lx,1);pdgh];

dis.rr_infection = [repmat(data.bmi_rr(1,1),lx,1); 1; 1; data.bmi_rr(1,1); data.bmi_rr(2,1)];

dis.ph([1:lx, lx+adInd]) = data.bmi_rr(1,2) * dis.ph([1:lx, lx+adInd]);
dis.pd([1:lx, lx+adInd]) = data.bmi_rr(1,3) * dis.pd([1:lx, lx+adInd]);

dis.ph(lx + data.adInd + 1) = data.bmi_rr(2,2) * dis.ph(lx + data.adInd + 1);
dis.pd(lx + data.adInd + 1) = data.bmi_rr(2,3) * dis.pd(lx + data.adInd + 1);

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

%% vaccination

%Vaccination: Broadly protective sarbecovirus vaccine (BPSV)
dis.hrv1 = 1/14;                       %time to develop v-acquired immunity
dis.scv1 = 0.35;                       %infection-blocking efficacy
heff1 = 0.80;                       %severe-disease-blocking efficacy
dis.hv1  = 1-((1-heff1)/(1-dis.scv1)); 
dis.trv1 = 0;%.52;                       %transmission-blocking efficacy
dis.nuv1 = 1/365/5;                      %duration of v-acquired immunity

dis.Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr)  +((1-dis.hv1)*dis.ph.*dis.Tsh);
dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./dis.Ts_v1;
dis.h_v1  = (1-dis.hv1)*dis.ph./dis.Ts_v1;

% SARS-X specific
% dis.hrv2 = 1/14;                       %time to develop v-acquired immunity
dis.scv2 = 0.55;                       %infection-blocking efficacy
heff2 = 0.90;                       %severe-disease-blocking efficacy
dis.hv2  = 1-((1-heff2)/(1-dis.scv2)); 
dis.trv2 = 0;%.52;                       %transmission-blocking efficacy
dis.nuv2 = 1/365/5;                      %duration of v-acquired immunity

dis.Ts_v2 = ((1-(1-dis.hv2)*dis.ph).*dis.Tsr)  +((1-dis.hv2)*dis.ph.*dis.Tsh);
dis.g2_v2 = (1-(1-dis.hv2)*dis.ph)./dis.Ts_v2;
dis.h_v2  = (1-dis.hv2)*dis.ph./dis.Ts_v2;

%% Transmission

% Deff  = data.basic_contact_matrix .* repmat(dis.rr_infection,1,data.ntot) .* repmat(data.NNs,1,data.ntot)./repmat(data.NNs',data.ntot,1);
% onesn = ones(data.ntot,1);
% F     = zeros(3*data.ntot,3*data.ntot);
% F(1:data.ntot,data.ntot+1:end)=[dis.red*Deff,Deff];
% 
% vvec = [(dis.sig1+dis.sig2).*onesn;      dis.g1.*onesn;       (dis.g2+dis.h).*onesn];%g2 and h are vectors
% V    = diag(vvec);
% V(data.ntot+1:2*data.ntot,1:data.ntot)   = diag(-dis.sig1.*onesn);
% V(2*data.ntot+1:3*data.ntot,1:data.ntot) = diag(-dis.sig2.*onesn);
% 
% GD = F/V;
% d = eigs(GD,1);%largest in magnitude (+/-) 
% R0a = max(d); 

zs = zeros(size(data.NNs));
dis.CI = get_R(data.ntot, dis, dis.h, dis.g2, data.NNs, zs, zs,...
    data.NNs, data.contacts.basic_contact_matrix, 1, 1, 0, 0); %R0a;


R0beta = R0betafun(dis);
dis.R0 = R0beta(1);
dis.beta = R0beta(2);

dis.Td = get_doubling_time(data,dis);

end
