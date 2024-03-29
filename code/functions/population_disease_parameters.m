% function to combine population and disease parameters to get
% within-country pathogen parameters
%
% data: struct of general model parameters
% dis: struct of pathogen parameters
% R0betafun: the function that computes beta from R0
%
% dis: struct of pathogen parameters

function dis = population_disease_parameters(data,dis,R0betafun)

%% COUNTRY PARAMETERS

%% INITIAL DISEASE PARAMETERS
%Population by Age
nn     = data.Npop';
nn     = [nn(1:16),sum(nn(17:end))];%last age range for disease is 80+
nntot  = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
ranges = [1,3,9,4];
nntot  = repelem(nntot,ranges);
nnprop = nn./nntot;
subs   = 1:4;
subs   = repelem(subs,ranges);
nSectors = data.nSectors;
adInd = data.adInd;


%Probabilities
phgs    = dis.ihr./dis.ps;
pdgh    = dis.ifr./dis.ihr;
dis.hfr = pdgh;
phgs    = accumarray(subs',phgs.*nnprop);
dis.ph  = [repmat(phgs(adInd),nSectors,1);phgs];
nnh     = nn.*dis.ihr;
nnhtot  = [nnh(1),sum(nnh(2:4)),sum(nnh(5:13)),sum(nnh(14:end))];
nnhtot  = repelem(nnhtot,ranges);
nnhprop = nnh./nnhtot;
pdgh    = accumarray(subs',pdgh.*nnhprop);
dis.pd  = [repmat(pdgh(adInd),nSectors,1);pdgh];

dis.rr_infection = 1; % [repmat(data.bmi_rr(1,1),nSectors,1); 1; 1; data.bmi_rr(1,1); data.bmi_rr(2,1)];

dis.ph([1:nSectors, nSectors+adInd]) = data.bmi_rr(1,2) * dis.ph([1:nSectors, nSectors+adInd]);
dis.pd([1:nSectors, nSectors+adInd]) = data.bmi_rr(1,3) * dis.pd([1:nSectors, nSectors+adInd]);

dis.ph(nSectors + adInd + 1) = data.bmi_rr(2,2) * dis.ph(nSectors + adInd + 1);
dis.pd(nSectors + adInd + 1) = data.bmi_rr(2,3) * dis.pd(nSectors + adInd + 1);

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

%% vaccines

%Vaccine: Broadly protective sarbecovirus vaccine (BPSV)
dis.hrv1 = 1/21;                       %time to develop v-acquired immunity
dis.scv1 = 0.35;                       %infection-blocking effectiveness
heff1 = 0.80;                       %severe-disease-blocking effectiveness
dis.hv1  = 1-((1-heff1)/(1-dis.scv1)); 
dis.trv1 = 0;%.52;                       %transmission-blocking effectiveness
dis.nuv1 = 1/365000000; %365/5;                      %duration of v-acquired immunity

Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr)  +((1-dis.hv1)*dis.ph.*dis.Tsh);
dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./Ts_v1;
dis.h_v1  = (1-dis.hv1)*dis.ph./Ts_v1;

% SARS-X specific
dis.hrv2 = 1/21;                       %time to develop v-acquired immunity
dis.scv2 = 0.55;                       %infection-blocking effectiveness
heff2 = 0.90;                       %severe-disease-blocking effectiveness
dis.hv2  = 1-((1-heff2)/(1-dis.scv2)); 
dis.trv2 = 0;                       %transmission-blocking effectiveness
dis.nuv2 = 1/365000000;                     %duration of v-acquired immunity

Ts_v2 = ((1-(1-dis.hv2)*dis.ph).*dis.Tsr) + ((1-dis.hv2)*dis.ph.*dis.Tsh);
dis.g2_v2 = (1-(1-dis.hv2)*dis.ph)./Ts_v2;
dis.h_v2  = (1-dis.hv2)*dis.ph./Ts_v2;

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

NNs = data.NNs;
zs = zeros(size(NNs));
dis.CI = get_candidate_infectees(length(NNs), dis, NNs,zs, zs, 0, 0, NNs, data.contacts.basic_contact_matrix);

R0beta = R0betafun(dis);
dis.R0 = R0beta(1);
dis.beta = R0beta(2);

dis.Td = get_doubling_time(data,dis);
% dis.time_to_5 = log2(5) * dis.Td;

dis.generation_time = log(dis.R0) / (log(2) / dis.Td);

% % get upper Rs in different configurations
% configurations = [data.x_elim, data.x_schc, data.x_econ];
% NNbar                = data.NNs;
% R0s = [];
% zs = zeros(size(data.NNs));
% for i=1:size(configurations,2)
%     Dtemp   = p2MakeDs(data,NNbar,configurations(:,i),data.wfh(1,:));
%     R0s(i) = get_R(data.nStrata, dis, NNbar, zs, zs, dis.beta, 0, 0, 0, data, i, Dtemp);
% end
% dis.R0s = R0s;

end
