function [data,dis,p2] = p2Params(data,inp2,inp4,inp5,inp6)
%% COUNTRY PARAMETERS:

%Population by Age
nn     = data.Npop';
nn     = [nn(1:16),sum(nn(17:end))];%last age range is 80+
nntot  = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
ranges = [1,3,9,4];
nntot  = repelem(nntot,ranges);
nnprop = nn./nntot;
subs   = 1:4;
subs   = repelem(subs,ranges);

%Population by Sector
ntot                  = size(data.NNs,1);
data.NNs(data.NNs==0) = 1;

%Age-Sector Breakdown
adInd = 3;
lx    = length(data.obj);

data.alp = 1;

[Dout,data] = p2MakeDs(data,data.NNs,ones(lx,1),zeros(1,lx));

%% INITIAL DISEASE PARAMETERS:
    
if strcmp(inp2,'Influenza 2009');
    dis = p2Params_Flu2009;
elseif strcmp(inp2,'Influenza 1957');
    dis = p2Params_Flu1957;
elseif strcmp(inp2,'Influenza 1918');
    dis = p2Params_Flu1918;
elseif strcmp(inp2,'Covid Wildtype');
    dis = p2Params_CovidWT;    
elseif strcmp(inp2,'Covid Omicron');
    dis = p2Params_CovidOM;    
elseif strcmp(inp2,'Covid Delta');
    dis = p2Params_CovidDE;    
elseif strcmp(inp2,'SARS');
    dis = p2Params_SARS;
else
    error('Unknown Disease!');
end   

dis.R0  = inp5*dis.R0;
dis.ps  = min(inp6*dis.ps ,1);
dis.ihr = min(inp6*dis.ihr,dis.ps);
dis.ifr = min(inp6*dis.ifr,dis.ihr);

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

%Calculations
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
dis.beta=dis.R0/R0a;%beta scales actual R0 to fitted R0

%Vaccination
dis.hrv1 = 1/28;                       %time to develop v-acquired immunity
dis.scv1 = 0.60;                       %infection-blocking efficacy
dis.heff = 0.87;                       %severe-disease-blocking efficacy
dis.hv1  = 1-((1-dis.heff)/(1-dis.scv1)); 
dis.trv1 = 0.52;                       %transmission-blocking efficacy
dis.nuv1 = 1/365;                      %duration of v-acquired immunity

dis.Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr)  +((1-dis.hv1)*dis.ph.*dis.Tsh);
dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./dis.Ts_v1;
dis.h_v1  = (1-dis.hv1)*dis.ph./dis.Ts_v1;

%% PREPAREDNESS PARAMETERS:

p2 = struct;

if strcmp(inp4,'CURRENT');
    p2.Tres  = data.Tres;                       %Response Time
    p2.t_tit = data.t_tit;                      %Test-Isolate-Trace Time
    p2.trate = data.trate;                      %Test-Isolate-Trace Rate
    p2.sdl   = data.sdl;                        %Social Distancing Asymptote
    p2.sdb   = data.sdb;                        %Social Distancing Steepness
    p2.Hmax  = data.Hmax*sum(data.Npop)/(10^5); %Hospital Capacity
    t_vax    = data.t_vax;                      %Vaccine Administration Time
    arate    = data.arate*sum(data.Npop/10^5);  %Vaccine Administration Rate
    puptake  = data.puptake;                    %Vaccine Uptake
    
elseif strcmp(inp4,'RANDOM');
    p2.Tres  =      random(truncate(makedist('normal','mu',66.1964,'sigma',12.3309),-75,Inf));
    p2.t_tit =      random(truncate(makedist('lognormal','mu',4.12346,'sigma',0.565271),-75,Inf));
    p2.trate =      random(truncate(makedist('lognormal','mu',5.16249,'sigma',1.17374),0,Inf));   
    p2.sdl   =      random(truncate(makedist('uniform','lower',0.1,'upper',0.4),0.1,0.4)); 
    p2.sdb   =  exp(random(truncate(makedist('normal','mu',-12.8411,'sigma',10.8119),-Inf,Inf)));%%%
    p2.Hmax  =      random(truncate(makedist('lognormal','mu',4.12546,'sigma',0.836187),0,Inf))*sum(data.Npop)/(10^5);
    t_vax    =      random(truncate(makedist('normal','mu',450.916,'sigma',57.0053),-75,Inf));
    arate    =      random(truncate(makedist('normal','mu',329.448,'sigma',112.91),0,Inf))*sum(data.Npop)/(10^5);
    puptake  =      random(truncate(makedist('normal','mu',0.70855,'sigma',0.142042),0,1));
    %
    p2.t_vax   = t_vax;
    p2.arate   = arate;
    p2.puptake = puptake;
    %    
    % clist      = dir(fullfile('*.mat'));
    % clist      = {clist.name};
    % clist(end) = [];
    % X          = zeros(length(clist),9);
    % for i = 1:length(clist);
    %     load(clist{i},'data');  
    %     X(i,1) = data.Tres;
    %     X(i,2) = data.t_tit;
    %     X(i,3) = data.trate;
    %     X(i,4) = data.sdl;
    %     X(i,5) = data.sdb;
    %     X(i,6) = data.Hmax;
    %     X(i,7) = data.t_vax;
    %     X(i,8) = data.arate;
    %     X(i,9) = data.puptake;   
    % end
    % 
    % var = X(:,9);
    % lim = [0,1];
    % figure;
    % histogram((var),length(clist));
    % 
    % p   = fitdist((var),'normal');
    % figure;
    % hold on;
    % sup = linspace(min(var),max(var),1000);
    % plot(sup,pdf(p,sup),'r');
    % 
    % pt = truncate(p,lim(1),lim(2));
    % plot(sup,pdf(pt,sup),'g');
    
elseif strcmp(inp4,'LEVEL1');
    p2.Tres  = 79.000000000000000;
    p2.t_tit = 245.3750;
    p2.trate = 4.2808;   
    p2.sdl   = 0.399999999999975;
    p2.sdb   = 1.011754387252751;
    p2.Hmax  = 10.5000*sum(data.Npop)/(10^5);
    t_vax    = 603.9147444001774;
    arate    = 49.9622236620684*sum(data.Npop/10^5);
    puptake  = 0.11494999999999999;
    
elseif strcmp(inp4,'LEVEL2');
    p2.Tres  = 74.750000000000000;
    p2.t_tit = 88.0000;
    p2.trate = 37.4892;   
    p2.sdl   = 0.399999999940825;
    p2.sdb   = 0.003603254435597;
    p2.Hmax  = 22.0500*sum(data.Npop)/(10^5);
    t_vax    = 524.0138662541277;
    arate    = 165.8753601340584*sum(data.Npop/10^5);
    puptake  = 0.40202500000000001;   
    
elseif strcmp(inp4,'LEVEL3');
    p2.Tres  = 71.000000000000000;
    p2.t_tit = 62.0000;
    p2.trate = 126.9593;   
    p2.sdl   = 0.317269065393003;
    p2.sdb   = 0.000008975136278;
    p2.Hmax  = 46.5750*sum(data.Npop)/(10^5);
    t_vax    = 461.0510305161198;
    arate    = 280.2617158912323*sum(data.Npop/10^5);
    puptake  = 0.66097499999999997;
    
elseif strcmp(inp4,'LEVEL4');
    p2.Tres  = 61.000000000000000;
    p2.t_tit = 43.0000;
    p2.trate = 373.7261;   
    p2.sdl   = 0.100000000021886;
    p2.sdb   = 0.000000000000038;
    p2.Hmax  = 137.7730*sum(data.Npop)/(10^5);
    t_vax    = 417.7319167406767;
    arate    = 437.5090847154494*sum(data.Npop/10^5);
    puptake  = 0.82464999999999989;

else
    error('Unknown Preparedness!');    
    
end

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
p2.qg2_v1 = (1-(1-dis.hv1)*dis.ph)./(dis.Ts_v1-p2.dur);
p2.qh     = dis.ph./(dis.Ts-p2.dur);
p2.qh_v1  = (1-dis.hv1)*dis.ph./(dis.Ts_v1-p2.dur);

%Hospital Capacity
p2.thl   = max(1,0.25*p2.Hmax);%lower threshold can't be less than 1 occupant
p2.Hmax  = max(4*p2.thl,p2.Hmax);
p2.SHmax = 2*p2.Hmax;

%Vaccine Uptake
Npop    = data.Npop;
NNage   = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
puptake = min(0.99*(1-NNage(1)/sum(NNage)),puptake);%population uptake cannot be greater than full coverage in non-pre-school age groups
up3fun  = @(u3) puptake*sum(NNage) - u3*(NNage(2)/2 + NNage(3)) - min(1.5*u3,1)*NNage(4);
if up3fun(0)*up3fun(1)<=0;
    u3  = fzero(up3fun,[0 1]);
else
    u3  = fminbnd(up3fun,0,1);
end
u4      = min(1.5*u3,1);
u1      = 0;
up2fun  = @(u2) u2*NNage(2) + u3*NNage(3) + u4*NNage(4) - puptake*sum(NNage);
u2      = fzero(up2fun,[0 1]);
uptake  = [u1,u2,u3,u4];
% if ((uptake*NNage'/sum(NNage))-puptake)~=0;
%     error('Vaccine uptake error!');
% end

%Vaccine Administration Rate
t_ages     = min((uptake.*NNage)/arate,Inf);%arate may be 0
if strcmp(inp2,'Influenza 1918');
    t_ages     = [t_ages(3),t_ages(4),t_ages(2),t_ages(1)];
    p2.aratep1 = [0;0;arate;0];%Period 1 - working-age%to be split across all economic sectors in heSimCovid19vax.m
    p2.aratep2 = [0;0;0;arate];%Period 2 - retired-age
    p2.aratep3 = [0;arate;0;0];%Period 3 - school-age
    p2.aratep4 = [0;0;0;0];    %Period 4 - pre-school-age
else
    t_ages     = [t_ages(4),t_ages(3),t_ages(2),t_ages(1)];
    p2.aratep1 = [0;0;0;arate];%Period 1 - retired-age
    p2.aratep2 = [0;0;arate;0];%Period 2 - working-age%to be split across all economic sectors in heSimCovid19vax.m
    p2.aratep3 = [0;arate;0;0];%Period 3 - school-age
    p2.aratep4 = [0;0;0;0];    %Period 4 - pre-school-age
end    
tpoints    = cumsum([t_vax,t_ages]);
p2.startp1 = tpoints(1);
p2.startp2 = tpoints(2);
p2.startp3 = tpoints(3);
p2.startp4 = tpoints(4);
p2.end     = tpoints(5);%End of Rollout

end