function [CIa, beta, R0] = p2getCandidateInfectees(data,dis)

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

%Probabilities
phgs    = dis.ihr./dis.ps;
phgs    = accumarray(subs',phgs.*nnprop);
dis.ph  = [repmat(phgs(adInd),lx,1);phgs];

%Durations
dis.Ts = ((1-dis.ph).*dis.Tsr)   + (dis.ph.*dis.Tsh);

%Rates
dis.sig1 = (1-dis.ps)/dis.Tlat;
dis.sig2 = dis.ps/dis.Tlat;
dis.g1   = 1/dis.Tay;
dis.g2   = (1-dis.ph)./dis.Ts;
dis.h    = dis.ph./dis.Ts;

%Transmission
Deff  = Dout.*repmat(data.NNs,1,ntot)./repmat(data.NNs',ntot,1);
onesn = ones(ntot,1);
F     = zeros(3*ntot,3*ntot);
F(1:ntot,ntot+1:end)=[dis.red*Deff,Deff];

vvec = [(dis.sig1+dis.sig2).*onesn;      dis.g1.*onesn;       (dis.g2+dis.h).*onesn];%g2 and h are vectors
V    = diag(vvec);
V(ntot+1:2*ntot,1:ntot)   = diag(-dis.sig1.*onesn);
V(2*ntot+1:3*ntot,1:ntot) = diag(-dis.sig2.*onesn);

GD = F/V;
d = eigs(GD,1);%largest in magnitude (+/-) 
CIa = max(d); 

R0 = normrnd(dis.R0,dis.R0/10,1,1);
beta = R0/CIa;

end