function [data,f,g]=p2Run(data,dis,inp3,int,Xit,p2);

adInd = 3;
lx    = length(data.obj);
ln    = length(data.NNs);

NNbar                = data.NNs;
XitMat               = reshape(Xit,lx,int);
WitMat               = XitMat.^(1/data.alp);
WitMat(data.EdInd,:) = XitMat(data.EdInd,:);
NNvec                = repmat(NNbar(1:lx),1,int).*WitMat;
NNworkSum            = sum(NNvec,1);
NNvec(lx+1:ln,:)     = repmat(NNbar(lx+1:ln),1,int);
NNvec(lx+adInd,:)    = sum(NNbar([1:lx,lx+adInd]))-NNworkSum;
data.NNvec           = NNvec;

Dvec = zeros(ln,ln,int);
for i = 1:int;
    [Dtemp,~]   = p2MakeDs(data,NNvec(:,i),XitMat(:,i),data.hw(i,:));
    Dvec(:,:,i) = Dtemp;
end
data.Dvec = Dvec;

[data,f,g] = p2SimVax(data,NNvec,Dvec,dis,NNvec(:,1),inp3,WitMat,p2);

end

%%

function [data,f,g]=p2SimVax(data,NNvec,Dvec,dis,S0,inp3,WitMat,p2)               
%% PARAMETERS:
ntot          = size(data.NNs,1);
adInd         = 3;
lx            = ntot-4;
NNbar         = NNvec(:,1);
sumWorkingAge = sum(NNbar([1:lx,lx+3]));

nc = max(struct2array(data.compindex));
    
zn = zeros(ntot,1);

t0 = data.tvec(1);
y0 = [S0;repmat(zn,6,1);NNbar-S0;repmat(zn,nc-9,1);S0];

Tout       = t0;
Iout       = zn';
Isaout     = zn';
Issout     = zn';
Insout     = zn';
Hout       = zn';
Dout       = zn';
Wout       = [];
hwout      = [];
poutout    = 0;
betamodout = 1;
Vout       = zn';
rout       = 0;

%% LOOP:

i = 1;

tend = data.tvec(end);

while Tout(end)<tend; 

    Wit               = WitMat(:,i);    
    NNfeed            = NNvec(:,i);
    NNfeed(NNfeed==0) = 1;
    D                 = Dvec(:,:,i);

    %Vaccination Rollout by Sector
    NNnext              = NNvec(:,i);
    NNnext(lx+[1,2])    = 1;
    NNnext([1:lx,lx+3]) = NNnext([1:lx,lx+3])/sumWorkingAge;
    NNnext(end)         = 1;
    
%     p2.ratep1 = NNnext.*[repmat(p2.aratep1(3),lx,1);p2.aratep1];    
%     p2.ratep2 = NNnext.*[repmat(p2.aratep2(3),lx,1);p2.aratep2];
%     p2.ratep3 = NNnext.*[repmat(p2.aratep3(3),lx,1);p2.aratep3];
%     p2.ratep4 = NNnext.*[repmat(p2.aratep4(3),lx,1);p2.aratep4];

    [tout,Iclass,Isaclass,Issclass,Insclass,Hclass,Dclass,pout,betamod,Vclass,y0,inext]=...
     integr8(data,NNfeed,D,i,t0,tend,dis,y0,inp3,p2);
    
    Tout       = [Tout;tout(2:end)];  
    Iout       = [Iout;Iclass(2:end,:)];
    Isaout     = [Isaout;Isaclass(2:end,:)];
    Issout     = [Issout;Issclass(2:end,:)];
    Insout     = [Insout;Insclass(2:end,:)];
    Hout       = [Hout;Hclass(2:end,:)];
    Dout       = [Dout;Dclass(2:end,:)]; W   = Wit'.*ones(length(tout),lx);
    Wout       = [Wout;W(1:end-1,:)];    hw  = data.hw(i,:).*ones(length(tout),lx);
    hwout      = [hwout;hw(1:end-1,:)];
    poutout    = [poutout;pout(2:end)];
    betamodout = [betamodout;betamod(2:end)];
    Vout       = [Vout;Vclass(2:end,:)];

%         if hospInc==0
%             Rt(j)=heComputeEigs(pr,beta,D,NNfeed,ntot,Sclass(end,:)');
%         end

    if Tout(end)<tend;
        
        data.tvec = [data.tvec(1:end-1),Tout(end),tend];

        t0 = Tout(end);

        Xh2w                   = NNvec(1:lx,inext)-NNvec(1:lx,i);%Addition to each wp next intervention step
        Xw2h                   = -Xh2w; 
        Xw2h(Xw2h<0)           = 0;
        Xw2h                   = Xw2h./NNvec(1:lx,i);
        Xw2h(NNvec(1:lx,i)==0) = 0;

        if NNvec(lx+adInd,i)>0;%when would this not be the case?
            Xh2w(Xh2w<0) = 0;
            Xh2w         = Xh2w/NNvec(lx+adInd,i);
        else
            Xh2w         = 0;
        end

        %Move all infection statuses:

        y0    = reshape(y0,[ntot,nc]);%IC
        y0w2h = y0(1:lx,:).*repmat(Xw2h,1,nc);%IC%number of people to be put at home (+)
        y0w2h = [-y0w2h;sum(y0w2h,1)];

        y0h2w = y0(lx+adInd,:);
        y0h2w = kron(y0h2w,Xh2w);
        y0h2w = [y0h2w;-sum(y0h2w,1)];
        
        y0([1:lx,lx+adInd],:) = y0([1:lx,lx+adInd],:)+y0w2h+y0h2w;
        y0                    = reshape(y0,ntot*nc,1);
        i                     = inext;

    end   

end
    
%% OUTPUTS:  

Wout  = [Wout;Wout(end,:)];
hwout = [hwout;hwout(end,:)];
g     = [Tout,Wout,hwout,Isaout,Issout,Insout,Hout,Dout,Vout,betamodout];
f     = [Tout,...
         sum(Iout,2),...
         sum(Hout,2),...
         sum(Dout,2),...
         poutout,...
         betamodout,...  
         sum(Vout(:,lx+1),2),...
         sum(Vout(:,lx+2),2),...
         sum(Vout(:,[1:lx,lx+3]),2),...
         sum(Vout(:,lx+4),2),...
         sum(Dout(:,lx+1),2),...
         sum(Dout(:,lx+2),2),...
         sum(Dout(:,[1:lx,lx+3]),2),...
         sum(Dout(:,lx+4),2)];
  
end

%%

function [tout,Iclass,Isaclass,Issclass,Insclass,Hclass,Dclass,pout,betamod,Vclass,y0new,inext]=...
          integr8(data,NN0,D,i,t0,tend,dis,y0,inp3,p2)
%% CALL:

ntot = size(data.NNs,1);
fun  = @(t,y)ODEs(data,NN0,D,i,t,dis,y,p2);

if strcmp(inp3,'Elimination');
	options = odeset('Events',@(t,y)elimination(t,y,data,NN0,D,ntot,dis,i,p2));
elseif strcmp(inp3,'Economic Closures');
    options = odeset('Events',@(t,y)reactive_closures(t,y,data,NN0,D,ntot,dis,i,p2));
elseif strcmp(inp3,'School Closures');
    options = odeset('Events',@(t,y)reactive_closures(t,y,data,NN0,D,ntot,dis,i,p2));
elseif strcmp(inp3,'No Closures');
    options = odeset('Events',@(t,y)unmitigated(t,y,data,NN0,D,ntot,dis,i,p2));
else
    error('Unknown Mitigation Strategy!');
end

% try 
    [tout,yout,~,~,ie] = ode45(fun,[t0 tend],y0,options);
% catch
%     options.RelTol     = 1e-5;
%     options.AbsTol     = 1e-8;
%     [tout,yout,~,~,ie] = ode45(fun,[t0 tend],y0,options);
% end

y0new     = yout(end,:)'; 
if tout(end)<tend;
    inext = data.inext(ie(end));
else
    inext = NaN;
end

%% OUTPUT VARIABLES:

indices = 1:ntot;
compindex = data.compindex;
S     = yout(:,(compindex.S_index(1)-1)*ntot + indices);
Sn    = yout(:,(compindex.S_index(2)-1)*ntot + indices);
Ina   = yout(:,(compindex.I_index(1)-1)*ntot + indices);
Isa   = yout(:,(compindex.I_index(2)-1)*ntot + indices);
Ins   = yout(:,(compindex.I_index(3)-1)*ntot + indices);
Iss   = yout(:,(compindex.I_index(4)-1)*ntot + indices);
H     = yout(:,(compindex.H_index(1)-1)*ntot + indices);
D     = yout(:,(compindex.D_index(1)-1)*ntot + indices);


Iclass   = Ina + Isa + Ins + Iss; %yout(:, 2*ntot+1: 3*ntot) + yout(:, 3*ntot+1: 4*ntot) + ...
%            yout(:, 4*ntot+1: 5*ntot) + yout(:, 5*ntot+1: 6*ntot) + ...
%            yout(:,11*ntot+1:12*ntot) + yout(:,12*ntot+1:13*ntot) + ...
%            yout(:,13*ntot+1:14*ntot) + yout(:,14*ntot+1:15*ntot);
Isaclass = Isa; %yout(:, 3*ntot+1: 4*ntot) + yout(:,12*ntot+1:13*ntot);
Issclass = Iss; %yout(:, 5*ntot+1: 6*ntot) + yout(:,14*ntot+1:15*ntot);
Insclass = Ins; %yout(:, 4*ntot+1: 5*ntot) + yout(:,13*ntot+1:14*ntot); 
Hclass   = H; %yout(:, 6*ntot+1: 7*ntot) + yout(:,15*ntot+1:16*ntot);
Dclass   = D;%yout(:,17*ntot+1:18*ntot);
Vclass   = zeros(size(D));%yout(:,18*ntot+1:19*ntot);

%% TIME-DEPENDENT PARAMETERS:

occ   = max(1,sum(Hclass,2));
Hmax  = p2.Hmax;
SHmax = p2.SHmax;
th0   = max(1,1+1.87*((occ-Hmax)/(SHmax-Hmax)));

pd  = min(th0.*dis.pd',1);
Th  = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
mu  = pd./Th;
ddk = 10^5*sum(mu.*Hclass,2)/sum(NN0);

sd_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);

if i==1;
    betamod = ones(size(occ));
elseif any(i==data.imand);
    betamod = min(max(p2.sdl,sd_fun(p2.sdl,p2.sdb,ddk)), max(p2.sdl,sd_fun(p2.sdl,p2.sdb,2)));
else
    betamod = max(p2.sdl,sd_fun(p2.sdl,p2.sdb,ddk));
end


amp   = (Sn+(1-dis.heff).*(S-Sn))./S;
ph    = amp.*dis.ph';
Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g3    = (1-pd)./Th;
h     = ph./Ts;
dur   = p2.dur;
qh    = ph./(Ts-dur);

Hdot   = h.*Ins      +qh.*Iss        -(g3+mu).*H;
occdot = sum(Hdot,2);%+Hv1dot
r      = occdot./occ;


Ip    = 10^5*sum(Ina+Ins+Isa+Iss,2)/sum(NN0);%+Inav1+Insv1+Isav1+Issv1
trate = p2.trate;
b0    = 2.197;
b1    = 0.1838;
b2    = -1.024;

if i~=5;
    pout = (Ip<trate) .*   (1./(1+exp(b0+b1*Ip+b2*log10(trate))))/dur + ...
           (Ip>=trate).*min(1./(1+exp(b0+b1*Ip+b2*log10(trate))),trate/10^5)/dur;
    pout = pout.*(tout>p2.t_tit).*(tout<p2.end);    
else
    pout = zeros(size(tout));
end

end

%%

function [f,g]=ODEs(data,NN0,D,i,t,dis,y,p2)

ntot=size(data.NNs,1);

y_mat = reshape(y,ntot,[]); 


%% IC:
compindex = data.compindex;
% S=      y(0*ntot+1:1*ntot);
% E=      y(1*ntot+1:2*ntot);
% Ina=    y(2*ntot+1:3*ntot);
% Isa=    y(3*ntot+1:4*ntot);
% Ins=    y(4*ntot+1:5*ntot);
% Iss=    y(5*ntot+1:6*ntot);
% H=      y(6*ntot+1:7*ntot);
% R=      y(7*ntot+1:8*ntot);
% Sn=     y(19*ntot+1:20*ntot);

S=      y_mat(:,compindex.S_index(1));
E=      y_mat(:,compindex.E_index(1));
Ina=    y_mat(:,compindex.I_index(1));
Isa=    y_mat(:,compindex.I_index(2));
Ins=    y_mat(:,compindex.I_index(3));
Iss=    y_mat(:,compindex.I_index(4));
H=      y_mat(:,compindex.H_index(1));
R=      y_mat(:,compindex.R_index(1));
Sn=     y_mat(:,compindex.S_index(2));

%% HOSPITAL OCCUPANCY:

occ   = max(1,sum(H)); %+Hv1
Hmax  = p2.Hmax;
SHmax = p2.SHmax;

%% TIME-DEPENDENT DISEASE PARAMETERS:

%Amplitudes
amp = (Sn+(1-dis.heff).*(S-Sn))./S;
th0 = max(1,1+1.87*((occ-Hmax)/(SHmax-Hmax)));

%Probabilities
ph = amp.*dis.ph;
pd = min(th0*dis.pd,1);

%Calculations
Ts = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
Th = ((1-pd).*dis.Threc)+(pd.*dis.Thd);

sig1 = dis.sig1;
sig2 = dis.sig2;
g1   = dis.g1;
g2   = (1-ph)./Ts;
g3   = (1-pd)./Th;
h    = ph./Ts;
mu   = pd./Th;
nu   = dis.nu;

%Transmission
red  = dis.red;
beta = dis.beta;

%Preparedness
dur    = p2.dur;
qg1    = p2.qg1;
qg2    = (1-ph)./(Ts-dur);
qh     = ph./(Ts-dur);


%% FOI:

phi = 1 .* dis.rr_infection;  %+data.amp*cos((t-32-data.phi)/(365/2*pi));

ddk    = 10^5*sum(mu.*(H))/sum(NN0);%+Hv1
sd_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);

if i==1;
    betamod = 1;
elseif any(i==data.imand);
    betamod = min(max(p2.sdl,sd_fun(p2.sdl,p2.sdb,ddk)), max(p2.sdl,sd_fun(p2.sdl,p2.sdb,2)));
else
    betamod = max(p2.sdl,sd_fun(p2.sdl,p2.sdb,ddk));
end

I       = (red*Ina+Ins); %+(1-trv1)*(red*Inav1+Insv1);%Only non-self-isolating compartments
foi     = phi.*beta.*betamod.*(D*(I./NN0));

seedvec = 10^-15*sum(data.Npop)*ones(ntot,1);
seed    = phi.*beta.*betamod.*(D*(seedvec./NN0));

%% SELF-ISOLATION:
    
if t<p2.end && i~=5 && t>=p2.t_tit;
    Ip    = 10^5*sum(Ina+Ins+Isa+Iss)/sum(NN0);%+Inav1+Insv1+Isav1+Issv1
    trate = p2.trate;
    b0    = 2.197;
    b1    = 0.1838;
    b2    = -1.024;
    p3    = p2.self_isolation_compliance .* ((Ip<trate) .*   (1./(1+exp(b0+b1*Ip+b2*log10(trate))))/dur + ...
            (Ip>=trate).*min(1./(1+exp(b0+b1*Ip+b2*log10(trate))),trate/10^5)/dur);
    p4    = p3;
else       
    p3=0;
    p4=0;
        
end

%% VACCINATION:

% %Uptake
% uptake=uptake-V./NNage;
% uptake(uptake<0)=0;
% uptake(uptake>0)=1;
% uptake=[repmat(uptake(3),numSectors,1);uptake];

% nonVax=NN0-V;
%S (and R) ./nonVax accounts for the (inefficient) administration of vaccines to exposed, infectious and hospitalised people
%nonVax approximates S+E+I+H+R (unvaccinated) and D (partially)
%S (or R) ./nonVax is approximately 1 (or 0) when prevalence is low but is closer to 0 (or 1) when prevalence is high
%nonVax is non-zero as long as uptake is less than 100%

% if t>=pend
%     v1rates=zeros(ntot,1);
%     v1rater=zeros(ntot,1);
%     Vdot=   zeros(ntot,1);  
% elseif t>=p2.startp5
%     v1rates=p2.ratep5.*S./nonVax;
%     v1rater=p2.ratep5.*R./nonVax;
%     Vdot=   p2.ratep5;
      
% elseif t>=startp4
%     v1rates=ratep4.*S./nonVax;
%     v1rater=ratep4.*R./nonVax;
%     Vdot=   ratep4;
% elseif t>=startp3
%     v1rates=ratep3.*S./nonVax;
%     v1rater=ratep3.*R./nonVax;
%     Vdot=   ratep3;
% elseif t>=startp2
%     v1rates=ratep2.*S./nonVax;
%     v1rater=ratep2.*R./nonVax;
%     Vdot=   ratep2;
% elseif t>=startp1
%     v1rates=ratep1.*S./nonVax;
%     v1rater=ratep1.*R./nonVax;
%     Vdot=   ratep1;
% else
%     v1rates=zeros(ntot,1);
%     v1rater=zeros(ntot,1);
%     Vdot=   zeros(ntot,1);
% end

Sndot=      -Sn.*(foi+seed) ; %         -v1rates.*Sn./S;

%% EQUATIONS:

Sdot=       -S.*(foi+seed)  +nu.*R ; %  -v1rates                +nuv1.*Sv1;
% Shv1dot=    v1rates     -hrv1*Shv1   -Shv1.*foi;
% Sv1dot=                  hrv1*Shv1   -Sv1.*(1-scv1).*foi  -nuv1.*Sv1;  %+nu.*Rv1; 

Edot=        S.*(foi+seed)   -(sig1+sig2).*E;% +Shv1.*foi
% Ev1dot=                                  Sv1.*(1-scv1).*foi  -(sig1+sig2).*Ev1;

Inadot=     sig1.*E     -g1.*Ina                -p3*Ina;
Insdot=     sig2.*E     -(g2+h).*Ins            -p4*Ins;
% Inav1dot=   sig1.*Ev1   -g1.*Inav1              -p3*Inav1;
% Insv1dot=   sig2.*Ev1   -(g2_v1+h_v1).*Insv1    -p4*Insv1;

Isadot=     p3*Ina      -qg1.*Isa;
Issdot=     p4*Ins      -(qg2+qh).*Iss;
% Isav1dot=   p3*Inav1    -qg1.*Isav1;
% Issv1dot=   p4*Insv1    -(qg2_v1+qh_v1).*Issv1;

Hdot=       h.*Ins      +qh.*Iss        -(g3+mu).*H;
% Hv1dot=     h_v1.*Insv1 +qh_v1.*Issv1   -(g3+mu).*Hv1;

Rdot=       g1.*Ina     +qg1.*Isa   +g2.*Ins        +qg2.*Iss       +g3.*H      -nu.*R  ;%    -v1rater;
% Rv1dot=     g1.*Inav1   +qg1*Isav1  +g2_v1.*Insv1   +qg2_v1.*Issv1  +g3.*Hv1    +v1rater;   %-nu.*Rv1;

DEdot=      mu.*H    ; %           +mu.*Hv1;     

% Shv1dot = zeros(size(Rdot));
% Sv1dot = zeros(ntot,1);
% Ev1dot = zeros(ntot,1);
% Inav1dot = zeros(ntot,1);
% Insv1dot = zeros(ntot,1);
% Isav1dot = zeros(ntot,1);
% Issv1dot = zeros(ntot,1);
% Hv1dot = zeros(ntot,1);
% Rv1dot = zeros(ntot,1);
% Vdot = zeros(ntot,1);

% Hdot=       h.*(Ins+Iss)    -(g3+mu).*(min(occ,SHmax)/occ).*H   -(g3_oc+mu_oc).*(max(0,occ-SHmax)/occ).*H;
% Hv1dot=     h_v1.*(Insv1+Issv1)  -(g3+mu).*(min(occ,SHmax)/occ).*Hv1   -(g3_oc+mu_oc).*(max(0,occ-SHmax)/occ).*Hv1;
% 
% Rdot= g1.*(Ina+Isa)+g2.*(Ins+Iss) + g3.*(min(occ,SHmax)/occ).*H  + g3_oc.*(max(0,occ-SHmax)/occ).*H    -nu.*R   -v1rater;
% Rv1dot= g1.*(Inav1+Isav1)+g2_v1.*(Insv1+Issv1)+g3.*(min(occ,SHmax)/occ).*Hv1 +g3_oc.*(max(0,occ-SHmax)/occ).*Hv1+v1rater;%-nu.*Rv1;
% 
% DEdot= mu.*(min(occ,SHmax)/occ).*H+mu_oc .*(max(0,occ-SHmax)/occ).*H+mu.*(min(occ,SHmax)/occ).*Hv1  +mu_oc.*(max(0,occ-SHmax)/occ).*Hv1;     

%% OUTPUT:

f= [Sdot;Edot;...
    Inadot;Isadot;Insdot;Issdot;...
    Hdot;Rdot;...
%     Shv1dot;Sv1dot;Ev1dot;...
%     Inav1dot;Isav1dot;Insv1dot;Issv1dot;...
%     Hv1dot;Rv1dot;...
    DEdot;...%Vdot;
    Sndot];
% f= [Sdot;Edot;...
%     Iadot;Ipdot;...
%     Inmdot;Ismdot;...
%     Insdot;Issdot;...
%     Qmdot;Qsdot;...
%     Hdot;DEdot;Rdot];
f(y<eps) = max(0,f(y<eps)); 

g=h.*(Ins+Iss); %+h_v1.*(Insv1+Issv1);%Hin

end

%%

function [value,isterminal,direction] = elimination(t,y,data,N,D,ntot,dis,i,p2)
    
    ymat = reshape(y,ntot,[]);
    compindex = data.compindex;
    
    S    = ymat(:,compindex.S_index(1));
    H    = ymat(:,compindex.H_index(1));
    Sn   = ymat(:,compindex.S_index(2));
    
    Ina   = ymat(:,compindex.I_index(1));
    Isa   = ymat(:,compindex.I_index(2));
    Ins   = ymat(:,compindex.I_index(3));
    Iss   = ymat(:,compindex.I_index(4));
    occ   = max(1,sum(H)); %+Hv1
    
    amp   = (Sn+(1-dis.heff).*(S-Sn))./S;
    ph    = amp.*dis.ph;
    Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    g2    = (1-ph)./Ts;
    h     = ph./Ts;
     dur   = p2.dur;
    
    if t<p2.t_tit;   
        p3 = 0;
        p4 = 0;
    else
        Ip    = 10^5*sum(Ina+Ins+Isa+Iss)/sum(N); %+Inav1+Insv1+Isav1+Issv1
        trate = p2.trate;
        b0    = 2.197;
        b1    = 0.1838;
        b2    = -1.024;
        p3    = p2.self_isolation_compliance .* ((Ip<trate) .*   (1./(1+exp(b0+b1*Ip+b2*log10(trate))))/dur + ...
                (Ip>=trate).*min(1./(1+exp(b0+b1*Ip+b2*log10(trate))),trate/10^5)/dur);
        p4    = p3;
    end
    
    Rt1 = get_R(ntot,dis,h,g2,S,data.NNvec(:,3),data.Dvec(:,:,3),dis.beta,1,p3,p4);
    Rt2 = get_R(ntot,dis,h,g2,S,data.NNvec(:,5),data.Dvec(:,:,5),dis.beta,1,0,0);
    
    %% Event 1: Early Lockdown
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0);
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: Late Lockdown
    
    value(2)      = - abs(i-1) + min(occ-0.95*p2.Hmax,0);
    direction(2)  = 1;
    isterminal(2) = 1;
    
    %% Event 3: Reopening
    
    value(3)      = - abs(i-2) + min(t-(data.tvec(end-1)+7),0) + min(1.0000-Rt1,0);
    direction(3)  = 0;
    isterminal(3) = 1;
    
    %% Event 4: Relockdown
    
    value(4)      = - abs(i-3) + min(t-(data.tvec(end-1)+0.1),0) + min(Rt1-1.2000,0);
    direction(4)  = 0;
    isterminal(4) = 1;
    
    %% Event 5: End
    
    value(5)      = floor(i/5) + abs(min(t-(data.tvec(end-1)+0.1),0)) + min(t-p2.end,0)*min(1.00-Rt2,0);
    %measures can be removed at any stage if (Rt<1) or (after end of vaccination campaign)
    direction(5)  = 0;
    isterminal(5) = 1;
    
end 

function [value,isterminal,direction] = reactive_closures(t,y,data,N,D,ntot,dis,i,p2)
    
    ymat = reshape(y,ntot,[]);
    compindex = data.compindex;

    
    S    = ymat(:,compindex.S_index(1));
    H    = ymat(:,compindex.H_index(1));
    Sn   = ymat(:,compindex.S_index(2));
    
    Ins   = ymat(:,compindex.I_index(3));
    Iss   = ymat(:,compindex.I_index(4));
    occ   = max(1,sum(H)); %+Hv1
    
    Hmax  = p2.Hmax;
    SHmax = p2.SHmax;
    amp   = (Sn+(1-dis.heff).*(S-Sn))./S;
    th0   = max(1,1+1.87*((occ-Hmax)/(SHmax-Hmax)));
    ph    = amp.*dis.ph;
    pd    = min(th0*dis.pd,1);
    Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    Th    = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
    g2    = (1-ph)./Ts;
    g3    = (1-pd)./Th;
    h     = ph./Ts;
    mu    = pd./Th;
%     h_v1  = dis.h_v1;
    dur   = p2.dur;
    qh    = ph./(Ts-dur);


    Hdot   = h.*Ins      +qh.*Iss        -(g3+mu).*H;
    occdot = sum(Hdot);%+Hv1dot
    r      = occdot/occ;
    Tcap   = t + log(p2.Hmax/occ)/r;
    Tcap   = Tcap-4;
    Rt2    = get_R(ntot,dis,h,g2,S,data.NNvec(:,5),data.Dvec(:,:,5),dis.beta,1,0,0);
    
    %% Event 1: Response Time
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0) ;
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: Early Lockdown
    
    value(2)     = - abs((i-2)*(i-4)) + min(t-(data.tvec(end-1)+0.1),0) + min(t-Tcap,0);
    direction(2) = 1;
    if r>0.025
        isterminal(2) = 1;
    else
        isterminal(2) = 0;
    end
    
    %% Event 3: Late Lockdown
    
    value(3)      = - abs((i-1)*(i-2)*(i-4)) + min(t-(data.tvec(end-1)+0.1),0) + min(occ-0.95*p2.Hmax,0);
    direction(3)  = 1;
    isterminal(3) = 1;
    
    %% Event 4: Reopening
    
    value(4)      = abs(i-3) + abs(min(t-(data.tvec(end-1)+7),0)) + max(0,occ-p2.thl);
    direction(4)  = -1;
    isterminal(4) = 1;
    
    %% Event 5: End
    
    value(5)      = abs((i-1)*(i-2)*(i-4)) + abs(min(t-(data.tvec(end-1)+0.1),0)) + (min(0.025-r,0)*max(0,occ-p2.thl) + min(t-p2.end,0))*min(1.00-Rt2,0);
    %measures can be removed if (not in hard lockdown) and ((Rt<1) or (after end of vaccination campaign and below 25% occupancy or low growth rate))
    direction(5)  = 0;
    isterminal(5) = 1;
    
end

function [value,isterminal,direction] = unmitigated(t,y,data,N,D,ntot,dis,i,p2)
    

    ymat = reshape(y,ntot,[]);
    compindex = data.compindex;

    S    = ymat(:,compindex.S_index(1));
    H    = ymat(:,compindex.H_index(1));
    Sn   = ymat(:,compindex.S_index(2));
    occ  = max(1,sum(H)); %+Hv1
        
    amp  = (Sn+(1-dis.heff).*(S-Sn))./S;
    ph   = amp.*dis.ph;
    Ts   = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    g2   = (1-ph)./Ts;
    h    = ph./Ts;
    
    Rt2  = get_R(ntot,dis,h,g2,S,data.NNvec(:,5),data.Dvec(:,:,5),dis.beta,1,0,0);

    %% Event 1: Response Time
    
    value(1)      = - abs(i-1) + min(t-p2.Tres,0);
    direction(1)  = 1;
    isterminal(1) = 1;
    
    %% Event 2: Late Lockdown
    
    value(2)      = - abs(i-1) + min(occ-0.95*p2.Hmax,0);
    direction(2)  = 1;
    isterminal(2) = 1;

    %% Event 3: End
    
    value(3)      = floor(i/5) + abs(min(t-(data.tvec(end-1)+0.1),0)) + min(t-p2.end,0)*min(1.00-Rt2,0);
    %measures can be removed at any stage if (Rt<1) or (after end of vaccination campaign)
    direction(3)  = 0;
    isterminal(3) = 1;
    
end

% function Rt = rep_num(ntot,dis,h,g2,S,N,D,betamod,p3,p4)
%         
%     FOIu = repmat(S,1,ntot).*dis.beta.*dis.rr_infection.*betamod.*D./repmat(N',ntot,1);%+Shv1
%     
%     F                                = zeros(3*ntot,3*ntot);
%     F(1:ntot, 1*ntot+1:3*ntot) = [dis.red*FOIu,  FOIu]; %,  dis.red*(1-dis.trv1)*FOIu,  (1-dis.trv1)*FOIu
%     
%     onesn                            = ones(ntot,1);
%     vvec                             = [(dis.sig1+dis.sig2).*onesn;  (dis.g1+p3).*onesn;...
%                                         (g2+h+p4).*onesn;  ];  %    (dis.sig1+dis.sig2).*onesn;      (dis.g1+p3).*onesn;         (dis.g2_v1+dis.h_v1+p4).*onesn
%     V                                = diag(vvec);
%     V(1*ntot+1:2*ntot,1:ntot)        = diag(-dis.sig1.*onesn);
%     V(2*ntot+1:3*ntot,1:ntot)        = diag(-dis.sig2.*onesn);
%     
%     NGM = F/V;
%     Rt  = eigs(NGM,1);
%     
% end

