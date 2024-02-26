function [data,dis,p2] = p2Params(data,dis,vaccine_day,bpsv)

%% PREPAREDNESS PARAMETERS:

p2 = struct;

p2.self_isolation_compliance = data.self_isolation_compliance;
% p2.Tres  = data.Tres;                       %Response Time
% p2.t_tit = data.t_tit;                      %Test-Isolate-Trace Time
p2.trate = data.trate;                      %Test-Isolate-Trace Rate
p2.sdl   = data.sdl;                        %Social Distancing Asymptote
p2.sdb   = data.sdb;                        %Social Distancing Steepness
p2.Hmax  = data.Hmax*sum(data.Npop)/10^5;   %Hospital Capacity
t_vax    = data.t_vax;                      %Vaccine Administration Time

if bpsv==0 && vaccine_day==365
    vaccination_rate_pc = 0.5/100;
    vaccine_uptake  = .4;
else
    vaccination_rate_pc    = data.vaccination_rate_pc;  %Vaccine Administration Rate
    vaccine_uptake  = data.vaccine_uptake;                    %Vaccine Uptake
end


% p2.Tres = data.tvec(1) + (-data.tvec(1) + data.Tres)*dis.Td/data.Td_CWT;

%Test-Isolate-Trace
% p2.t_tit  = data.tvec(1) + (-data.tvec(1) + data.t_tit)*dis.Td/data.Td_CWT;

p2.Tres = data.rts;
p2.t_tit = data.rts;
p2.dur    = 1;
p2.frac_sym_infectiousness_remaining = min(1,p2.dur./(dis.Tsr+dis.Tsh)*2);
p2.frac_asym_infectiousness_remaining = min(1,p2.dur./dis.Tay);

% p2.qg1    = 1/(dis.Tay-p2.dur);
% p2.qg2    = (1-dis.ph)./(dis.Ts-p2.dur);
% p2.qg2_v1 = (1-(1-dis.hv1)*dis.ph)./(dis.Ts_v1-p2.dur);
% p2.qg2_v2 = (1-(1-dis.hv2)*dis.ph)./(dis.Ts_v2-p2.dur);
% p2.qh     = dis.ph./(dis.Ts-p2.dur);
% p2.qh_v1  = (1-dis.hv1)*dis.ph./(dis.Ts_v1-p2.dur);
% p2.qh_v2  = (1-dis.hv2)*dis.ph./(dis.Ts_v2-p2.dur);

%Hospital Capacity
p2.thl   = max(1,0.25*p2.Hmax);%lower threshold can't be less than 1 occupant
p2.Hmax  = max(4*p2.thl,p2.Hmax);
p2.SHmax = 2*p2.Hmax;

%% Vaccine Uptake
Npop    = data.Npop;
NNage   = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
over14frac = Npop(4)/sum(Npop(2:4));
% vaccine_uptake = min(0.99*(1-NNage(1)/sum(NNage)),vaccine_uptake); %population uptake cannot be greater than full coverage in non-pre-school age groups
% up3fun  = @(u3) vaccine_uptake*sum(NNage) - u3*(NNage(2)/2 + NNage(3)) - min(1.5*u3,1)*NNage(4);
% if up3fun(0)*up3fun(1)<=0
%     u3  = fzero(up3fun,[0 1]);
% else
%     u3  = fminbnd(up3fun,0,1);
% end
% u4      = min(1.5*u3,1);
% u1      = 0;
% up2fun  = @(u2) u2*NNage(2) + u3*NNage(3) + u4*NNage(4) - vaccine_uptake*sum(NNage);
% u2      = fzero(up2fun,[0 1]);


p2.t_vax2 = p2.Tres + 7 + vaccine_day;  

if bpsv == 1
    t_vax = p2.Tres;
end

vaccination_rate = vaccination_rate_pc*sum(data.Npop);

t_vax2 = p2.t_vax2;
t_bspv = t_vax2 - t_vax;
% if vaccine_uptake*sum(NNage(2:4)) < t_bspv*arate
%     arate = sum(vaccine_uptake.*NNage)/t_bspv;
% end
uptake  = vaccine_uptake*[0 over14frac 1 1]; %[u1,u2,u3,u4];
% if abs((uptake*NNage'/sum(NNage))-vaccine_uptake)>1e-10
%     error('Vaccine uptake error!');
% end

%Vaccine Administration Rate
t_ages     = min((uptake.*NNage)/vaccination_rate,Inf);%vaccination_rate may be 0

if bpsv==1
    % if primer starts before booster, there are t_bspv days of primer.
    % these people must then be boosted.
    if t_bspv < t_ages(4)
        t_ages = [t_ages(4)-t_bspv, t_ages(4), t_ages(3), t_ages(2)];
        p2.group_order = [4,4,3,2];
    else
        t_ages = [t_ages(4), t_bspv-t_ages(4), t_ages(4), t_ages(3), t_ages(2)];
        p2.group_order = [4,0,4,3,2];
    end
else
    % if booster starts before primer, just work down the ages
    t_ages     = [t_ages(4),t_ages(3),t_ages(2)];
    p2.group_order = [4,3,2];
end    
tpoints    = cumsum([min(t_vax, p2.t_vax2), t_ages]);
p2.tpoints = tpoints;
p2.arate = vaccination_rate;


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