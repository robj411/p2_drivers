function [cost,ccost_t] = p2Cost(data,dis,p2,g)

t  = g(:,1);
lx = length(data.obj);
ln = lx+4;

%% VLYL

deaths    = g(end,1+2*lx+4*ln+1:1+2*lx+5*ln);
cost(1,:) = deaths;

lyl       = deaths.*data.lgh;
cost(2,:) = lyl;

vlyl      = lyl*data.vly;
cost(3,:) = vlyl;

deaths          = g(:,1+2*lx+4*ln+1:1+2*lx+5*ln);
ccost_t(:,1:ln) = deaths.*data.lgh.*data.vly;

%% VSYL

Stu              = lx+2;
students         = data.NNs(Stu);
cost(4,lx+[1,2]) = students;

%Student Supply
isoasy       = g(:,1+2*lx+0*ln+Stu).*14/(dis.Tay-p2.dur);%.*(1-(1/3));%sfh still contribute; 14-day isolation period
isosym       = g(:,1+2*lx+1*ln+Stu);
isorec       = g(:,1+2*lx+1*ln+Stu).*(14-dis.Ts(Stu)+p2.dur)./(dis.Ts(Stu)-p2.dur);%.*(1-(1/3));
nissym       = g(:,1+2*lx+2*ln+Stu);
hospts       = g(:,1+2*lx+3*ln+Stu);
deaths       = g(:,1+2*lx+4*ln+Stu);
abs          = isoasy + isosym + isorec + nissym + hospts + deaths;%numbers of students
absint       = trapz(t,abs)/365;
cost(5,lx+1) = absint;
vsyl_sts     = absint*data.vsy;
cost(6,lx+1) = vsyl_sts;

%Student Demand
pres         = students-abs;
presl        = pres.*(1-g(:,1+data.EdInd));%.*(1-(1/3));%numbers of students
preslint     = trapz(t,presl)/365;%= (diff(t)'*presl)/365;
cost(5,lx+2) = preslint;
vsyl_std     = preslint*data.vsy;
cost(6,lx+2) = vsyl_std;

ccost_t(:,ln+lx+1) = cumtrapz(t,abs,1)./365.*data.vsy;
ccost_t(:,ln+lx+2) = cumtrapz(t,presl,1)./365.*data.vsy;

%% SGDPL

notEd = [1:(data.EdInd-1),(data.EdInd+1):lx];

%Labour Supply
hw            = g(:,1+1*lx+notEd);
isoasy        = g(:,1+2*lx+0*ln+notEd).*(1-hw).*14/(dis.Tay-p2.dur);%hw still contribute; 14-day isolation period
isosym        = g(:,1+2*lx+1*ln+notEd);
isorec        = g(:,1+2*lx+1*ln+notEd).*(1-hw).*(14-dis.Ts(notEd)'+p2.dur)./(dis.Ts(notEd)'-p2.dur);
nissym        = g(:,1+2*lx+2*ln+notEd);
hospts        = g(:,1+2*lx+3*ln+notEd);
deaths        = g(:,1+2*lx+4*ln+notEd);%number of workers absent
abspc         = max((isoasy + isosym + isorec + nissym + hospts + deaths)./data.NNs(notEd)',0);%percentage of workers absent
prespc        = 1-abspc;%percentage of workers present
presx         = prespc.^data.alp;%percentage of gdp output%the alpha relationship only holds for present workers!!!
absx          = 1-presx;%percentage of gdp lost
absxint       = trapz(t,absx);
gdpl_lbs      = absxint.*data.obj(notEd)';
cost(7,notEd) = gdpl_lbs;

%Labour Demand
w             = g(:,1+notEd);
x             = w.^data.alp;
xint          = diff(t)'*(1-x(1:end-1,:));
gdpl_lbd      = xint.*data.obj(notEd)';
cost(8,notEd) = gdpl_lbd;

%Consumer Demand
% betamod       = g(:,1+2*lx+6*ln+1);
% conloss       = min((1-betamod).*data.hconsl(notEd)',1).*data.hcon(notEd)';
% prdloss       = (absx+1-x).*data.obj(notEd)';
% difloss       = max(0,conloss-prdloss);%to avoid double counting
% gdpl_crd      = trapz(t,difloss);
cost(9,notEd) = 0;%gdpl_crd;
%zero during ld: betamod(sum(x,2)<44)=1;

%Medium-Term
cost(10,notEd) = 0;

ccost_t(:,2*ln+notEd) = cumtrapz(t,absx,1).*data.obj(notEd)';
ccost_t(:,3*ln+notEd) = cumtrapz(t,(1-x),1).*data.obj(notEd)';

%% IMPC

% tstart       = min(p2.Tres,data.tvec(2));%response time or late lockdown time
% w            = g(:,1+notEd);
% if sum(w(end,:))==44;    
%     tend     = data.tvec(end-1);%lifting time or simulation end time
% else
%     tend     = data.tvec(end);
%     error('Implementation Cost Error!');
% end
% betamod      = g(:,1+2*lx+6*ln+1);
% hospts       = min(sum(g(:,(1+2*lx+3*ln+1):(1+2*lx+4*ln)),2),p2.Hmax);
% vaxxed       = sum(g(:,(1+2*lx+5*ln+1):(1+2*lx+6*ln)),2);%number of people vaccinated
% units        = [max(0,tend-tstart),...
%                 data.trate*sum(data.Npop/10^5)*max(0,tend-p2.t_tit),...
%                 sum(data.Npop)*trapz(t,1-betamod),...
%                 trapz(t,hospts),...
%                 2*vaxxed(end)];
% impcost      = data.pppf*(data.impcost(:,1)' + units.*data.impcost(:,2)');
% cost(11,1:5) = impcost;
% 
% cunits                = [cumtrapz(t,(tstart<t)&(t<tend)),...
%                          data.trate*sum(data.Npop/10^5)*cumtrapz(t,(p2.t_tit<t)&(t<tend)),...
%                          sum(data.Npop)*cumtrapz(t,1-betamod),...
%                          cumtrapz(t,hospts),...
%                          2*vaxxed];
% ccost_t(:,4*ln+[1:5]) = data.pppf*(data.impcost(:,1)' + cunits.*data.impcost(:,2)');

end