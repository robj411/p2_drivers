close all;
clear all;

%% ons estimates
onsprev = readtable('onsprevalence.csv');
starts = arrayfun(@(i)sprintf('%d/%d/%d',onsprev.Var1(i),month(onsprev.Var2(i),"mmmm"),onsprev.Var3(i)),1:size(onsprev,1),'UniformOutput',false);
starts = starts(:);
ends = arrayfun(@(i)sprintf('%d/%d/%d',onsprev.Var5(i),month(onsprev.Var6(i),"mmmm"),onsprev.Var7(i)),1:size(onsprev,1),'UniformOutput',false);
ends = ends(:);
startdays = caldays(arrayfun(@(i) between(datetime('01/01/2020','Format','dd/M/yyyy'),datetime(starts(i),'Format','dd/M/yyyy'),'Days'),1:length(starts)));
enddays = caldays(arrayfun(@(i) between(datetime('01/01/2020','Format','dd/M/yyyy'),datetime(ends(i),'Format','dd/M/yyyy'),'Days'),1:length(ends)));
cases_per_day = onsprev.Var8 ./ (enddays' - startdays') .* 4;



%% infections

prevtable = readtable('10.UKprevalence.csv');
prevt    = prevtable.t;
prevI     = prevtable.I;
t1    = 1:244; 
all_cases_raw     = interp1(prevt,prevI,t1)'; 

scatter(t1,all_cases_raw);

%% compare to ons

prevIsums = arrayfun(@(i) sum(prevI(prevt<=enddays(i)&prevt>=startdays(i))),1:length(enddays));
plot(prevIsums(prevIsums>0)',cases_per_day(prevIsums>0))


%% tests

testtable = readtable('10.UKtesting.csv');
testt2    = flipud(testtable.date);
daily_tests_raw     = flipud(testtable.newVirusTestsByPublishDate);
dvec  = datevec(testt2);
dvec  = dvec(:,1)-2020;
t2    = day(testt2,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);

scatter(t2,daily_tests_raw);

%% infections identified

casetable = readtable('10.UKcases.csv');
caset3    = flipud(casetable.date);
positive_tests_raw     = flipud(casetable.newCasesBySpecimenDate);
dvec  = datevec(caset3);
dvec  = dvec(:,1)-2020;
t3    = day(caset3,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);

scatter(t3,positive_tests_raw);

%% proportion of infections identified

[newt,i1,i3] = intersect(t1,t3);
all_cases_tmp         = all_cases_raw(i1);
positive_tests_tmp         = positive_tests_raw(i3);
fraction_cases_identified_raw         = positive_tests_tmp./all_cases_tmp;

scatter(newt,fraction_cases_identified_raw);

%% collation

[t,i,i2] = intersect(newt,t2);
all_cases        = all_cases_tmp(i);
daily_tests        = daily_tests_raw(i2);
positive_tests        = positive_tests_tmp(i);
fraction_cases_identified        = fraction_cases_identified_raw(i);

[t,all_cases,positive_tests]
scatter(all_cases,positive_tests)
scatter(daily_tests,fraction_cases_identified)

Npop = 675.29;

%figure;
%scatter(I/Npop,p);
%figure;
%scatter(T/Npop,p);

%% logistic fit

%not log,not log
tit_fun1        = @(b0,b1,b2,x,y) 1./(1+exp(b0+b1*x+b2*y));
[tit_fit1,gof1] = fit([all_cases/Npop,daily_tests/Npop],fraction_cases_identified,tit_fun1,...
                      'StartPoint',[2.5,2.5,-2.5],'Lower',[log(9),0,-10],'Upper',[10,10,0],...
                      'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3);

%not log,log (original)
tit_fun2        = @(b0,b1,b2,x,y) 1./(1+exp(b0+b1*x+b2*log10(y)));
[tit_fit2,gof2] = fit([all_cases/Npop,daily_tests/Npop],fraction_cases_identified,tit_fun2,...
                      'StartPoint',[2.5,2.5,-2.5],'Lower',[log(9),0,-10],'Upper',[10,10,0],...
                      'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3);

%log,not log
tit_fun3        = @(b0,b1,b2,x,y) 1./(1+exp(b0+b1*log10(x)+b2*y));
[tit_fit3,gof3] = fit([all_cases/Npop,daily_tests/Npop],fraction_cases_identified,tit_fun3,...
                      'StartPoint',[2.5,2.5,-2.5],'Lower',[log(9),0,-100],'Upper',[100,100,0],...
                      'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3);

%log,log
tit_fun4        = @(b0,b1,b2,x,y) 1./(1+exp(b0+b1*log10(x)+b2*log10(y)));
[tit_fit4,gof4] = fit([all_cases/Npop,daily_tests/Npop],fraction_cases_identified,tit_fun4,...
                      'StartPoint',[2.5,2.5,-2.5],'Lower',[log(9),0,-10],'Upper',[10,10,0],...
                      'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3);

Ip      = linspace(0,1000,500);
Tp      = logspace(-2,3,500);
[Ip,Tp] = meshgrid(Ip,Tp);
                
figure;
hold on;
scatter3(all_cases/Npop,daily_tests/Npop,fraction_cases_identified,'r*');
surf(Ip,Tp,tit_fun2(tit_fit2.b0,tit_fit2.b1,tit_fit2.b2,Ip,Tp));
xlim([0 1000]);
ylim([0 1000]);
colormap(flipud(bone));
shading interp;
alpha 0.5;
view(-10,20);

% f  = figure('Units','centimeters','Position',[0 0 15 10]);
% set(f,'defaulttextInterpreter','latex');
% set(f,'defaultAxesTickLabelInterpreter','latex');
% set(f,'defaultLegendInterpreter','latex');
% set(f,'defaultColorbarTickLabelInterpreter','latex');
% set(f,'DefaultAxesFontSize',12);
% fs = 12;
% hold on;
% scatter3(I/Npop,T/Npop,p,'ko','filled');
% surf(Ip,Tp,tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip,Tp));
% %set(gca,'YScale','log');
% xlim([0 300]);
% ylim([0 300]);
% view(-20,30);
% %axis square;
% grid on;
% box on;
% %xticks([]);    
% %yticks([0.1,0.3,1,3,10]);
% %set(gca,'xticklabels',{'Year 1 - Jan','Jul',...
% %                       'Year 2 - Jan','Jul',...
% %                       'Year 3 - Jan','Jul'});       
% %set(gca,'yticklabels',{'0.1','0.3','1','3','10'});
% xlabel('Prevalence per 100k');
% ylabel('Daily Tests per 100k');
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',-45,'VerticalAlignment','middle')
% vec_pos = get(get(gca,'ylabel'),'Position');
% set(hYLabel,'position',vec_pos+[10 0 -0.12])
% %zlabel('Proportion');
% title('Proportion of Infections Identified');
% set(gca,'FontSize',fs);
% colormap(summer);
% shading interp;
% alpha 0.5;

figure;
hold on;
plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),4));
plot(4*ones(1,11),[0:0.1:1],'r');
plot(Ip(1,:),(4/10^5)*ones(size(Ip(1,:))),'r');
plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),40));
plot(40*ones(1,11),[0:0.1:1],'r');
plot(Ip(1,:),(40/10^5)*ones(size(Ip(1,:))),'r');
plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),130));
plot(130*ones(1,11),[0:0.1:1],'r');
plot(Ip(1,:),(130/10^5)*ones(size(Ip(1,:))),'r');
plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),340));
plot(340*ones(1,11),[0:0.1:1],'r');
plot(Ip(1,:),(340/10^5)*ones(size(Ip(1,:))),'r');
plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),840));
plot(840*ones(1,11),[0:0.1:1],'r');
plot(Ip(1,:),(840/10^5)*ones(size(Ip(1,:))),'r');
xlim([0 100]);

%% linear fit

% %figure;
% %scatter3(I/Npop,T/Npop,p,'r*');
% %hold on;
% Ip      = linspace(0,500,100);
% Tp      = linspace(0,1000,100);
% [Ip,Tp] = meshgrid(Ip,Tp);
% %surf(Ip,Tp,1./(1+10.^(1+0.2*Ip-0.02*Tp)));
% %colormap(flipud(bone));
% %shading interp;
% %alpha 0.5;
% 
% % IND = (p<=1)&(T/Npop<=1800)&(I/Npop<=200);
% % p   = p(IND);
% % T   = T(IND);
% % I   = I(IND);
% 
% tit_fun = @(b0,b1,b2,x,y) 1./(1+10.^(b0+b1*x+b2*y));
% tit_fit = fit([I/Npop,T/Npop],p,tit_fun,...
%           'StartPoint',[3,0.2,-0.02],'Lower',[2,0,-10],'Upper',[10,10,0],...
%           'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3)
% 
% figure;
% scatter3(I/Npop,T/Npop,p,'r*');
% hold on;
% surf(Ip,Tp,tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip,Tp));
% colormap(flipud(bone));
% shading interp;
% alpha 0.5;
% zlim([0 1]);
% 
% figure;
% hold on;
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),0));
% plot(0*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(0/10^5)*ones(size(Ip(1,:))),'r');
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),25));
% plot(25*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(25/10^5)*ones(size(Ip(1,:))),'r');
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),110));
% plot(110*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(110/10^5)*ones(size(Ip(1,:))),'r');
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),340));
% plot(330*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(340/10^5)*ones(size(Ip(1,:))),'r');

%% log10-scaled fits

% figure;
% scatter3(log10(I/(Npop)),log10(T/Npop),p,'r*')
% hold on;
% Ip=-1:0.1:5;
% Tp=-1:0.1:5;
% [Ip,Tp] = meshgrid(Ip,Tp);
% %surf(Ip,Tp,1./(1+10.^(-2*(-Ip+0.7))))
% %surf(Ip,Tp,1./(1+10.^(-3*(Tp-2.45))))
% surf(Ip,Tp,1./(1+10.^(6+2*Ip-3*Tp)))
% alpha 0.5;
% %xlim([0 500])
% %ylim([0 500])
% %view(0,0);
% colormap(flipud(bone));
% shading interp;
% 
% tit_fun = @(b0,b1,b2,x,y) 1./(1+exp(b0+b1*x+b2*y));
% tit_fit = fit([log10(I/(Npop)),log10(T/Npop)],p,tit_fun,...
%               'StartPoint',[5,2,-3],'Lower',[4.6,0,-10],'Upper',[10,10,0],...
%               'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3);
%                  
% figure;
% scatter3(log10(I/(Npop)),log10(T/Npop),p,'r*')
% hold on;
% surf(log10(Ip),log10(Tp),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,log10(Ip),log10(Tp)))
% alpha 0.5;
% %xlim([0 500])
% %ylim([0.01 500])
% %view(0,0);
% colormap(flipud(bone));
% shading interp;
% 
% figure;
% hold on;
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,log10(Ip(1,:)),log10(5)))
% plot(0*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(0/10^5)*ones(size(Ip(1,:))),'r');
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),25))
% plot(25*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(25/10^5)*ones(size(Ip(1,:))),'r');
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),110))
% plot(110*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(110/10^5)*ones(size(Ip(1,:))),'r');
% plot(Ip(1,:),tit_fun(tit_fit.b0,tit_fit.b1,tit_fit.b2,Ip(1,:),330))
% plot(330*ones(1,11),[0:0.1:1],'r');
% plot(Ip(1,:),(330/10^5)*ones(size(Ip(1,:))),'r');