f  = figure('Units','centimeters','Position',[0 0 25 10]);
set(f,'DefaultAxesFontSize',10);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
fs = 12;

% g = gca;
% g.Position=[0.0 0.28 0.93 0.625];
subplot(1,2,1);
hold on;
subplot(1,2,2);
hold on;
% subplot(2,2,3);
% hold on;
% subplot(2,2,4);
% hold on;

%% FLU2009

dis.ps = 0.669;
ihr    = dis.ps*[0.00697 0.00274 0.00274 0.00274 ...
                 0.00274 0.00561 0.00561 0.00561 ...
                 0.00561 0.00561 0.01060 0.01060 ...
                 0.01060 0.01546 0.01546 0.01546 0.01546];
ifr    = dis.ps*[0.000276 0.00011 0.00011 0.00012 ...
                 0.00012 0.00030 0.00030 0.00030 ...
                 0.00030 0.00065 0.00065 0.00065 ...
                 0.00065 0.00980 0.00980 0.00980 0.00980];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color','green');
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color','green');
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color','green');
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color','green');

%% FLU1957

dis.ps = 0.669;
ihr    = dis.ps*13.5*[0.0001 0.0001 0.0001 0.0001 ...
                 0.0001 0.0001 0.0001 0.0001 ...
                 0.0001 0.0025 0.0025 0.0025 ...
                 0.0025 0.0200 0.0200 0.0200 0.0200];    
ifr    = dis.ps*[0.0001 0.0001 0.0001 0.0001 ...
                 0.0001 0.0001 0.0001 0.0001 ...
                 0.0001 0.0025 0.0025 0.0025 ...
                 0.0025 0.0200 0.0200 0.0200 0.0200];   
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color',[0.93 0.69 0.13]);
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color',[0.93 0.69 0.13]);
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color',[0.93 0.69 0.13]);
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color',[0.93 0.69 0.13]);

%% FLU1918

dis.ps = 0.669;
ihr    = dis.ps*8*[0.02284 0.00398 0.00478 0.00983 ...
                 0.01700 0.02922 0.02470 0.02205 ...
                 0.01647 0.01195 0.01647 0.01169 ...
                 0.03081 0.04144 0.04941 0.04941 0.04941];
ifr    = dis.ps*[0.02284 0.00398 0.00478 0.00983 ...
                 0.01700 0.02922 0.02470 0.02205 ...
                 0.01647 0.01195 0.01647 0.01169 ...
                 0.03081 0.04144 0.04941 0.04941 0.04941];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color',[1 0.5 0]);
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color',[1 0.5 0]);
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color',[1 0.5 0]);
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color',[1 0.5 0]);

%% COVIDOM

dis.ps = 0.595;
ihr    = 1.85*[0.000016 0.000016 0.000408 0.000408 ...	
          0.010400 0.010400 0.034300 0.034300 ...	
          0.042500 0.042500 0.081600 0.081600 ...	
          0.118000 0.118000 0.166000 0.166000 0.184000].*...
         [1.10 1.10 0.78 0.78 ...	
          0.43 0.43 0.31 0.31 ...	
          0.20 0.20 0.14 0.14 ...	
          0.14 0.14 0.20 0.20 0.33]; 
ifr    = 1.85*[0.000016 0.000016 0.000070 0.000070 ...
          0.000309 0.000309 0.000844 0.000844 ...
          0.001610 0.001610 0.005950 0.005950 ...
          0.019300 0.019300 0.042800 0.042800 0.078000].*...
         [1.10 1.10 0.78 0.78 ...	
          0.43 0.43 0.31 0.31 ...	
          0.20 0.20 0.14 0.14 ...	
          0.14 0.14 0.20 0.20 0.33];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color','red');
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color','red');
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color','red');
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color','red');

%% COVIDWT

dis.ps = 0.595;
ihr    = [0.000016 0.000016 0.000408 0.000408 ...	
          0.010400 0.010400 0.034300 0.034300 ...	
          0.042500 0.042500 0.081600 0.081600 ...	
          0.118000 0.118000 0.166000 0.166000 0.184000];
ifr    = [0.000016 0.000016 0.000070 0.000070 ...
          0.000309 0.000309 0.000844 0.000844 ...
          0.001610 0.001610 0.005950 0.005950 ...
          0.019300 0.019300 0.042800 0.042800 0.078000];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color','magenta');
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color','magenta');
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color','magenta');
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color','magenta');

%% COVIDDE

dis.ps = 0.595;
ihr    = 1.85*[0.000016 0.000016 0.000408 0.000408 ...	
          0.010400 0.010400 0.034300 0.034300 ...	
          0.042500 0.042500 0.081600 0.081600 ...	
          0.118000 0.118000 0.166000 0.166000 0.184000];
ifr    = 1.85*[0.000016 0.000016 0.000070 0.000070 ...
          0.000309 0.000309 0.000844 0.000844 ...
          0.001610 0.001610 0.005950 0.005950 ...
          0.019300 0.019300 0.042800 0.042800 0.078000];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color','blue');
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color','blue');
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color','blue');
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color','blue');

%% SARS

dis.ps = 0.867;
ihr    = [0.0578 0.0578 0.0578 0.0578 ...	
          0.0816 0.0816 0.0816 0.0816 ...	
          0.3026 0.3026 0.3026 0.3026 ...	
          0.8670 0.8670 0.8670 0.8670 0.6018];
ifr    = dis.ps*[0.017 0.017 0.017 0.017 ...
                 0.024 0.024 0.024 0.024 ...
                 0.089 0.089 0.089 0.089 ...
                 0.255 0.255 0.255 0.255 0.177];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

subplot(1,2,1);
plot(100*ihr,'-','linewidth',1.5,'color','black');
subplot(1,2,2);
plot(100*ifr,'-','linewidth',1.5,'color','black');
% subplot(2,2,3);
% plot(100*phgs,'-','linewidth',1.5,'color','black');
% subplot(2,2,4);
% plot(100*pdgh,'-','linewidth',1.5,'color','black');

% opdgh = max(pdgh,0.60);
% spdgh = min(2*pdgh,opdgh);
% plot(100*opdgh,'--','linewidth',1.5,'color','black');
% plot(100*spdgh,'-.','linewidth',1.5,'color','black');

%%

subplot(1,2,1);
xlim([1 17]);
ylim([0 100]);
xticks(1:17);
xticklabels({'0-4','5-9','10-14','15-19',...
             '20-24','25-29','30-34','35-39',...
             '40-44','45-49','50-54','55-59',...
             '60-64','65-69','70-74','75-79','80+'});
xtickangle(45);
grid on;
box on;
ylabel('IHR (\%)');
set(gca,'FontSize',fs);
legend('Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Wildtype','Covid Delta','SARS',...
       'location','northwest');

subplot(1,2,2);
xlim([1 17]);
ylim([0 30]);
xticks(1:17);
xticklabels({'0-4','5-9','10-14','15-19',...
             '20-24','25-29','30-34','35-39',...
             '40-44','45-49','50-54','55-59',...
             '60-64','65-69','70-74','75-79','80+'});
xtickangle(45);
grid on;
box on;
ylabel('IFR (\%)');
set(gca,'FontSize',fs);
legend('Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Wildtype','Covid Delta','SARS',...
       'location','northwest');

% subplot(2,2,3);
% xlim([1 17]);
% ylim([0 100]);
% xticks(1:17);
% xticklabels({'0-4','5-9','10-14','15-19',...
%              '20-24','25-29','30-34','35-39',...
%              '40-44','45-49','50-54','55-59',...
%              '60-64','65-69','70-74','75-79','80+'});
% xtickangle(45);
% grid on;
% box on;
% ylabel('SIHR (\%)');
% set(gca,'FontSize',fs);
% legend('Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Wildtype','Covid Delta','SARS',...
%        'location','northwest');
% 
% subplot(2,2,4);
% xlim([1 17]);
% ylim([0 100]);
% xticks(1:17);
% xticklabels({'0-4','5-9','10-14','15-19',...
%              '20-24','25-29','30-34','35-39',...
%              '40-44','45-49','50-54','55-59',...
%              '60-64','65-69','70-74','75-79','80+'});
% xtickangle(45);
% grid on;
% box on;
% ylabel('HFR (\%)');
% set(gca,'FontSize',fs);
% legend('Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Wildtype','Covid Delta','SARS',...
%        'location','northwest');