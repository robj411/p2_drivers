% parpool;
% 
diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918',...
              'Covid Omicron','Covid Wildtype','Covid Delta',...
              'SARS'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
% 
% for j = 1:length(diseases);
% for k = 1:length(strategies);
%     
%     inp2 = diseases(j);
%     inp3 = strategies(k);
%     T    = p2SimRand(NaN,inp2,inp3);
% 
% end
% end
% 
% delete(gcp);

%%

for j = 1:length(diseases);
for k = 4%:length(strategies);
    
    
    f  = figure('Units','centimeters','Position',[0 0 27 9]);
    set(f,'defaulttextInterpreter','latex');
    set(f,'defaultAxesTickLabelInterpreter','latex');
    set(f,'defaultLegendInterpreter','latex');
    set(f,'DefaultAxesFontSize',12);
    fs = 12;
    lw = 2;

    % ax = gca;
    % ax.Position = [0.05 0.20 0.90 0.70];
    hold on;
    
    inp2 = diseases(j);
    inp3 = strategies(k);
    T    = readtable(strcat('VOI_',string(inp2),'_',string(inp3),'.csv'));

    gdp  = 365*sum(table2array(T(:,331:375)),2);
    
    swarmchart(1*ones(1,1000),(100*T.SEC./gdp),'k.');
    swarmchart(2*ones(1,1000),(100*T.VLYL./gdp),'r.');
    swarmchart(3*ones(1,1000),(100*T.VSYL./gdp),'g.');
    swarmchart(4*ones(1,1000),(100*T.GDPL./gdp),'b.');
    plot(linspace(0.5,4.5,100),40*ones(1,100),'k-');
    %set(gca,'YScale','log');
    xlim([0.5 4.5]);
    %ylim([0 150]);
    grid on;
    grid minor;
    xticks([1:4]);
    xticklabels({'SEC','VLYL','VSYL','GDPL'});
    xtickangle(45);    
    
    %subplot(length(diseases),length(strategies),(j-1)*length(strategies)+k);
    
end
end

%%

% filename = 'VOI_ECD.csv';
% T        = readtable(filename);
% 
% Tres     = T.Var1;
% t_tit    = T.Var2;     
% trate    = T.Var3;     
% sdl      = T.Var4;      
% sdb      = T.Var5; 
% Hmax     = T.Var6;     
% t_vax    = T.Var7;     
% arate    = T.Var8;     
% puptake  = T.Var10;
% 
% SEC      = T.Var11;
% VLYL     = T.Var12;
% VSYL     = T.Var13;
% GDPL     = T.Var14;
% 
% figure;
% subplot(331);
% scatter(Tres,SEC,'.');
% R = corrcoef(Tres,SEC);
% title(round(R(1,2),2));
% %
% subplot(332);
% scatter(t_tit,SEC,'.');
% R = corrcoef(t_tit,SEC);
% title(round(R(1,2),2));
% %
% subplot(333);
% scatter(trate,SEC,'.');
% R = corrcoef(trate,SEC);
% title(round(R(1,2),2));
% %
% subplot(334);
% scatter(sdl,SEC,'.');
% R = corrcoef(sdl,SEC);
% title(round(R(1,2),2));
% %
% subplot(335);
% scatter(log(sdb),SEC,'.');
% R = corrcoef(log(sdb),SEC);
% title(round(R(1,2),2));
% %
% subplot(336);
% scatter(Hmax,SEC,'.');
% R = corrcoef(Hmax,SEC);
% title(round(R(1,2),2));
% %
% subplot(337);
% scatter(t_vax,SEC,'.');
% R = corrcoef(t_vax,SEC);
% title(round(R(1,2),2));
% %
% subplot(338);
% scatter(arate,SEC,'.');
% R = corrcoef(arate,SEC);
% title(round(R(1,2),2));
% %
% subplot(339);
% scatter(puptake,SEC,'.');
% R = corrcoef(puptake,SEC);
% title(round(R(1,2),2));
% 
% 
% 
% 
% 
% 
% figure;
% subplot(331);
% scatter(Tres,VLYL,'.');
% R = corrcoef(Tres,VLYL);
% title(round(R(1,2),2));
% %
% subplot(332);
% scatter(t_tit,VLYL,'.');
% R = corrcoef(t_tit,VLYL);
% title(round(R(1,2),2));
% %
% subplot(333);
% scatter(trate,VLYL,'.');
% R = corrcoef(trate,VLYL);
% title(round(R(1,2),2));
% %
% subplot(334);
% scatter(sdl,VLYL,'.');
% R = corrcoef(sdl,VLYL);
% title(round(R(1,2),2));
% %
% subplot(335);
% scatter(log(sdb),VLYL,'.');
% R = corrcoef(log(sdb),VLYL);
% title(round(R(1,2),2));
% %
% subplot(336);
% scatter(Hmax,VLYL,'.');
% R = corrcoef(Hmax,VLYL);
% title(round(R(1,2),2));
% %
% subplot(337);
% scatter(t_vax,VLYL,'.');
% R = corrcoef(t_vax,VLYL);
% title(round(R(1,2),2));
% %
% subplot(338);
% scatter(arate,VLYL,'.');
% R = corrcoef(arate,VLYL);
% title(round(R(1,2),2));
% %
% subplot(339);
% scatter(puptake,VLYL,'.');
% R = corrcoef(puptake,VLYL);
% title(round(R(1,2),2));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% figure;
% subplot(331);
% scatter(Tres,GDPL,'.');
% R = corrcoef(Tres,GDPL);
% title(round(R(1,2),2));
% %
% subplot(332);
% scatter(t_tit,GDPL,'.');
% R = corrcoef(t_tit,GDPL);
% title(round(R(1,2),2));
% %
% subplot(333);
% scatter(trate,GDPL,'.');
% R = corrcoef(trate,GDPL);
% title(round(R(1,2),2));
% %
% subplot(334);
% scatter(sdl,GDPL,'.');
% R = corrcoef(sdl,GDPL);
% title(round(R(1,2),2));
% %
% subplot(335);
% scatter(log(sdb),GDPL,'.');
% R = corrcoef(log(sdb),GDPL);
% title(round(R(1,2),2));
% %
% subplot(336);
% scatter(Hmax,GDPL,'.');
% R = corrcoef(Hmax,GDPL);
% title(round(R(1,2),2));
% %
% subplot(337);
% scatter(t_vax,GDPL,'.');
% R = corrcoef(t_vax,GDPL);
% title(round(R(1,2),2));
% %
% subplot(338);
% scatter(arate,GDPL,'.');
% R = corrcoef(arate,GDPL);
% title(round(R(1,2),2));
% %
% subplot(339);
% scatter(puptake,GDPL,'.');
% R = corrcoef(puptake,GDPL);
% title(round(R(1,2),2));