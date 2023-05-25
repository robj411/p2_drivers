function f=p2Video(data,f1,p2,g1,cost,ccost_t)

ln        = length(data.NNs);
lx        = length(data.obj);
tvec      = data.tvec;
thresh    = p2.Hmax;
limit     = p2.SHmax;
numThresh = length(thresh);

t0    = f1(:,1); 
s1    = f1(:,2); 
I1    = f1(:,3);
h1    = f1(:,4);
d1    = f1(:,5);
ddiff = diff(d1,1);
tdiff = diff(t0,1);
d1    = [0;ddiff./tdiff];
v1    = f1(:,6);
v2    = f1(:,7);
v3    = f1(:,8);
v4    = f1(:,9);
r1    = f1(:,14);

t1 = [t0(1):0.2:t0(end)];
I1 = interp1(t0,I1,t1)';
h1 = interp1(t0,h1,t1)';
d1 = interp1(t0,d1,t1)';
v1 = interp1(t0,v1,t1)';
v2 = interp1(t0,v2,t1)';
v3 = interp1(t0,v3,t1)';
v4 = interp1(t0,v4,t1)';

scal0 = sum(data.Npop)/(10^4);
scal1 = sum(data.Npop)/(10^5);
scal2 = sum(data.Npop)/(10^6);
scal3 = sum(data.Npop)/(10^7);
scal4 = sum(data.Npop)/(10^8);
scal5 = sum(data.Npop)/(10^9);
maxY  = max([100000,(5/4)*d1'/scal5,(5/4)*thresh/scal4,(5/4)*h1'/scal4,(5/4)*I1'/scal3]);
%maxY  = 250000;

T               = repmat(t1',lx+1,1);
S               = 0.5:1:lx+0.5;
w               = g1(:,1+[1:lx])';
w               = interp2(t0,1:45,w,t1',1:45);
x               = w.^data.alp;
x(data.EdInd,:) = w(data.EdInd,:);
X               = [x;ones(1,length(t1))];

ccost_t         = interp2(t0,1:size(ccost_t,2),ccost_t',t1',1:size(ccost_t,2))';

%% EPIDEMIC TRAJECTORY

f  = figure('Units','centimeters','Position',[0 0 25 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;
lw = 2;

ax = gca;
ax.Position = [0.05 0.20 0.90 0.70];
hold on;

%for i = 1:numThresh;
hh6 = plot([-80,tvec(end)],[thresh,thresh]/scal4,'--','linewidth',lw,'color',0.5*[1,1,1]);
%plot([-80,tvec(end)],[limit(j),limit(j)]/scal4,':','linewidth',lw,'color',.5*[1,1,1]);
%end

inc  = 20;
pl1  = 1;
plv1 = ones(size(tvec,2)-1);
plv2 = ones(size(tvec,2)-1);
% anin = [2*inc,39*inc,46*inc,57*inc,69*inc,81*inc,...%UK,CV,EC
%         96*inc,107*inc,112*inc,172*inc,224*inc];
% xx   = [0.65,0.51;...
%         0.77,0.91;...
%         0.77,0.91;...
%         0.80,0.91;...
%         0.77,0.91;...
%         0.80,0.91;...
%         0.77,0.91;...
%         0.82,0.93;...
%         0.82,0.93;...
%         0.82,0.93;...
%         0.82,0.93];
% yy   = [0.60,0.25;...
%         0.60,0.25;...
%         0.60,0.35;...
%         0.60,0.28;...
%         0.60,0.35;...
%         0.60,0.28;...
%         0.60,0.35;...
%         0.60,0.28;...
%         0.60,0.25;...
%         0.60,0.50;...
%         0.60,0.55];
% antx = {'disease seeded','government responds','lockdown 1','reopening 1','lockdown 2','reopening 2',...
%         'lockdown 3','reopening 3','vaccine rollout starts','vaccine rollout ends','exit wave'};
% anin = [2*inc,13*inc,18*inc,126*inc,185*inc,238*inc];%SI,SP,EL CURRENT PREPAREDNESS
% xx   = [0.65,0.51;...
%         0.70,0.84;...
%         0.70,0.84;...
%         0.75,0.89;...
%         0.78,0.92;...
%         0.75,0.89];
% yy   = [0.60,0.25;...
%         0.60,0.28;...
%         0.60,0.28;...
%         0.60,0.38;...
%         0.60,0.38;...
%         0.60,0.38];
% antx = {'disease seeded',sprintf(['border closure\n','lockdown\n']),sprintf(['lockdown\n','lifted\n']),...
%         'vaccine rollout starts',sprintf(['vaccine rollout ends\n','border reopens\n']),'exit wave'};
% maxY = 1.2*maxY;
% anin = [2*inc,20*inc,41*inc,89*inc,135*inc,146*inc,201*inc,213*inc];%SI,SP,EL LOWER PREPAREDNESS
% xx   = [0.65,0.51;...
%         0.70,0.84;...
%         0.75,0.89;...
%         0.75,0.93;...
%         0.78,0.92;...
%         0.80,0.93;...
%         0.80,0.93;...
%         0.85,0.91];
% yy   = [0.60,0.25;...
%         0.60,0.48;...
%         0.55,0.23;...
%         0.60,0.25;...
%         0.60,0.25;...
%         0.60,0.30;...
%         0.60,0.48;...
%         0.60,0.50];
% antx = {'disease seeded',sprintf(['late lockdown\n','capacity breached\n']),sprintf(['prolonged\n','lockdown lifted\n']),...
%         sprintf(['importation triggers\n','second lockdown']),sprintf(['slower vaccine\n','rollout']),sprintf(['importation triggers\n','third lockdown']),...
%         sprintf(['border\n','reopens']),sprintf(['low vaccine coverage\n','large exit wave'])};
% anin = [2*inc,42*inc,48*inc,132*inc,177*inc];%ET,SW,UN,L2
% xx   = [0.65,0.51;...
%         0.70,0.84;...
%         0.70,0.88;...
%         0.78,0.92;...
%         0.80,0.93];
% yy   = [0.60,0.25;...
%         0.60,0.30;...
%         0.55,0.25;...
%         0.60,0.25;...
%         0.60,0.30];
% antx = {'disease seeded',sprintf(['wave of infections\n','fewer hospitalisations']),sprintf(['hospital occupancy\n','within capacity\n']),...
%         sprintf(['second wave due \n','to waning immunity']),sprintf(['vaccine rollout'])};
anin = [2*inc,42*inc,48*inc,123*inc,177*inc,200*inc];%ET,SW,UN,L3
xx   = [0.65,0.51;...
        0.70,0.84;...
        0.70,0.88;...
        0.78,0.92;...
        0.80,0.91;...
        0.80,0.93];
yy   = [0.60,0.25;...
        0.60,0.30;...
        0.55,0.25;...
        0.60,0.30;...
        0.60,0.35;...
        0.60,0.425];
antx = {'disease seeded',sprintf(['slightly smaller\n','wave of infections']),sprintf(['hospital occupancy\n','well within capacity\n']),...
        sprintf(['smaller second wave']),sprintf(['faster vaccine\n','rollout']),sprintf(['higher vaccine\n','coverage'])};
% anin = [2*inc,27*inc,40*inc,48*inc,169*inc,205*inc,232*inc];%VN,CV,SC
% xx   = [0.65,0.51;...
%         0.70,0.84;...
%         0.75,0.90;...
%         0.73,0.91;...
%         0.78,0.92;...
%         0.78,0.92;...
%         0.75,0.93];
% yy   = [0.60,0.25;...
%         0.60,0.38;...
%         0.70,0.58;...
%         0.70,0.35;...
%         0.60,0.28;...
%         0.80,0.75;...
%         0.70,0.38];
% antx = {'disease seeded',sprintf(['home-working mandate\n','social distancing\n']),sprintf(['lockdown\n','schools closed\n']),'lockdown lifted',...
%         'vaccine rollout starts',sprintf(['vaccine rollout ends\n','schools reopen\n']),'exit wave'};

vd1           = VideoWriter('scenario_traj.avi');
vd1.FrameRate = 10;
open(vd1);

for i = [2*inc:inc:length(t1)-inc,length(t1)-inc];
    
    if (t1(i)>=tvec(2)) && (pl1==1);
        plot(tvec(2)*[1,1],[0,maxY],'k-','linewidth',0.01);
        pl1 = 0;
    end
    for j = 1:length(plv1);
    if (t1(i)>=tvec(j)) && (t1(i)<tvec(j+1)) && (plv1(j)==1);
        a       = [tvec(j) tvec(j) tvec(end) tvec(end)];
        b       = [-100 maxY+100 maxY+100 -100];
        pind    = tvec(j)<t1 & t1<tvec(j+1);
        px      = x(:,pind);
        alpha   = 1-mean(px,'all');
        shd     = fill(a,b,'yellow','linewidth',0.01,'facealpha',alpha);
        plv1(j) = 0;
    elseif (t1(i)>=tvec(j+1)) && (plv2(j)==1);
        set(shd,'Visible','off');
        a       = [tvec(j) tvec(j) tvec(j+1) tvec(j+1)];
        b       = [-100 maxY+100 maxY+100 -100];
        pind    = tvec(j)<t1 & t1<tvec(j+1);
        px      = x(:,pind);
        alpha   = 1-mean(px,'all');
        fill(a,b,'yellow','linewidth',0.01,'facealpha',alpha);
        plv2(j) = 0;
    end
    end
    hh5 = plot(t1(i-inc:i),(v1(i-inc:i)+v2(i-inc:i)+v3(i-inc:i)+v4(i-inc:i))/scal1,'-','linewidth',lw,'color','green');
    hh4 = plot(t1(i-inc:i),d1(i-inc:i)/scal5,'-','linewidth',lw,'color','black');
    hh3 = plot(t1(i-inc:i),h1(i-inc:i)/scal4,'-','linewidth',lw,'color','magenta');
    hh2 = plot(t1(i-inc:i),I1(i-inc:i)/scal3,'-','linewidth',lw,'color','red');
    mrk1 = plot(t1(i),(v1(i)+v2(i)+v3(i)+v4(i))/scal1,'o','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','green');
    mrk2 = plot(t1(i),d1(i)/scal5,'o','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','black');
    mrk3 = plot(t1(i),h1(i)/scal4,'o','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','magenta');
    mrk4 = plot(t1(i),I1(i)/scal3,'o','MarkerSize',7,'MarkerEdgeColor','black','MarkerFaceColor','red');

    ylmt = max(5*[(v1(1:i)+v2(1:i)+v3(1:i)+v4(1:i))'/scal1,d1(1:i)'/scal5,h1(1:i)'/scal4,I1(1:i)'/scal3]);   
    ylmt = max(ylmt,100);
    ylmt = min(ylmt,maxY);
    axis([t1(inc),t1(i)+t1(1+inc)-t1(1),0,ylmt]);
    %xlabel('Time','FontSize',fs);
    %ylabel('Number','FontSize',fs);%yvar
    %vec_pos=get(get(gca,'ylabel'),'Position');
    %set(get(gca,'ylabel'),'Position',vec_pos+[-10 0 0]);
    set(gca,'xtick',[[1,91,182,274],...
                 365+[1,91,182,274],...
               2*365+[1,91,182,274]]);
    set(gca,'xticklabels',{'Year 1 - Jan','Apr','Jul','Oct',...
                           'Year 2 - Jan','Apr','Jul','Oct',...
                           'Year 3 - Jan','Apr','Jul','Oct'});
    xtickangle(45);
    ax                = gca;
    ax.YAxis.Exponent = 3;
    grid on;
    box on;
    set(gca,'FontSize',fs);
    
    legend([hh2,hh3,hh6,hh4,hh5],'Prevalence (per 10m)','Hospital Occupancy (per 100m)','Hospital Capacity (per 100m)',...
                                 'Daily Deaths (per 1b)','Vaccinated (per 100k)','location','northwest');%'Position',[-0.29 0.27 1 1]);

    if ismember(i,anin);
        ind   = find(i==anin);
        a     = annotation('textarrow',xx(ind,:),yy(ind,:),'LineWidth',1,'String',antx{ind},'FontSize',fs,'interpreter','latex');
        drawnow;%pause();
        frame = getframe(gcf);
        for k = 1:30
            writeVideo(vd1,frame);
        end
        set(a,'Visible','off');
    else    
        drawnow;
        frame = getframe(gcf);
        writeVideo(vd1,frame);
    end
    set(mrk1,'Visible','off');
    set(mrk2,'Visible','off');
    set(mrk3,'Visible','off');
    set(mrk4,'Visible','off');

end

close(vd1);

%% MITIGATION MEASURES

% f = figure('Units','centimeters','Position',[0 0 10 10]);
% set(f,'defaulttextInterpreter','latex');
% set(f,'defaultAxesTickLabelInterpreter','latex');
% set(f,'defaultLegendInterpreter','latex');
% set(f,'defaultColorbarTickLabelInterpreter','latex');
% set(f,'DefaultAxesFontSize',12);
% 
% ax = gca;
% ax.Position = [0.175 0.20 0.80 0.70];
% 
% h = pcolor(T,S,X);
% 
% axis([-30,3*365,0.5,lx+0.5]);
% %axis square;
% set(gca,'layer','top');
% set(gca,'xtick',[[1,182],...
%              365+[1,182],...
%            2*365+[1,182]]);
% set(gca,'ytick',0.5+[5:5:lx]);
% grid on;
% box on;
% set(gca,'xticklabels',{'Year 1 - Jan','Jul',...
%                        'Year 2 - Jan','Jul',...
%                        'Year 3 - Jan','Jul'});       
% set(gca,'yticklabels',{'5','10','15','20','25','30','35','40','45'});
% xtickangle(45);
% %xlabel('Time');
% ylabel('Sector')
% set(h,'EdgeColor','none');
% set(gca,'FontSize',fs);
%
% colormap(hot);
% caxis([0,1]);
% colorbar;

%% COSTS

f = figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'DefaultAxesFontSize',12);

ax          = gca;
ax.Position = [0.15 0.20 0.80 0.70];

labs = categorical(["Lost Lives","Business Closures","School Closures","Direct Mitigation"]);
labs = reordercats(labs,[3 1 4 2]);
ymax = max([sum(ccost_t(end,1:ln)),sum(ccost_t(end,2*ln+1:4*ln)),sum(ccost_t(end,ln+1:2*ln)),sum(ccost_t(end,4*ln+1:end))]);
expo = floor(log10(ymax));
ylmt = ceil(ymax/10^expo)*10^expo;

vd2           = VideoWriter('scenario_cost.avi');
vd2.FrameRate = 10;
open(vd2);

for i = [2*inc:inc:length(t1)-inc,length(t1)];

    y    = [ccost_t(i,lx+1)             ccost_t(i,lx+2)             sum(ccost_t(i,[1:lx,lx+3])) ccost_t(i,ln) 0;...
            sum(ccost_t(i,2*ln+[1:ln])) sum(ccost_t(i,3*ln+[1:ln])) 0                           0             0;...
            ccost_t(i,ln+lx+[1:4])      0;...
            ccost_t(i,4*ln+[1:5])];
    b    = bar(labs,y,'stacked','FaceColor','flat');

    ylim([0 ylmt]);
    %axis square;
    xtickangle(28);
    grid on;
    grid minor;
    box on;
    ylabel('Socio-Economic Costs (\$, millions)');
    set(gca,'FontSize',fs);

    b(1).CData = [0.65 0.00 0.13;...
                  0.00 0.00 0.50;...
                  0.80 0.63 0.21;...
                  0.00 0.00 0.00];
    b(2).CData = [0.90 0.00 0.15;...
                  0.00 0.00 1.00;...
                  1.00 0.83 0.00;...
                  0.30 0.30 0.30];
    b(3).CData = [1.00 0.00 0.00;...
                  0.00 0.50 1.00;...
                  1.00 0.83 0.00;...
                  0.50 0.50 0.50];
    b(4).CData = [1.00 0.39 0.28;...
                  0.00 1.00 1.00;...
                  1.00 0.83 0.00;...
                  0.70 0.70 0.70];
    b(5).CData = [0.00 0.00 0.00;...
                  0.00 0.00 0.00;...
                  0.00 0.00 0.00;...
                  0.90 0.90 0.90];
    sec = sum(ccost_t(i,:));
    ylm = ylim;
    title(['\textbf{Socio-Economic Cost:} ' num2str(round(100*sec/sum(365*data.obj))) '\% of annual GDP']);
    set(get(gca,'title'),'Position',[2.25 1.065*ylm(2) 0]);
    
    drawnow;
    frame = getframe(gcf);
    if ismember(i,anin);
        for k = 1:30
            writeVideo(vd2,frame);
        end
    else    
        writeVideo(vd2,frame);
    end
    
end

close(vd2);

%%

vid1 = vision.VideoFileReader('scenario_traj.avi');
vid2 = vision.VideoFileReader('scenario_cost.avi');
vidP = vision.VideoPlayer;

video_object           = VideoWriter('scenario.avi');
video_object.FrameRate = 10;
open(video_object);

while ~isDone(vid1);
   frame1 = step(vid1);
   frame2 = step(vid2);
   frame  = horzcat(frame1, frame2);
   writeVideo(video_object,frame);
   step(vidP,frame);
end

release(vid1);
release(vid2);
release(vidP);
close(video_object);
delete 'scenario_traj.avi';
delete 'scenario_cost.avi';

end