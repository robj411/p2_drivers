function f=p2Plot(data,f1,p2,g1,cost,ccost_t,sec,inp1,inp2,inp3)

ln        = length(data.NNs);
lx        = length(data.obj);
tvec      = data.tvec;
thresh    = p2.Hmax;
limit     = p2.SHmax;
numThresh = length(thresh);

t1    = f1(:,1); 
I1    = f1(:,2);
h1    = f1(:,3);
d1    = f1(:,4);
ddiff = diff(d1,1);
tdiff = diff(t1,1);
d1    = [0;ddiff./tdiff];
p3    = f1(:,5);
beta  = f1(:,6);
v1    = f1(:,7);
v2    = f1(:,8);
v3    = f1(:,9);
v4    = f1(:,10);

scal4 = sum(data.Npop)/(10^4);
scal5 = sum(data.Npop)/(10^5);
scal6 = sum(data.Npop)/(10^6);
scal7 = sum(data.Npop)/(10^7);
scal8 = sum(data.Npop)/(10^8);
scal9 = sum(data.Npop)/(10^9);
%maxY  = max([100000,(5/4)*d1'/scal5,(5/4)*thresh/scal4,(5/4)*h1'/scal4,(5/4)*I1'/scal3]);
maxY  = 60000;%ceil(max([d1'/scal5,thresh/scal4,h1'/scal4,I1'/scal3])/10000)*10000;

T               = repmat(t1',lx+1,1);
S               = 0.5:1:lx+0.5;
w               = g1(:,1+[1:lx])';
x               = w.^data.alp;
x(data.EdInd,:) = w(data.EdInd,:);
X               = [x;ones(1,length(t1))];

%% EPIDEMIC TRAJECTORY

f  = figure('Units','centimeters','Position',[0 0 21 9]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultAxesColorOrder',[[0 0 0];[0 0 0]]);
set(f,'DefaultAxesFontSize',12);
fs = 12;
lw = 2;

ax = gca;
ax.Position = [0.05 0.20 0.90 0.70];
hold on;

for i = 1:length(tvec)-1; 
    a     = [tvec(i) tvec(i) tvec(i+1) tvec(i+1)];
    b     = [0 maxY maxY 0];
    pind  = tvec(i)<t1 & t1<tvec(i+1);
    px    = x(:,pind);
    alpha = 1-mean(px.^3,'all');
    if alpha == 0;
        continue;
    end
    fill(a,b,'yellow','facealpha',alpha,'linewidth',0.01,'EdgeColor','k');
end

yyaxis left;
hh3 = plot([-30,tvec(end)],[thresh,thresh]/scal8,'--','linewidth',lw,'color',0.5*[1,1,1]);
hh4 = plot(t1,d1/scal9,'-','linewidth',lw,'color','black');
hh2 = plot(t1,h1/scal8,'-','linewidth',lw,'color','magenta');
% hh1 = plot(t1,I1/scal3,'-','linewidth',lw,'color','red');

yyaxis right;
hh5 = plot(t1,100*(v1+v2+v3+v4)/sum(data.Npop),':','linewidth',lw,'color','cyan');

xlim([-30 tvec(end)]);
set(gca,'xtick',[[1,91,182,274],...
             365+[1,91,182,274],...
           2*365+[1,91,182,274]]);
set(gca,'xticklabels',{'Year 1','Apr','Jul','Oct',...
                       'Year 2','Apr','Jul','Oct',...
                       'Year 3','Apr','Jul','Oct'});
xtickangle(45);
yyaxis left;
ylim([0 60000]);
ax                   = gca;
ax.YAxis(1).Exponent = 3;
yyaxis right;
ylim([0 100]);
ylabel('\%');
set(get(gca,'ylabel'),'Rotation',0);
vec_pos = get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos + [-60 60 0]);

grid on;
box on;
legend([hh2,hh3,hh4,hh5],'Hospital Occupancy (per 100m)','Hospital Capacity (per 100m)',... %%%%%
                         'Daily Deaths (per 1b)','Vaccine Coverage (\%)','location','north');
t = text(-22.3,86.6,strvcat(inp1,inp2,inp3),'FontSize',fs);
t.BackgroundColor = [1 1 1];
t.EdgeColor       = [0 0 0];
set(gca,'FontSize',fs);

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

f = figure('Units','centimeters','Position',[0 0 6 9]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'DefaultAxesFontSize',12);

ax          = gca;
ax.Position = [0.25 0.20 0.65 0.70];

labs = categorical(["VLYL","GDPL","VSYL"]);
labs = reordercats(labs,[2 1 3]);
y    = 100*[cost(3,lx+1) cost(3,lx+2)   sum(cost(3,[1:lx,lx+3])) cost(3,ln);...
        sum(cost(7,:))   sum(cost(8,:)) sum(cost(9,:))           sum(cost(10,:));...
        cost(6,lx+1:ln)]/sum(365*data.obj);
b    = bar(labs,y,'stacked','FaceColor','flat');

xtickangle(45);
ymax = 100*max([sum(cost(3,:)),sum(cost(6,:)),sum(cost(7:10,:),'all')])/sum(365*data.obj);
ylmt = ceil(ymax/10)*10;
ylmt = ylmt + 10*((ylmt-ceil(ymax))<(0.25*ylmt));%%%%%
ylim([0 ylmt]);
ylabel('Societal Costs');
vec_pos = get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos + [-0.3 0 0]);
text(0.6,ylmt*1.05,'\%','FontSize',fs);

grid on;
grid minor;
box on;
b(1).CData = [0.00 0.00 0.00;...
              1.00 0.00 1.00;...
              1.00 0.00 1.00];
b(2).CData = [0.00 0.00 0.00;...
              1.00 1.00 0.00;...
              1.00 1.00 0.00];
b(3).CData = [0.00 0.00 0.00;...
              0.00 0.00 0.00;...
              0.00 0.00 0.00];
b(4).CData = [0.00 0.00 0.00;...
              0.00 0.00 0.00;...
              0.00 0.00 0.00];
set(gca,'FontSize',fs);

ax.YTick(end) = [];%%%%%
axes('Position',[.22 .81 .05 .05]);
px     = [1 5];
py1    = [1 2];
height = 1;
py2    = py1 + height;
plot(px,py1,'k','LineWidth',2);hold all;
plot(px,py2,'k','LineWidth',2);hold all;
fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','none');
box off;
axis off;

end