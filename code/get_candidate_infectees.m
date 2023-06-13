
addpath('VOI');

load('Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
fields    = fieldnames(data);
ikeep     = [6,7,8,13,14,16,17,18];
data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
data.tvec = [-75 365*3+1];
CD        = readtable('country_data.csv');
nsamples  = 10000;

inp2s = {'Influenza 2009','Influenza 1957','Influenza 1918','Covid Wildtype','Covid Omicron','Covid Delta','SARS'};
    
for j = 1:length(inp2s) 
    
    dis = get_dis_params(inp2s{j});   

    CIs = [];

    for i = 1:nsamples;

        ldata     = data;
        ldata     = p2RandCountry(ldata,CD,'all');

        [ldata,~,~]    = p2Params(ldata,'Covid Wildtype',dis);%to define wnorm and Td_CWT
        CI = p2getCandidateInfectees(ldata,dis);
        CIs(i) = CI;
    end
%     inp2s{j}
%     min(R0s)
%     max(R0s)
    [qR0,R0] = ksdensity(CIs,1:0.1:100,'Function','cdf');
    uqR0 = find(unique(qR0));
    qR0 = qR0(uqR0);
    R0 = R0(uqR0);
%     interp1(qR0,R0,.5)
    save(sprintf('%sR0.mat',inp2s{j}),'qR0','R0')

end

tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
nexttile   

[~,qi] = min(abs(qR0-0.5));
[~,ql] = min(abs(qR0-0.1));
[~,qh] = min(abs(qR0-0.9));
plot(R0,[0 diff(qR0)],'linewidth',3)
set(gca,'YTickLabel',[]);
ylabel('','FontSize',16);
xlabel('Candidate infectees','FontSize',16);
xline(R0(qi),'--','linewidth',3,'color','#BABCBC');
xline(R0(ql),'--','linewidth',3,'color','#BABCBC');
xline(R0(qh),'--','linewidth',3,'color','#BABCBC');

nexttile
Rs = 0:0.01:6;
plot(Rs,normpdf(Rs,3,.3),'linewidth',3,'color','#FF0038')
set(gca,'YTickLabel',[]);
ylabel('','FontSize',16);
xlabel('Basic reproduction number','FontSize',16);
xline(3,'--','linewidth',3,'color','#BABCBC');
xline(norminv(0.1,3,.3),'--','linewidth',3,'color','#BABCBC');
xline(norminv(0.9,3,.3),'--','linewidth',3,'color','#BABCBC');


recovery_times = [];
incubation_times = [];
for i = 1:nsamples

    ldata     = data;
    ldata     = p2RandCountry(ldata,CD,'all');

    [ldata,~,~]    = p2Params(ldata,'Covid Wildtype',dis);%to define wnorm and Td_CWT
    recovery_time = normrnd(2.5,0.25);
    incubation_time = normrnd(2.5,0.25);
    recovery_times(i) = recovery_time;
    incubation_times(i) = incubation_time;
    dis.Tlat = incubation_time;
    dis.Tay = incubation_time;
    dis.Tsh = recovery_time;
    dis.Tsr = recovery_time;
    CI = p2getCandidateInfectees(ldata,dis);
    CIs(i) = CI;
end
scatter(incubation_times,CIs)
scatter(recovery_times,CIs)
scatter3(incubation_times,recovery_times,CIs)
scatter(incubation_times,recovery_times,10,CIs,'filled')
