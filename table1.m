%DON'T PLOT FIGURES!!!

close all;
clear all;

countries  = {'Rwanda',...
              'India','Egypt','Philippines','Indonesia',...
              'South Africa','Brazil',...
              'United Kingdom','Australia','United States'};
% countries  = cell(1,56);
% clist      = dir(fullfile('*.mat'));
% clist      = {clist.name};
% for j = 1:length(clist)-1;
%     name         = clist{j};
%     countries{j} = name(1:end-4);
% end
diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918',...
              'Covid Omicron','Covid Wildtype','Covid Delta',...
              'SARS'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};

gdp   = zeros(1,length(countries));
i1max = length(diseases);
i2max = length(strategies);

% xs = zeros(26,length(countries));
for j = 1:length(countries);
    inp1     = countries{j};
    load(strcat(inp1,'.mat'),'data');
    gdp(j)   = sum(365*data.obj);
%     %demography
%     xs(1,j)  = dot([2:5:102]',data.Npop)/sum(data.Npop);%average age
%     xs(2,j)  = sum(data.Npop(2:4))/sum(data.Npop);% % in each age group
%     xs(3,j)  = sum(data.Npop(5:13))/sum(data.Npop);
%     xs(4,j)  = sum(data.Npop(14:end))/sum(data.Npop);
%     xs(5,j)  = dot(data.la,[data.Npop(1:17);sum(data.Npop(18:end))])/sum(data.Npop);%remaining life expectancy, not life expectancy at birth
%     %employment
%     xs(6,j)  = sum(data.NNs(1:45))/sum(data.Npop(5:13));% % of adults in workforce
%     xs(7,j)  = sum(data.NNs(1:5))/sum(data.NNs(1:45));% % of workforce in primary
%     xs(8,j)  = sum(data.NNs(6:25))/sum(data.NNs(1:45));
%     xs(9,j)  = sum(data.NNs(26:45))/sum(data.NNs(1:45));
%     %mixing
%     xs(10,j) = data.comm/(dot(sum(data.CM,2),[data.Npop(1:15);sum(data.Npop(16:end))])/sum(data.Npop));% % of household contacts
%     xs(11,j) = data.travelA3/(dot(sum(data.CM,2),[data.Npop(1:15);sum(data.Npop(16:end))])/sum(data.Npop));
%     xs(12,j) = data.schoolA2/(dot(sum(data.CM,2),[data.Npop(1:15);sum(data.Npop(16:end))])/sum(data.Npop));
%     xs(13,j) = data.workp/(dot(sum(data.CM,2),[data.Npop(1:15);sum(data.Npop(16:end))])/sum(data.Npop));
%     %economics
%     xs(14,j) = sum(365*data.obj)/sum(data.Npop);%gva per capita
%     inds     = find(data.NNs(1:5));
%     NNs      = data.NNs(0+inds);
%     obj      = data.obj(0+inds);  
%     xs(15,j) = dot(365*obj./NNs,NNs)/sum(NNs);%gva per worker in primary
%     inds     = find(data.NNs(6:25));
%     NNs      = data.NNs(5+inds);
%     obj      = data.obj(5+inds);
%     xs(16,j) = dot(365*obj./NNs,NNs)/sum(NNs);
%     inds     = find(data.NNs(26:45));
%     NNs      = data.NNs(25+inds);
%     obj      = data.obj(25+inds);
%     xs(17,j) = dot(365*obj./NNs,NNs)/sum(NNs);
%     xs(18,j) = data.wpc;
%     %preparedness
%     xs(19,j) = data.Tres;
%     xs(20,j) = data.t_tit;
%     xs(21,j) = data.trate;
%     sd_fun   = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
%     xs(22,j) = integral(@(x) sd_fun(data.sdl,data.sdb,x),0,4);
%     xs(23,j) = data.Hmax;
%     xs(24,j) = data.t_vax;
%     xs(25,j) = data.arate;
%     xs(26,j) = data.puptake;
end
% % feature selection:
% % remaining life expectancy vs at birth?
% % adults in workforce or workforce in population?
% % David's PCA sectors
% % contact rates by age?
% xsz = zscore(xs,0,2);
% XS  = array2table(xs);
% XS.Properties.VariableNames  = countries;
% writetable(XS,'xs_raw.csv');
% XSZ = array2table(xsz);
% XSZ.Properties.VariableNames = countries;
% writetable(XSZ,'xs_zscore.csv');

%%

SEC  = zeros(length(countries),i1max,i2max);
VLYL = zeros(length(countries),i1max,i2max);
VSYL = zeros(length(countries),i1max,i2max);
GDPL = zeros(length(countries),i1max,i2max);

parfor j = 1:length(countries);   
    for i1 = 1:i1max;
        for i2 = 1:i2max;
        
        inp1 = countries{j};
        inp2 = diseases{i1};
        inp3 = strategies{i2}; 
        sec  = 100*p2Sim(inp1,inp2,inp3,'CURRENT',1,1)/gdp(j);
        
        SEC(j,i1,i2)  = round(sec(1),1);
        VLYL(j,i1,i2) = round(sec(2),1);
        VSYL(j,i1,i2) = round(sec(3),1);
        GDPL(j,i1,i2) = round(sec(4),1);
        
        end        
    end    
end

%%

SECt  = zeros(i1max*i2max,length(countries));
VLYLt = zeros(i1max*i2max,length(countries));
VSYLt = zeros(i1max*i2max,length(countries));
GDPLt = zeros(i1max*i2max,length(countries));

for j = 1:length(countries);   
    for i1 = 1:i1max;
        for i2 = 1:i2max;

        SECt(i2+(i1-1)*i2max,j)  = SEC(j,i1,i2);
        VLYLt(i2+(i1-1)*i2max,j) = VLYL(j,i1,i2);
        VSYLt(i2+(i1-1)*i2max,j) = VSYL(j,i1,i2);
        GDPLt(i2+(i1-1)*i2max,j) = GDPL(j,i1,i2);

        end
    end
end

writeformattedtable(SECt, 'Table1.csv', countries,diseases,strategies);
writeformattedtable(VLYLt,'TableA1.csv',countries,diseases,strategies);
writeformattedtable(VSYLt,'TableA2.csv',countries,diseases,strategies);
writeformattedtable(GDPLt,'TableA3.csv',countries,diseases,strategies);

function writeformattedtable(matrix,tablename,countries,diseases,strategies);
    
    celmat = cell(size(matrix));
    for i = 1:size(matrix,1);
        for j = 1:size(matrix,2);
            celmat(i,j) = cellstr(sprintf('%.5f',matrix(i,j)));%%%%%
        end
    end
    
%     SEC1B = cellstr(strcat('*',char(celmat(1+(1-1)*length(strategies),1))));
%     SEC1D = cellstr(strcat('*',char(celmat(2+(7-1)*length(strategies),5))));
%     SEC1F = cellstr(strcat('*',char(celmat(3+(5-1)*length(strategies),8))));
%     SEC1H = cellstr(strcat('*',char(celmat(4+(3-1)*length(strategies),9))));
%     celmat(1+(1-1)*4,1) = SEC1B;
%     celmat(2+(7-1)*4,5) = SEC1D;
%     celmat(3+(5-1)*4,8) = SEC1F;
%     celmat(4+(3-1)*4,9) = SEC1H;
    
    celmat = cell2table(celmat);
    celmat = [repelem(diseases',length(strategies),1),repmat(strategies',length(diseases),1),celmat];
    celmat.Var1(mod(1:length(diseases)*length(strategies),length(strategies))~=1) = {''};
    celmat.Properties.VariableNames(1)     = {'Var'};
    celmat.Properties.VariableNames(3:end) = countries;

    writetable(celmat,tablename);

end