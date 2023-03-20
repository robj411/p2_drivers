%% inputs

% Single Country
% countries = {'United Kingdom'};

% Ten Representative Countries
% countries = {'Rwanda',...
%              'India','Egypt','Philippines','Indonesia',...
%              'South Africa','Brazil',...
%              'United Kingdom','Australia','United States'};

% Countries with Most Data
%63 of which have population by sector
%56 of which have current preparedness data: no testing in brunei, no social distancing in brunei, china, cyprus, ethiopia, iceland, laos(deaths), russia
%some of which have population < 1m
filename  = '../../Data/Preparedness/5.gva_sector.xlsx';
T         = readtable(filename);
countries = T.ToIndustry_Sector;

% Countries with Missing Data
% filename  = '../../Data/Preparedness/1.population_age.xlsx';
% T         = readtable(filename);
% countries = T.Location;

% ARE CONTINUE STATEMENTS BEING USED???

j = 1;

%%

for i = 1:length(countries)%[115,124,120,116,119]
    
country = countries{i};

%% outputs

warning('off');
data = struct;

%% income group

filename = '../../Data/Preparedness/0.income_group.xlsx';
T        = readtable(filename);
kr       = find(strcmp(T.Var1,country));

%if isempty(kr); error(['Could not find income group of country: ',country]); end
if isempty(kr); 
    continue; 
end

year   = T(kr,:).Var2;
igroup = char(T(kr,:).Var3);
gnipc  = T(kr,:).Var40;

data.igroup = igroup;

%disp([country,' is a ',igroup,' and its GNI per capita is $',num2str(gnipc),' (World Bank, ',num2str(year),')']);

%% population by age

filename = '../../Data/Preparedness/1.population_age.xlsx';
T        = readtable(filename);
kr       = find(strcmp(T.Location,country));

if ~isempty(kr);
    kc   = find(strcmp(T.Properties.VariableNames,'x0_4'));
    Npop = 1000*table2array(T(kr,kc:kc+20))';
else
    continue;
    % Npop = zeros(21,1);
end

data.Npop = Npop;

%disp(['The population is ',num2str(sum(Npop)),' (United Nations, 2019)']);

%% population by sector
  
filename = '../../Data/Preparedness/2.population_sector.xlsx';
T        = readtable(filename);
kr       = find(strcmp(T.Country_,country));

if ~isempty(kr);
    kc     = find(strcmp(T.Properties.VariableNames,'Agriculture_Hunting_Forestry'));        
    year   = T(kr,:).Year_;
    NNs    = table2array(T(kr,kc:kc+44))';
    NEC    = T(kr,:).NotElsewhereClassified;
    source = 'OECD/ILO';
    %% informal employment    
    filename = '../../Data/Preparedness/2.population_sector_informal.xlsx';
    T        = readtable(filename);
    kr       = find(strcmp(T.Country,country));
    
    if ~isempty(kr);
        kc   = find(strcmp(T.Properties.VariableNames,'x4_A_Agriculture_ForestryAndFishing'));
        Ninf = table2array(T(kr,kc:kc+20))';
        ind  = [1,1,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,...
                4,5,6,7,8,8,8,8,8,9,10,10,10,11,12,13,14,15,16,17,18,19,20]';
        Ntot = accumarray(ind,NNs);
        disc = max(0,Ninf(1:20)-Ntot);
        for k = 1:45;
            indj   = ind(k);
            if Ntot(indj)>0;
                prop = NNs(k)/Ntot(indj);
            else
                prop = 1/sum(ind==indj);
            end
            NNs(k) = NNs(k) + disc(indj)*prop; 
        end
        NNs = NNs + max(NEC,Ninf(21))*NNs/sum(NNs);       
    else                
        NNs = NNs + NEC*NNs/sum(NNs); 
    end
    %% scaling
    if year<2019;
        filename = '../../Data/Preparedness/2.population_sector_change.xlsx';
        T        = readtable(filename);
        T.Properties.VariableNames = {'a','b','c','d','x2009','x2010','x2011','x2012','x2013',...
                                      'x2014','x2015','x2016','x2017','x2018','x2019'};
                                  
        kr   = find(strcmp(T.b,country));
        kc   = find(strcmp(T.Properties.VariableNames,['x',num2str(year)]));
        dpop = T(kr,:).x2019/table2array(T(kr,kc));
        NNs  = dpop*NNs;
    end
    %% capping
    if sum(NNs)>sum(Npop(5:13));
        NNs = sum(Npop(5:13))*NNs/sum(NNs);
    end
else
    continue;
    % year   = 0;
    % NNs    = zeros(45,1);
    % source = 'Estimated';
end

NNs = [NNs;Npop(1);sum(Npop(2:4));sum(Npop(5:13))-sum(NNs);sum(Npop(14:end))];
NNs = round(NNs);

data.NNs = NNs;

%disp(['The workforce accounts for ',num2str(round(100*sum(NNs(1:45))/sum(Npop(5:13)))),...
%      '% of the adult population, across ',num2str(nnz(NNs(1:45))),' sectors (',source,', ',num2str(year),')']);

%% contact matrix

filename = '../../Data/Preparedness/3.contact_matrices_all.xlsx';

if any(strcmp(sheetnames(filename),country));
    opts       = detectImportOptions(filename);
    opts.Sheet = country;
    T          = readtable(filename,opts);
    
    CM     = table2array(T);
    source = 'Prem et al., 2021';
else
    continue;
    % CM     = zeros(16,16);
    % source = 'Estimated, 0';
end

data.CM = CM;

%disp(['The average number of contacts per person per day is ',...
%      num2str(round(dot(sum(CM,2),[Npop(1:15);sum(Npop(16:end))])/sum(Npop),2)),' (',source,')']);

%% contact rates

if any(strcmp(sheetnames(filename),country));
    filename = '../../Data/Preparedness/4.contact_rates_home.xlsx';
    T        = table2array(readtable(filename,opts));
    comm     = dot(sum(T,2),[Npop(1:15);sum(Npop(16:end))])/sum(Npop);
    
    filename = '../../Data/Preparedness/4.contact_rates_other.xlsx';
    T        = table2array(readtable(filename,opts));
    travelA3 = dot(sum(T(5:13,5:13),2),Npop(5:13))/sum(Npop(5:13));
    
    filename = '../../Data/Preparedness/4.contact_rates_school.xlsx';
    T        = table2array(readtable(filename,opts));
    schoolA1 = T(1,1);
    schoolA2 = dot(sum(T(2:4,2:4),2),Npop(2:4))/sum(Npop(2:4));
    
    filename = '../../Data/Preparedness/4.contact_rates_work.xlsx';
    T        = table2array(readtable(filename,opts));
    workp    = dot(sum(T(5:13,:),2),Npop(5:13))/sum(Npop(5:13));
    
    source = 'Prem et al., 2021';
else
    comm     = 0;
    travelA3 = 0;
    schoolA1 = 0;
    schoolA2 = 0;
    workp    = 0;
    source   = 'Estimated, 0';
end

filename   = '../../Data/Preparedness/4.contact_rates_Beraud.xlsx';
opts       = detectImportOptions(filename);
opts.Sheet = 'Hospitality';
T          = readtable(filename,opts).Var2;
hospA2     = T(2);
hospA3     = T(3);
hospA4     = T(4);

opts.Sheet = 'B45';
T          = readtable(filename,opts).Var2;
B          = T(1:45);
    
opts.Sheet = 'C45';
T          = readtable(filename,opts).Var2;
C          = T(1:45);

data.comm     = comm;
data.hospA2   = hospA2;
data.hospA3   = hospA3;
data.hospA4   = hospA4;
data.travelA3 = travelA3;
data.schoolA1 = schoolA1;
data.schoolA2 = schoolA2;
data.workp    = workp;
data.B        = B';
data.C        = C';

%disp(['There are on average ',num2str(round(comm,2)),' contacts per person per day at home, ',num2str(round(travelA3,2)),..
%      ' during transport, ',num2str(round(schoolA2,2)),' at school and ',num2str(round(workp,2)),' at work (',source,')']);

%% gva by sector

filename = '../../Data/Preparedness/5.gva_sector.xlsx';
T        = readtable(filename);

kr = find(strcmp(T.ToIndustry_Sector,country));
if ~isempty(kr); 
    kc = find(strcmp(T.Properties.VariableNames,'D01T02_Agriculture_Hunting_Forestry'));
    
    obj    = table2array(T(kr,kc:kc+44))';
    source = 'OECD, 2018';
    %% scaling
    filename = '../../Data/Preparedness/5.gva_sector_change.xlsx';
    T        = readtable(filename);
    
    kr = find(strcmp(T.CountryName,country));
    
    dgva = T(kr,:).x2019/T(kr,:).x2018;
    if (dgva>0 && dgva<Inf);
        obj = dgva*obj;
    end
else
    continue;
    % obj    = zeros(45,1);
    % source = 'Estimated, 0';
end

data.obj = obj;

%disp(['The GDP is $',num2str(sum(365*obj)),' million (',source,')']);

%% seasonality

% filename = '../../Data/Preparedness/6.latitude.xlsx';
% T        = readtable(filename);
% 
% kr = find(strcmp(T.Country,country));
% if isempty(kr); error(['Could not find population-centroid latitude of country: ',country]); end
% 
% lat = T(kr,:).Latitude;
% amp = abs(-0.0214+0.1377*sin(lat/(180/(2*pi))));
% phi = (lat<0)*181;
% 
% data.amp = amp;
% data.phi = phi;
% 
% %disp(['The amplitude of seasonal forcing is ',num2str(amp),' (Hall et al., 2019; Douglas et al., 1997)']);

%% economic closures

filename   = '../../Data/Preparedness/6.economic_closures.xlsx';
opts       = detectImportOptions(filename);
% opts.Sheet = 'Alpha';
% T          = readtable(filename,opts);
% 
% kr = find(strcmp(T.country,country) & T.year==2019);
% if ~isempty(kr);
%     alp    = T(kr,:).labsh;
%     source = 'Feenstra et al., 2015';
%     
%     if isnan(alp);
%         alp    = 0;
%         source = 'Estimated, 0';
%     end
% else
%     alp    = 0;
%     source = 'Estimated, 0';
% end

opts.Sheet = 'Australia 45';
T          = readtable(filename,opts);
x_elim     = T.Var2;
x_elim(41) = 1.00;

opts.Sheet   = 'United Kingdom 45';
T            = readtable(filename,opts);
x_econ(:,1)  = T.Aug;
x_econ(:,2)  = T.Apr;
x_econ(41,1) = 1.00;
x_econ(41,2) = 0.10;

opts.Sheet   = 'Indonesia 45';
T            = readtable(filename,opts);
x_schc(:,1)  = T.TriwulanIV;
x_schc(:,2)  = T.TriwulanII;
x_schc(41,1) = 0.10;
x_schc(41,2) = 0.10;

% data.alp    = alp;
data.x_elim = x_elim;
data.x_econ = x_econ;
data.x_schc = x_schc;

%disp(['The labour share of productivity is ',num2str(alp),' (',source,')']);

%% home-working

filename   = '../../Data/Preparedness/7.home_working.xlsx';
opts       = detectImportOptions(filename);
opts.Sheet = 'LD vs. FO';
T          = readtable(filename,opts);
w05        = T.Var2(2:end);
w50        = T.Var3(2:end);
w95        = T.Var4(2:end);

opts.Sheet = 'Gottlieb';
T          = readtable(filename,opts);
kr         = find(strcmp(T.Var1,country));
if ~isempty(kr);
    w      = T(kr,:).Var2;
    scal   = w*sum(NNs(1:45))/dot(w50,NNs(1:45));
    wfh    = min(scal*[w05';w95'],1);
    source = 'Gottlieb et al., 2021';
else
    opts.Sheet         = 'GDPpcppp';
    opts.VariableNames = arrayfun(@num2str,1:68,'UniformOutput',0);
    opts.DataRange     = 'A4:BP270';
    T                  = readtable(filename,opts);
    kr                 = find(strcmp(T.x1,country));
    gdppcppp           = str2double(cell2mat(T(kr,:).x68));
    w                  = max(0,-0.2444+0.0791*log10(gdppcppp));
    scal               = w*sum(NNs(1:45))/dot(w50,NNs(1:45));
    wfh                = min(scal*[w05';w95'],1);
    source             = 'Estimated, 0';
end

data.wpc = w;
data.wfh = wfh;

%disp(['The share of home-working is ',num2str(w),' (',source,')']);

%% vaccination

filename = '../../Data/Preparedness/8.vaccination.csv';
T        = readtable(filename);

kr = strcmp(T.location,country);
if any(kr);    
    date = T(kr,:).date;
    dvec = datevec(date);
    dvec = dvec(:,1)-2020;
    date = day(date,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);
    Tts  = T(kr,:).people_vaccinated;
    Tts  = max(fillmissing(Tts,'linear'),0);
    v2   = T(kr,:).people_fully_vaccinated;
    v2   = max(fillmissing(v2,'linear'),0);
    Vts  = (Tts+v2)/2/sum(Npop/10^5);%this is for administration rate only!!!
    
    rng default;
    m       = (Vts(end)-Vts(find(Vts>0,1)))/(date(end)-date(find(Vts>0,1)));
    lb      = [300,                           350,       0];
    x0      = [max(300,date(end)-Vts(end)/m), date(end), Vts(end)];
    ub      = [800,                           date(end), 100000];
    %x0      = [500, 650,       60000];
    %x0      = [date(find(Vts>0,1)),date(end),Vts(end)];
    fun     = @(par,date)funky(par,date);
    options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1000000);
    problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',date,'ydata',Vts,...
                                 'lb',lb,'ub',ub,'options',options);
    ms      = MultiStart;
    poptim  = run(ms,problem,1); 
    
    %plot(date,Vts);
    %hold on;
    %plot(date,funky(poptim,date));
    
    t_vax   = poptim(1);
    arate   = poptim(3)/(poptim(2)-poptim(1));
    puptake = max(T(kr,:).people_fully_vaccinated_per_hundred)/100;
    source  = 'Our World in Data, 2022';
else
    t_vax   = 0;
    arate   = 0;
    puptake = 0;
    source  = 'Estimated, 0';
end

data.t_vax   = t_vax;
data.arate   = arate;
data.puptake = puptake;

%disp(['The vaccine rollout begins on day ',num2str(t_vax),...
%      ' with an administration rate of ',num2str(arate),' vaccines per 100k per day',...
%      ' and the uptake is ',num2str(100*puptake),'% (',source,')']);

%% hospital capacity

filename = '../../Data/Preparedness/9.hospital_capacity.xlsx';
opts     = detectImportOptions(filename);

opts.Sheet = 'Beds';
T          = readtable(filename,opts);
kr         = strcmp(T.Var1,country);
if any(kr);
    year   = T(kr,:).Var2;
    hcap   = 100*T(kr,:).Var3;
    source = 'World Bank/OECD';
else
    year   = 0;
    hcap   = 0;
    source = 'Estimated';
end

opts.Sheet = 'BOR';
T          = readtable(filename,opts);
kr         = strcmp(T.Var1,country);
if any(kr);
    BOR = T(kr,:).Var3/100;
else
    BOR = 0.85;%assumption
end

opts.Sheet = 'Covid';
T          = readtable(filename,opts);
kr         = strcmp(T.Var1,country);
if any(kr);
    Cmax = max(T(kr,:).Var3)/10;
else
    Cmax = 0;
end

data.Hmax  = hcap*(1-BOR);
%data.SHmax = hcap*(1-BOR)*2;

%disp(['The hospital capacity is ',num2str(hcap*(1-BOR)),' beds per 100k (',source,',',num2str(year),')']);

%% testing

filename = '../../Data/Preparedness/10.testing.csv';
T        = readtable(filename);

kr = strcmp(T.Entity,country);
if any(kr);    
    date = T(kr,:).Date;
    dvec = datevec(date);
    dvec = dvec(:,1)-2020;
    date = day(date,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);
    Tts  = 100*T(kr,:).CumulativeTotalPerThousand;
    Tts  = max(fillmissing(Tts,'linear'),0);
    
    fdat       = 1:850';
    fTts       = nan(size(fdat));
    fTts(date) = Tts;
    fTts       = max(fillmissing(fTts,'linear'),0);
    
    daily = movmean(diff(fTts),[0 30]);
    i_tit = find(daily>1,1);
    if ~isempty(i_tit);
        i_tit = i_tit+1;
        t_tit = fdat(i_tit);
    else
        i_tit = 367;
        t_tit = 367; 
    end
    t_end = 731;
    i_end = find(fdat==t_end,1);
    tspan = t_end-t_tit;
    Tspan = fTts(i_end)-fTts(i_tit-1);
    trate = Tspan/tspan;
    
    %     figure;
    %     hold on;
    %     scatter(date,Tts);
    %     plot(fdat,fTts,'r');
    % 
    %     figure;
    %     hold on;
    %     scatter(fdat(2:end),diff(fTts),'go');
    %     scatter(fdat(2:end),daily,'bo');
    %     scatter(t_tit,daily(i_tit-1),'r*');
    %     plot(t_tit:t_end,trate*ones(1,tspan+1));
    
    %     rng default;
    %     m       = (Tts(end)-Tts(find(Tts>0,1)))/(date(end)-date(find(Tts>0,1)));
    %     lb      = [1,                           50,        0];
    %     x0      = [max(1,date(end)-Tts(end)/m), date(end), Tts(end)];
    %     ub      = [800,                         date(end), 1000000];
    %     %x0      = [90,  date(end), 100000];
    %     %x0      = [date(find(Tts>0,1)),date(end),Tts(end)];
    %     fun     = @(par,date)funky(par,date);
    %     options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',500);
    %     problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',fun,'xdata',date,'ydata',Tts,...
    %                                  'lb',lb,'ub',ub,'options',options);
    %     ms      = MultiStart;
    %     poptim  = run(ms,problem,1); 
    %     
    %     figure;
    %     plot(date,Tts,'linewidth',2);
    %     hold on;
    %     plot([1:850],funky(poptim,[1:850]));
    %     plot([1:850],funky(x0,[1:850]));
    %     
    %     t_tit  = poptim(1);
    %     trate  = poptim(3)/(poptim(2)-poptim(1));
    %     source = 'Our World in Data, 2022';
else
    t_tit  = 0;
    trate  = 0;
    source = 'Estimated, 0';
end

% dur   = 1;
% Ip    = linspace(0,1000,500);
% trate = 110;
% b0    = 2.197;
% b1    = 0.1838;
% b2    = -1.024;
% p3    = (Ip<trate) .*   (1./(1+exp(b0+b1*Ip+b2*log10(trate))))/dur + ...
%         (Ip>=trate).*min(1./(1+exp(b0+b1*Ip+b2*log10(trate))),trate/10^5)/dur;
% p4    = p3;
% plot(Ip,p3);
% hold on;
% plot(trate*ones(1,11),[0:0.1:1],'r');
% plot(Ip,(trate/10^5)*ones(size(Ip)),'r');
% plot(Ip,(Ip.*p3)/trate,'g');%test positivity rate

data.t_tit = t_tit;
data.trate = trate;

%disp(['Mass testing begins on day ',num2str(t_tit),...
%      ' with an administration rate of ',num2str(trate),' tests per 100k per day','% (',source,')']);

%% response time

filename = '../../Data/Preparedness/11.response.csv';
T        = readtable(filename);

kc = find(strcmpi(T.Properties.VariableNames,strrep(strrep(strrep(country,' ',''),'-','_'),'''','_')));
if any(kc);    
    date   = T.country_name;
    dvec   = datevec(date);
    dvec   = dvec(:,1)-2020;
    date   = day(date,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);
    Bts    = table2array(T(:,kc));
    i_r    = find(Bts>=20,1);
    Tres   = date(i_r);
    source = 'Blavatnik, 2022';
else
    Tres   = 0;
    source = 'Estimated, 0';
end

data.Tres = Tres;

%disp(['The response time is on the ',num2str(Tres),' day of the first year (',source,')']);

%% social distancing

filename = '../../Data/Preparedness/12.mobility.csv';
T        = readtable(filename);
kr       = strcmp(T.Entity,country);

if any(kr);
    date = T(kr,:).Day;
    dvec = datevec(date);
    dvec = dvec(:,1)-2020;
    t2   = day(date,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);  
    Rts  = min(T(kr,:).retail_and_recreation,0);
    Gts  = min(T(kr,:).grocery_and_pharmacy,0);
    Mts  = 1+((Rts+Gts)/200);

    filename = '../../Data/Preparedness/12.excess_deaths.csv';
    T        = readtable(filename);
    kr       = strcmp(T.location_name,country);
    exd      = max(1,T(kr,:).mean_value);
    
    filename = '../../Data/Preparedness/12.deaths.csv';
    T        = readtable(filename);
    kc       = find(strcmpi(T.Properties.VariableNames,strrep(strrep(strrep(country,' ',''),'-','_'),'''','_')));
    date     = T.date;
    dvec     = datevec(date);
    dvec     = dvec(:,1)-2020;
    t1       = day(date,'dayofyear')+366*min(1,dvec)+365*max(0,dvec-1);
    Dts      = movmean(table2array(T(:,kc))/10,7);
    t1       = t1(~isnan(Dts));
    Dts      = exd*Dts(~isnan(Dts));
    
    [t,i1,i2] = intersect(t1,t2);
    Dts       = Dts(i1);
    Mts       = Mts(i2);
    
    low     = 1;
    upp     = 366;%+31+28+31;
    lu      = min(Mts);   
    tit_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
    tit_fit = fit([Dts(t>low&t<upp)],Mts(t>low&t<upp),tit_fun,...
              'StartPoint',[0.25,0.00000001],'Lower',[0.10,0],'Upper',[min(lu,0.4),10],...
              'Robust','LAR','MaxIter',10*10^3,'MaxFunEvals',10*10^3);

    % figure;
    % hold on;
    % scatter(Dts(t>low&t<upp),Mts(t>low&t<upp),'bo');
    % D = linspace(0,10,1000);      
    % plot(D,tit_fun(tit_fit.l,tit_fit.b,D),'b-');
    % xlim([0 10]);
    % ylim([0 1]);
    
    sdl    = tit_fit.l;
    sdb    = tit_fit.b;
    source = 'Our World in Data/Wang et al., 2022';
else
    continue;
    sdl    = 0;
    sdb    = 0;
    source = 'Estimated, 0';
end

data.sdl = sdl;
data.sdb = sdb;

%disp(['The transmission modifier asymptote is ',num2str(sdl),'(',source,')']);

%% costing - vlyl
 
%value of statistical life
filename = '../../Data/Preparedness/13.vly.csv';
T        = readtable(filename);
kr       = strcmp(T.DataSource,country);

if any(kr);
    year     = T(kr,:).Var3;
    gnipcppp = T(kr,:).Var69;
    vsl      = 160*gnipcppp/10^6;%100*gnipcppp/10^6;%in millions
    
    %discounted population-weighted life expectancy (disease independent)
    filename   = '../../Data/Preparedness/13.life_expectancy.csv';
    T          = readtable(filename);
    na         = [Npop(1:17)',sum(Npop(18:end))];%length is 18 to match life table
    for k = 1:18;
        label  = strcat('AGE',num2str(5*(k-1)),'-',num2str(5*(k-1)+4));
        index  = strcmp(T.Location,country) & strcmp(T.Dim2ValueCode,label);
        la(k)  = mean(T(index,:).Value,1); 
    end
    lg         = [dot(la(1),na(1))/sum(na(1)),...
                  dot(la(2:4),na(2:4))/sum(na(2:4)),...
                  dot(la(5:13),na(5:13))/sum(na(5:13)),...
                  dot(la(14:end),na(14:end))/sum(na(14:end))];
    for k = 1:length(lg); 
        lgh(k) = sum(1./((1+0.03).^[1:lg(k)]));
    end    
    
    %value of life-year
    %avage    = round(dot(Npop,[2:5:102]')/sum(Npop));
    %label    = strcat('AGE',num2str(avage-mod(avage,5)),'-',num2str(avage-mod(avage,5)+4));  
    %krv      = strcmp(T.Location,country)&strcmp(T.Dim2ValueCode,label);    
    %rlifex   = mean(T(krv,:).Value,1);
    %vly      = vsl/rlifex;
    vly    = vsl/(dot(lgh,[Npop(1);sum(Npop(2:4));sum(Npop(5:13));sum(Npop(14:end))])/sum(Npop));
    source = strcat('World Bank,',num2str(year));
else
    vly    = 0;
    la     = zeros(1,18);
    source = 'Estimated, 0';
end

data.vly = vly;
data.la  = la;

%disp(['The value of a life year is $',num2str(vly),' million (',source,')']);

%% costing - vsyl

% filename = '../../Data/15.vsy.xlsx';
% T        = readtable(filename);
% kr       = strcmp(T.CountryName,country);
% 
% if any(kr);
%     year   = cell2mat(T(kr,:).Year);
%     gdp    = T(kr,:).MOSTRECENT;
%     vsy    = 2.02*gdp/sum(Npop(2:4))/10^6;%in millions
%     source = strcat('World Bank,',num2str(year));
% else
%     vsy    = 0;
%     source = 'Estimated, 0';
% end

gdp = 365*sum(obj);
vsy = 0.5454*gdp/sum(Npop(2:4));;%2.02*gdp/sum(Npop(2:4));

if vsy~=0;
    source = 'OECD, 2018';
else
    source = 'Estimated, 0';
end

data.vsy = vsy;

%disp(['The value of a school year is $',num2str(vsy),' million (',source,')']);

%% costing - gdpl

% filename = '../../Data/Preparedness/15.consumption_sector.xlsx';
% T        = readtable(filename);
% kr       = strcmp(T.Var2,country);
% 
% if any(kr); 
%     kc     = find(strcmp(T.Properties.VariableNames,'TTL_01T02_Agriculture_Hunting_Forestry'));
%     hcon   = table2array(T(kr,kc:kc+44))';
%     source = 'OECD, 2018';
% else
%     hcon   = zeros(45,1);
%     source = 'Estimated, 0';
% end
% 
% filename = '../../Data/Preparedness/15.consumption_sector_sweden.xlsx';
% T        = readtable(filename);
% hconsl   = T(:,:).Var5;
% 
% data.hcon   = hcon;
% data.hconsl = hconsl;
% 
% %disp(['The annual household consumption expenditure is $',num2str(sum(365*hcon)),' million (',source,')']);

%% costing - implementation costs

% fixed = [0,          24.9333726, 0,              32.45842874, 484.2267634]';
% units = [0.09353734, 0.00006628, 0.000000150830, 0.00533906,  0.00007910]';
% 
% data.impcost = [fixed,units];

%% costing - PPP factor

% filename = '../../Data/Preparedness/14.ppp.xlsx';
% T        = readtable(filename);
% kr       = strcmp(T.CountryName,country);
% 
% if any(kr); 
%     year = T(kr,:).Year;
%     ppp  = T(kr,:).MOSTRECENT;
%     
%     filename = '../../Data/Preparedness/14.exchange_rate.xlsx';
%     T        = readtable(filename);
%     kr       = strcmp(T.CountryName,country);
%     
%     if any(kr);
%         %kc     = find(strcmp(T.Properties.VariableNames,strcat('x',year)));
%         %excr   = table2array(T(kr,kc));
%         excr   = T(kr,:).MOSTRECENT;
%         pppf   = ppp/excr;
%         source = strcat('World Bank, ',char(year));
%     else
%         pppf   = 0;
%         source = 'Estimated, 0';
%     end
% else
%     pppf   = 0;
%     source = 'Estimated, 0';
% end
% 
% data.pppf = pppf;
% 
% %disp(['The PPP conversion factor is ',num2str(pppf),' (',source,')']);

%% analysis

% if any(isnan([Npop;NNs;CM(:);comm;hospA2;hospA3;hospA4;travelA3;schoolA1;schoolA2;workp;B;C;obj])); 
%     error(['NaN in data for country: ',country]); 
% end

save(strcat(country,'.mat'),'data');

candidates{j,1}  = country;
candidates{j,2}  = igroup;
candidates{j,3}  = gnipc;
candidates{j,4}  = sum(Npop);%dot(Npop,[2:5:102]')/sum(Npop);
candidates{j,5}  = 100*sum(NNs(1:45))/sum(Npop(5:13));
candidates{j,6}  = 100*sum(NNs(1:5))/sum(NNs(1:45));
candidates{j,7}  = 100*sum(NNs(6:25))/sum(NNs(1:45));
candidates{j,8}  = 100*sum(NNs(26:45))/sum(NNs(1:45));
candidates{j,9}  = dot(sum(CM,2),[Npop(1:15);sum(Npop(16:end))])/sum(Npop);
candidates{j,10} = comm;
candidates{j,11} = travelA3;
candidates{j,12} = schoolA2;
candidates{j,13} = workp;
candidates{j,14} = 100*workp/(dot(sum(CM,2),[Npop(1:15);sum(Npop(16:end))])/sum(Npop));%this is NOT consistent
candidates{j,15} = 10^6*365*obj(1)/NNs(1);%this is NOT labour productivity
%candidates{j,16} = amp;
%candidates{j,17} = alp;
candidates{j,18} = w;
candidates{j,19} = t_vax;
candidates{j,20} = arate;
candidates{j,21} = 100*puptake;
candidates{j,22} = hcap;
candidates{j,23} = 100*BOR;
candidates{j,24} = hcap*(1-BOR);
candidates{j,25} = Cmax;
candidates{j,26} = t_tit;
candidates{j,27} = trate;
candidates{j,28} = Tres;
candidates{j,29} = sdl;
candidates{j,30} = sdb;
candidates{j,31} = vly;
candidates{j,32} = vsy;
%candidates{j,33} = 100*sum(hcon)/sum(obj);
%candidates{j,34} = pppf;

if mod(j,10)==0;
    display(j);
end
j = j+1;

end

%% postprocessing

for i = 1:size(candidates,1)
    
    hicc(i)  = strcmp(candidates{i,2},'HIC');
    umicc(i) = strcmp(candidates{i,2},'UMIC');
    lmicc(i) = strcmp(candidates{i,2},'LMIC');
    licc(i)  = strcmp(candidates{i,2},'LIC');
    
end

hics  = candidates(hicc,:);
umics = candidates(umicc,:);
lmics = candidates(lmicc,:);
lics  = candidates(licc,:);

[~,ihic]  = sort([hics{:,3}], 'ascend','MissingPlacement','first'); 
[~,iumic] = sort([umics{:,3}],'ascend','MissingPlacement','first'); 
[~,ilmic] = sort([lmics{:,3}],'ascend','MissingPlacement','first'); 
[~,ilic]  = sort([lics{:,3}], 'ascend','MissingPlacement','first'); 

hics  = hics(ihic,:);
umics = umics(iumic,:);
lmics = lmics(ilmic,:);
lics  = lics(ilic,:);

candidates = [lics;lmics;umics;hics];

%% plotting

names = categorical(candidates(:,1));
names = reordercats(names,cellstr(candidates(:,1)));

f  = figure('Units','centimeters','Position',[0 0 45 15]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',10);
fs = 12;

g = gca;
g.Position = [0.035 0.28 0.93 0.625];
hold on;

b = bar(names,cell2mat(candidates(:,3)));
%plot(names,70*ones(1,length(candidates)),'-','linewidth',1,'color','k');

%ylim([0 100000]);
%ylabel('GNI per capita (\$)');
g.XAxis.TickLength = [0 0];
xtickangle(45);
%g.YAxis.Exponent = 3;
grid on;
box on;
set(gca,'FontSize',fs);

b.FaceColor = 'flat';
b.CData(1:size(lics,1),:)                                                        = repmat([0 1 1],size(lics,1),1);
b.CData(size(lics,1)+1:size(lics,1)+size(lmics,1),:)                             = repmat([0 .5 1],size(lmics,1),1);
b.CData(size(lics,1)+size(lmics,1)+1:size(lics,1)+size(lmics,1)+size(umics,1),:) = repmat([0 0 .5],size(umics,1),1);
b.CData(size(lics,1)+size(lmics,1)+size(umics,1)+1:end,:)                        = repmat([0.41 0.16 0.38],size(hics,1),1);

%% levels

% %save P2levels.mat candidates;
% load('P2levels.mat');
% 
% pop  = cell2mat(candidates(:,4));
% test = cell2mat(candidates(:,21));
% test = test(pop>10^6 & test~=0);%is it appropriate to delete 0s?
% l    = quantile(test,0.125+[0.00,0.25,0.50,0.75]);
% 
% figure;
% hold on;
% histogram(test,100,'FaceColor','yellow','FaceAlpha',0.00);
% xline(l(1),'LineWidth',1.5,'Color','k');
% xline(l(2),'LineWidth',1.5,'Color','r');
% xline(l(3),'LineWidth',1.5,'Color','b');
% xline(l(4),'LineWidth',1.5,'Color','g');
% 
% % D       = linspace(0,10,1000);      
% % tit_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
% % figure;
% % hold on;
% % plot(D,tit_fun(m(4),l(4),D),'k-');
% % plot(D,tit_fun(m(3),l(3),D),'r-');
% % plot(D,tit_fun(m(2),l(2),D),'b-');
% % plot(D,tit_fun(m(1),l(1),D),'g-');

%% modelling

% if strcmp(igroup,'HIC');  
%     replacement = 'Italy';  
%     
% elseif strcmp(igroup,'UMIC');  
%     replacement = 'Guyana'; 
% 
% elseif strcmp(igroup,'LMIC');  
%     replacement = 'Angola'; 
% 
% elseif strcmp(igroup,'LIC'); 
%     replacement = 'Sudan';  
%     
% else
%     error(['Could not find replacement for country: ',country]);
%     
% end

%     filename = '../../Data/1.population_age.xlsx';
%     T        = readtable(filename);
%     kr       = find(strcmp(T.Location,replacement));
%     kc       = find(strcmp(T.Properties.VariableNames,'x0_4'));
% 
%     Nad = sum(1000*table2array(T(kr,kc+4:kc+12)));
% 
%     T  = readtable(filename_e);
%     kr = find(strcmp(T{:,1},replacement));
%     kc = find(strcmp(T.Properties.VariableNames,cname_e));
% 
%     pNs    = table2array(T(kr,kc:kc+44))';
%     if width(T)~=49;
%         pNs    = pNs/Nad;
%     else
%         pNs    = (pNs + T(kr,:).Var41*(pNs/sum(pNs)))/Nad;
%     end            
%     NNs    = sum(Npop(5:13))*pNs;

% b.CData(find(names=='Sudan'),:)  = [1 0 0];
% b.CData(find(names=='Angola'),:) = [1 0 0];
% b.CData(find(names=='Guyana'),:) = [1 0 0];
% b.CData(find(names=='Italy'),:)  = [1 0 0];

% b.CData(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Italy'),5)),2),:)...
% =repmat([1 1 1],sum(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Italy'),5)),2)),1);
% b.CData(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Guyana'),5)),2),:)...
% =repmat([1 1 1],sum(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Guyana'),5)),2)),1);
% b.CData(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Angola'),5)),2),:)...
% =repmat([1 1 1],sum(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Angola'),5)),2)),1);
% b.CData(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Sudan'),5)),2),:)...
% =repmat([1 1 1],sum(round(cell2mat(candidates(:,5)),2)==round(cell2mat(candidates(find(names=='Sudan'),5)),2)),1);

%% functions

function V=funky(params,dates)
    syms ddate;
    V = piecewise(ddate<params(1),0,...
                  params(1)<ddate<params(2),...
                  -params(1)*(params(3)/(params(2)-params(1))) + ddate*(params(3)/(params(2)-params(1))),...
                  ddate>params(2),params(3));
                          
    V = subs(V,ddate,dates);
    V = fillmissing(double(V),'nearest');
end