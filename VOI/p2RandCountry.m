function data = p2RandCountry(data,CD)

%Npop
nonempind = find(~isnan(CD.Npop1));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,4:24));
defivalue = 50*10^6*randvalue'/sum(randvalue);
data.Npop = defivalue;

%NNs
nonempind = find(~isnan(CD.NNs1));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,25:69));%number of workers by sector in real country
defivalue = randvalue/sum(table2array(CD(randindex,3+[5:13])));%proportion of adult population by sector in real country
defivalue = sum(data.Npop(5:13))*defivalue;%number of workers by sector in artificial country
defivalue = [defivalue,data.Npop(1),sum(data.Npop(2:4)),sum(data.Npop(5:13))-sum(defivalue),sum(data.Npop(14:end))]';
data.NNs  = defivalue;

%CM
nonempind = find(~isnan(CD.CMaa));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,70:325));
defivalue = reshape(randvalue,16,16);
data.CM   = defivalue;

%comm
nonempind = find(~isnan(CD.comm));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,326));
defivalue = randvalue;
data.comm = defivalue;

%travelA3
nonempind     = find(~isnan(CD.travelA3));
randindex     = nonempind(randi(numel(nonempind)));
randvalue     = table2array(CD(randindex,327));
defivalue     = randvalue;
data.travelA3 = defivalue;

%schoolA1&schoolA2
nonempind     = find(~isnan(CD.schoolA1));
randindex     = nonempind(randi(numel(nonempind)));
randvalue1    = table2array(CD(randindex,328));
randvalue2    = table2array(CD(randindex,329));
data.schoolA1 = randvalue1;
data.schoolA2 = randvalue2;

%workp
nonempind  = find(~isnan(CD.workp));
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,330));
defivalue  = randvalue;
data.workp = defivalue;

%obj
nonempind                   = find(~isnan(CD.obj1)&~isnan(CD.NNs1));
randindex                   = nonempind(randi(numel(nonempind)));
randvalue                   = table2array(CD(randindex,331:375));%gva by sector in real country
defivalue                   = randvalue./table2array(CD(randindex,25:69));%gva per worker by sector in real country
defivalue(isnan(defivalue)) = 0;
defivalue(isinf(defivalue)) = 0;
defivalue                   = data.NNs(1:45).*defivalue';%gva by sector in artificial country
data.obj                    = defivalue;

%wfh
nonempind = find(~isnan(CD.wfhl1));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,376:465));
defivalue = reshape(randvalue,45,2)';
data.wfh  = defivalue;

%Tres
nonempind = find(~isnan(CD.Tres));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,472));
defivalue = randvalue;
data.Tres = defivalue;

%t_tit
nonempind  = find(~isnan(CD.t_tit));
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,470));
defivalue  = randvalue;
data.t_tit = defivalue;

%trate
nonempind  = find(~isnan(CD.trate));
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,471));
defivalue  = randvalue;
data.trate = defivalue;

%sdl
nonempind = find(~isnan(CD.sdl));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,473));
defivalue = randvalue;
data.sdl  = defivalue;

%sdb
nonempind = find(~isnan(CD.sdb));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,474));
defivalue = randvalue;
data.sdb  = defivalue;

%Hmax
nonempind = find(~isnan(CD.Hmax));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,469));
defivalue = randvalue;
data.Hmax = defivalue;

%t_vax
nonempind  = find(~isnan(CD.t_vax));
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,466));
defivalue  = randvalue;
data.t_vax = defivalue;

%arate
nonempind  = find(~isnan(CD.arate));
randindex  = nonempind(randi(numel(nonempind)));
randvalue  = table2array(CD(randindex,467));
defivalue  = randvalue;
data.arate = defivalue;

%puptake
nonempind    = find(~isnan(CD.puptake));
randindex    = nonempind(randi(numel(nonempind)));
randvalue    = table2array(CD(randindex,468));
defivalue    = randvalue;
data.puptake = defivalue;

%la
nonempind = find(~isnan(CD.la1));
randindex = nonempind(randi(numel(nonempind)));
randvalue = table2array(CD(randindex,476:493));
defivalue = randvalue;
data.la   = defivalue;

%vly
na        = [data.Npop(1:17)',sum(data.Npop(18:end))];%length is 18 to match life table
gdp       = 365*sum(data.obj);
la        = data.la;
lg        = [dot(la(1),na(1))/sum(na(1)),...
             dot(la(2:4),na(2:4))/sum(na(2:4)),...
             dot(la(5:13),na(5:13))/sum(na(5:13)),...
             dot(la(14:end),na(14:end))/sum(na(14:end))];
for k = 1:length(lg); 
    lgh(k) = sum(1./((1+0.03).^[1:lg(k)]));
end 
vsl       = 160*gdp/sum(na);
defivalue = vsl/(dot(lgh,[na(1);sum(na(2:4));sum(na(5:13));sum(na(14:end))])/sum(na));
data.vly  = defivalue;

%vsy
defivalue = 0.5454*gdp/sum(data.Npop(2:4));
data.vsy  = defivalue;

end