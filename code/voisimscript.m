
if exist('p2SimRand')==0
    VOI_path = 'VOI';
    addpath(VOI_path);
end


diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918',...
              'Covid Omicron','Covid Delta',...
              'SARS','Covid Wildtype'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
income_levels = {'LLMIC','MIC','HIC'};

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

for i = 1:length(income_levels)
    for j = 1:length(diseases)
        for k = 1:length(strategies)
            inp2 = diseases{j};
            inp3 = strategies{k};
            income_level = income_levels{i};
            p2SimRand(inp2,inp3,income_level);
        end
    end
end
