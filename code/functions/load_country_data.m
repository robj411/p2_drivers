function [CD, country_parameter_distributions] = load_country_data()
    CD        = readtable('../data/country_data.csv');
    CD.popsum = sum(table2array(CD(:,4:24)),2);
    
    CMcols = find(contains(fieldnames(CD),'CM'));
    dim = sqrt(length(CMcols));
    Npopcols = find(contains(fieldnames(CD),'Npop'));
    average_contacts = NaN(size(CD,1),1);
    for i=1:size(CD,1)
        cm = table2array(CD(i,CMcols));
        cmm = reshape(cm,16,[]);
        Npop = table2array(CD(i,Npopcols));
        Npop(dim) = sum(Npop(dim:length(Npop)));
        Npop = Npop(1:dim);
        if ~isnan(cm(1))
            average_contacts(i) = sum(sum(cmm,2) .* Npop') / sum(Npop); % max(eig(cmm));%
        end
    end
    CD.average_contacts = average_contacts;
    
    country_parameter_distributions = readtable('../data/parameter_distributions.csv');
    
end
