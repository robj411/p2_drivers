% prepares final inputs for model run

% p2: struct of p2 intervention parameters
% data: struct of general model parameters
% strategy: name of one of four mitigation strategies

% data: struct of general model parameters
% configuration: matrix of economic configurations: numbers between 0 and 1
% giving the extent of closure by sector

function data = get_strategy_design(data,strategy,p2)

    % load in strategy parameters based on name
    % configuration includes all economic configurations to be used
    % hw includes three levels of home working
    % imand lists which configurations have a mandate applied
    % inext maps from events to configurations
    nSectors = data.nSectors;
    
    configuration = 0;
    x_unmit = data.x_unmit;
    x_elim = data.x_elim;
    x_schc = data.x_schc;
    x_econ = data.x_econ;
    wfh = data.wfh;
    if strcmp(strategy,'Elimination')
        configuration      = [ones(1*nSectors,1);x_econ(:,2);x_elim(:,1);ones(nSectors,1);ones(nSectors,1);ones(nSectors,1)];%;x_econ(:,1)];
        data.hw    = [zeros(1,nSectors);wfh(2,:);wfh(1,:);zeros(1,nSectors);zeros(1,nSectors);wfh(1,:)];%;zeros(1,nSectors)];%
        data.imand = [2];
        data.inext = [3,2,3,2,5,5]; % [2,2,3,2,5,5];
    elseif strcmp(strategy,'Economic Closures')
        configuration      = [ones(1*nSectors,1);x_unmit;x_econ(:,2);x_econ(:,1);ones(nSectors,1);ones(nSectors,1)];
        data.hw    = [zeros(1,nSectors);wfh(1,:);wfh(2,:);wfh(1,:);zeros(1,nSectors);wfh(1,:)];
        data.imand = [3];
        data.inext = [4,3,3,4,5,5]; % [2,3,3,4,5,5]
    elseif strcmp(strategy,'School Closures')
        configuration      = [ones(1*nSectors,1);x_unmit;x_schc(:,2);x_schc(:,1);ones(nSectors,1);ones(nSectors,1)];
        data.hw    = [zeros(1,nSectors);wfh(1,:);wfh(2,:);wfh(1,:);zeros(1,nSectors);wfh(1,:)];
        data.imand = [3];
        data.inext = [4,3,3,4,5,5]; % [2,3,3,4,5,5]
    elseif strcmp(strategy,'No Closures')
        configuration      = [ones(1*nSectors,1); repmat(x_unmit,3,1); ones(1*nSectors,1);ones(nSectors,1)];
        data.hw    = [zeros(5,nSectors);zeros(1,nSectors)];
        data.hw    = [zeros(1,nSectors);wfh(1,:);wfh(2,:);wfh(1,:);zeros(1,nSectors);wfh(1,:)];
        data.imand = [10];
        data.inext = [2,2,2,2,5,5];
    end
    
    data.tvec = [-0.1 data.tvec(end)];
    data.tvec = data.tvec + min(data.t_import,p2.Tres);
    data.configuration = configuration;

    if diff(p2.group_order(1:2))==-4
        % disp(data.tvec)
        data.tvec = [data.tvec(1), p2.tpoints(1) + 1e-5, data.tvec(2)];
    end
end
               
               
               
                
