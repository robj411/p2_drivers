function [ldata, configuration] = get_strategy_design(ldata,strategy,p2)

    % load in strategy parameters based on name
    % configuration includes all economic configurations to be used
    % hw includes three levels of home working
    % imand lists which configurations have a mandate applied
    % inext maps from events to configurations
    nSectors = ldata.nSectors;
    
    configuration = 0;
    if strcmp(strategy,'Elimination')
        configuration      = [ones(1*nSectors,1);ldata.x_econ(:,2);ldata.x_elim(:,1);ones(2*nSectors,1)];
        ldata.hw    = [zeros(1,nSectors);ldata.wfh(2,:);ldata.wfh(1,:);zeros(2,nSectors)];
        ldata.imand = [2];
        ldata.inext = [2,2,3,2,5];
    elseif strcmp(strategy,'Economic Closures')
        configuration      = [ones(1*nSectors,1);ldata.unmit;ldata.x_econ(:,2);ldata.x_econ(:,1);ones(nSectors,1)];
        ldata.hw    = [zeros(1,nSectors);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,nSectors)];
        ldata.imand = [3];
        ldata.inext = [2,3,3,4,5];
    elseif strcmp(strategy,'School Closures')
        configuration      = [ones(1*nSectors,1);ldata.unmit;ldata.x_schc(:,2);ldata.x_schc(:,1);ones(nSectors,1)];
        ldata.hw    = [zeros(1,nSectors);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,nSectors)];
        ldata.imand = [3];
        ldata.inext = [2,3,3,4,5];
    elseif strcmp(strategy,'No Closures')
        configuration      = [ones(1*nSectors,1); repmat(ldata.unmit,3,1); ones(1*nSectors,1)];
        ldata.hw    = [zeros(5,nSectors)];
        ldata.imand = [10];
        ldata.inext = [2,2,5];
    end
    
    ldata.tvec = [-0.1 ldata.tvec(end)];
    ldata.tvec = ldata.tvec + min(ldata.t_import,p2.Tres);

end
                