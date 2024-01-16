function [ldata, xoptim] = get_strategy_design(ldata,strategy,p2)
    lx = ldata.lx;
    xoptim = 0;
    if strcmp(strategy,'Elimination')
        xoptim      = [ones(1*lx,1);ldata.x_econ(:,2);ldata.x_elim(:,1);ones(2*lx,1)];
        ldata.hw    = [zeros(1,lx);ldata.wfh(2,:);ldata.wfh(1,:);zeros(2,lx)];
        ldata.imand = [2];
        ldata.inext = [2,2,3,2,5];
    elseif strcmp(strategy,'Economic Closures')
        xoptim      = [ones(1*lx,1);ldata.unmit;ldata.x_econ(:,2);ldata.x_econ(:,1);ones(lx,1)];
        ldata.hw    = [zeros(1,lx);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,lx)];
        ldata.imand = [3];
        ldata.inext = [2,3,3,4,5];
    elseif strcmp(strategy,'School Closures')
        xoptim      = [ones(1*lx,1);ldata.unmit;ldata.x_schc(:,2);ldata.x_schc(:,1);ones(lx,1)];
        ldata.hw    = [zeros(1,lx);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,lx)];
        ldata.imand = [3];
        ldata.inext = [2,3,3,4,5];
    elseif strcmp(strategy,'No Closures')
        xoptim      = [ones(1*lx,1); repmat(ldata.unmit,3,1); ones(1*lx,1)];
        ldata.hw    = [zeros(5,lx)];
        ldata.imand = [10];
        ldata.inext = [2,2,5];
    else
        error('Unknown Mitigation Strategy!');
    end
    
    ldata.tvec = [-0.1 ldata.tvec(end)];
%     if strcmp(countrytype,'Origin')
%         ldata.t_import = 0;
%     elseif strcmp(countrytype,'Secondary')
        ldata.tvec = ldata.tvec + min(ldata.t_import,p2.Tres);
%     end

end
                