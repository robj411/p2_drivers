function sec = p2Sim(inp1,inp2,inp3)
    
    load(strcat(inp1,'.mat'),'data');
    lx        = length(data.B);
    data.tvec = [-75 365*3+1];
    
    
    dis = get_dis_params(inp2);    
    
    [data,~,~]    = p2Params(data,'Covid Wildtype',dis);%to define wnorm and Td_CWT
    [data,dis,p2] = p2Params(data,inp2,dis);
    
    int = 5;
    if strcmp(inp3,'Elimination');
        xoptim     = [ones(1*lx,1);data.x_econ(:,2);data.x_elim(:,1);ones(2*lx,1)];
        data.hw    = [zeros(1,lx);data.wfh(2,:);data.wfh(1,:);zeros(2,lx)];
        data.imand = [2];
        data.inext = [2,2,3,2,5];
    elseif strcmp(inp3,'Economic Closures');
        xoptim     = [ones(2*lx,1);data.x_econ(:,2);data.x_econ(:,1);ones(lx,1)];
        data.hw    = [zeros(1,lx);data.wfh(1,:);data.wfh(2,:);data.wfh(1,:);zeros(1,lx)];
        data.imand = [3];
        data.inext = [2,3,3,4,5];
    elseif strcmp(inp3,'School Closures');
        xoptim     = [ones(2*lx,1);data.x_schc(:,2);data.x_schc(:,1);ones(lx,1)];
        data.hw    = [zeros(1,lx);data.wfh(1,:);data.wfh(2,:);data.wfh(1,:);zeros(1,lx)];
        data.imand = [3];
        data.inext = [2,3,3,4,5];
    elseif strcmp(inp3,'No Closures');
        xoptim     = [ones(5*lx,1)];
        data.hw    = [zeros(5,lx)];
        data.imand = [10];
        data.inext = [2,2,5];
    else
        error('Unknown Mitigation Strategy!');
    end
    
    [data,f,g] = p2Run(data,dis,inp3,int,xoptim,p2);
    [cost,~]   = p2Cost(data,dis,p2,g);
    sec(1)     = sum(cost([3,6,7:10],:),'all');
    sec(2)     = sum(cost([3],:),'all');
    sec(3)     = sum(cost([6],:),'all');
    sec(4)     = sum(cost([7:10],:),'all');

    p2Plot(data,f,p2,g,cost,NaN,sec(1),inp1,inp2,inp3);
    
end