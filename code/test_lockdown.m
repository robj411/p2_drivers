function R = test_lockdown(p2,ldata,dis2)

    lx = ldata.lx;
    NN = ldata.NNs;
    workpop = sum(NN(1:lx));
    NN(1:lx) = NN(1:lx).*ldata.x_econ(:,2);
    NN(lx+ldata.adInd) = NN(lx+ldata.adInd) + workpop - sum(NN(1:lx));

    betamod = p2.sdl;
    % assume cases per 100k
    Ip = 10000; % = 10%
    p3 = get_case_ID_rate(p2, Ip); 
    p4 = p3; 

    
    R = zeros(1,3);
    R(3) = p3;
    
    contact_matrix = p2MakeDs(ldata,NN,ldata.x_econ(:,2),ldata.wfh(2,:));

    R(1) = get_R(ldata.ntot, dis2, dis2.h, dis2.g2, NN, NN, contact_matrix, dis2.beta, betamod, p3, p4);
       
    
    contact_matrix = p2MakeDs(ldata,NN,ldata.x_schc(:,2),ldata.wfh(2,:));


    R(2) = get_R(ldata.ntot, dis2, dis2.h, dis2.g2, NN, NN, contact_matrix, dis2.beta, betamod, p3, p4);

end
