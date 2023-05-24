function dis = get_dis_params(inp2)

    if strcmp(inp2,'Influenza 2009')
        dis = p2Params_Flu2009;
    elseif strcmp(inp2,'Influenza 1957')
        dis = p2Params_Flu1957;
    elseif strcmp(inp2,'Influenza 1918')
        dis = p2Params_Flu1918;
    elseif strcmp(inp2,'Covid Wildtype')
        dis = p2Params_CovidWT;    
    elseif strcmp(inp2,'Covid Omicron')
        dis = p2Params_CovidOM;    
    elseif strcmp(inp2,'Covid Delta')
        dis = p2Params_CovidDE;    
    elseif strcmp(inp2,'SARS')
        dis = p2Params_SARS;
    else
        error('Unknown Disease!');
    end
end

