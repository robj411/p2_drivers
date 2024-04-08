function [namevec, vals] = tabulate_inputs(ldata,p2,dis2)

    % contacts
    leavenames = {'CM','basic_contact_matrix'};
    [contactname, contactvals] = get_struct_names_vals(ldata.contacts,leavenames);

    % dis
    keepnames = {'CI','beta','Td','generation_time','ps','Tlat','Tay','Tsr','Tsh','Threc','Thd','R0','frac_presymptomatic','ihr','ifr'};
    leavenames = setdiff(fieldnames(dis2),keepnames);
    [disname, disvals] = get_struct_names_vals(dis2,leavenames);


    % p2
    keepnames = {'Hmax','thl','trate','frac_sym_infectiousness_averted','frac_presym_infectiousness_averted','frac_asym_infectiousness_averted','Tres' };
    leavenames = setdiff(fieldnames(p2),keepnames);
    [p2name, p2vals] = get_struct_names_vals(p2,leavenames);


    % data
    keepnames = {'remote_quantile','response_time_quantile','remote_teaching_effectiveness','self_isolation_compliance',...
        'sd_baseline','sd_death_coef','sd_mandate_coef','labsh' ,'NNs','t_import', 'Hres' ,'la' ,'obj','gdp','vly','vsy','frac_tourism_international'};
    leavenames = setdiff(fieldnames(ldata),keepnames);
    ldata.obj = ldata.obj/ldata.gdp;
    ldata.vly = ldata.vly/ldata.gdp;
    ldata.vsy = ldata.vsy/ldata.gdp;
    
    [dataname, datavals] = get_struct_names_vals(ldata,leavenames);

    % join

    vals = [contactvals disvals p2vals datavals];
    namevec = [contactname disname p2name dataname];
    
end


function [namevec, vals] = get_struct_names_vals(stct,leavenames)

    namevec = {};
    vals = [];

    nms = fieldnames(stct);

    for i = 1:numel(nms)
        nm = nms{i};
        if sum(strcmp(nm, leavenames))==0
            item = stct.(nm);
            dims = size(item);
            if dims(1)==1
                dims = dims(2);
            end
            nentries = prod(dims);
            if nentries==1
                vals = [vals item];
                namevec{length(namevec)+1} = nm;
            else
                vals = [vals item(:)'];
                x = 1:nentries;
                for j = x
                    namevec{length(namevec)+1} = strcat(nm,num2str(j));
                end
            end
        end
    end

end
