function [namevec, vals] = tabulate_inputs(ldata,p2,dis2)

    % contacts
    leavenames = {'CM','basic_contact_matrix'};
    [contactname, contactvals] = get_struct_names_vals(ldata.contacts,leavenames);

    % dis
    keepnames = {'CI','beta','Td','generation_time','ps','Tlat','Tay','Tsr','Tsh','Threc','Thd','R0','frac_presymptomatic','ihr','ifr','hfr'};
    leavenames = setdiff(fieldnames(dis2),keepnames);
    [disname, disvals] = get_struct_names_vals(dis2,leavenames);


    % p2
    leavenames = {'arate','final_doubling_time_threshold','group_order','hosp_final_threshold','sdb','sdl','t_tit','t_vax2','time_to_test','tpoints'};
    [p2name, p2vals] = get_struct_names_vals(p2,leavenames);


    % data
    leavenames = {'EdInd','HospInd','Npop','adInd','compindex','contacts','educationloss_all_students','gdppc','lgh',...
        'nSectors','nStrata','response_time','t_vax','tvec','x_unmit','vaccination_rate_pc','vaccine_uptake','x_econ','x_elim','x_schc'};
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
