% updates parameters relating to hospital outcomes as a function of
% hospital occupancy
%
% occ: hospital occupancy
% p2: struct of p2 intervention parameters
% dis: struct of pathogen parameters
% t: current time
%
% dis2: updated struct of pathogen parameters

function dis2 = update_hosp_dis_parameters(occ, p2, dis, t)

    dis2 = dis;
    
    %if t<max(p2.tpoints)
    
        Hmax  = p2.Hmax;
        new_pd = (1 + dis.rel_FR.*max(0, occ - Hmax) / occ).*dis.pd;

        %Probabilities
        pd = min(new_pd,1);

        %Calculations
        Threc = dis.Threc;
        Thd = dis.Thd;

        Th = ((1-pd).*Threc) + (pd.*Thd);
        dis2.g3   = (1-pd)./Th;
        dis2.mu   = pd./Th;
    %end

end



