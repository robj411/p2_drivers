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
    
    % only apply before vaccination programme completes (prevents expensive
    % exit waves)
    %if t<max(p2.tpoints)
    
        Hmax  = p2.Hmax;
%         SHmax = 2*Hmax;
%         excess_H = max(0, occ - Hmax);
%         p_excess = excess_H ./ occ;
%         (1 + 0.87*max(0, occ - Hmax) / occ)*pd;
%         new_pd = ((1 - p_excess) + p_excess.*1.87).*dis.pd;
        new_pd = (1 + 9*max(0, occ - Hmax) / occ).*dis.pd;

        %Amplitudes
%         th0 = max(1,1+1.87*((occ-Hmax)./(SHmax-Hmax)));
%         new_pd = th0.*dis.pd;

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



