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
    if t<max(p2.tpoints)
    
        Hmax  = p2.Hmax;
        SHmax = 2*Hmax;

        %Amplitudes
        th0 = max(1,1+1.87*((occ-Hmax)./(SHmax-Hmax)));

        %Probabilities
        pd = min(th0*dis.pd,1);

        %Calculations
        Threc = dis.Threc;
        Thd = dis.Thd;

        Th = ((1-pd).*Threc) + (pd.*Thd);
        dis2.g3   = (1-pd)./Th;
        dis2.mu   = pd./Th;
    end

end



