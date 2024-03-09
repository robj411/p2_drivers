% updates parameters relating to hospitalisation as a function of
% waning vaccine immunity, and past infection

% dis: struct of pathogen parameters
% S: number susceptible by stratum
% Sn: number susceptible and not previously infected or vaccinated by stratum
% compindex: struct of compartment indices 
% y_mat: matrix of present values in all compartments

% dis2: updated struct of pathogen parameters
% V: amount of BPSV in population
% B: amount of SARS-X vaccine in population
% vaccine_pp: amount of BPSV per person
% booster_pp: amount of SARS-X vaccine per person

function [dis2, V, B, vaccine_pp, booster_pp] = update_vax_dis_parameters(dis, S, Sn, compindex, y_mat)

    V =    y_mat(:,compindex.V_index(1));
    B =    y_mat(:,compindex.V_index(2));
    vaccinated_people = sum(y_mat(:,compindex.vaccine),2) + 1e-16;
    boosted_people = sum(y_mat(:,compindex.booster),2) + 1e-16;
    vaccine_pp = V./vaccinated_people;
    booster_pp = B./boosted_people;
    
    dis2 = dis;
    
    ph = dis.ph;
    Tsr = dis.Tsr;
    Tsh = dis.Tsh;
    
    amp   = (Sn+(1-dis.hv2).*(S-Sn))./S;
    naive_adjusted_ph    = amp.*ph;
    
    Ts    = ((1-naive_adjusted_ph).*dis.Tsr) + (naive_adjusted_ph.*dis.Tsh);
    
    dis.g2    = (1-naive_adjusted_ph)./Ts;
    dis.h     = naive_adjusted_ph./Ts;
    
    dis2.trv1 = dis.trv1 .* vaccine_pp;
    dis2.trv2 = dis.trv2 .* booster_pp;
    dis2.scv1 = dis.scv1 .* vaccine_pp;
    dis2.scv2 = dis.scv2 .* booster_pp;
    
    hv1 = dis.hv1 .* vaccine_pp;
    hv2 = dis.hv2 .* booster_pp;
    Ts_v1 = ((1-(1-hv1).*ph).*Tsr)  +((1-hv1).*ph.*Tsh);
    Ts_v2 = ((1-(1-hv2).*ph).*Tsr)  +((1-hv2).*ph.*Tsh);
    dis2.Ts_v1 = Ts_v1;
    dis2.Ts_v2 = Ts_v2;
    
    dis2.g2_v1 = (1-(1-hv1).*ph)./Ts_v1;
    dis2.g2_v2 = (1-(1-hv2).*ph)./Ts_v2;
    dis2.h_v1  = (1-hv1).*ph./Ts_v1;
    dis2.h_v2 = (1-hv2).*ph./Ts_v2;
    
end