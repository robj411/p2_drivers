% function to estimate rates of change of infectious compartments and growth rate of infections
%
% dis2: struct of pathogen parameters
% compindex: struct of vectors of indices locating compartments
% y_mat: matrix of current disease states
%
% Iadot: rate of change of Ia 
% Isdot: rate of change of Is
% Iav1dot: rate of change of Iav1
% Isv1dot: rate of change of Isv1 
% Iav2dot: rate of change of Iav2
% Isv2dot: rate of change of Isv2
% IdotoverI: growth rate

function [Iadot, Isdot, Iav1dot, Isv1dot, Iav2dot, Isv2dot, IdotoverI] =  get_Idot(dis2, compindex, y_mat)

    E_index = compindex.E_index;
    I_index = compindex.I_index;

    E =      y_mat(:,E_index(1));
    Ev1 =      y_mat(:,E_index(2));
    Ev2 =      y_mat(:,E_index(3));
    
    Ia =    y_mat(:,I_index(1));
    Is =    y_mat(:,I_index(2));
    Iav1 =    y_mat(:,I_index(3));
    Isv1 =    y_mat(:,I_index(4));
    Iav2 =    y_mat(:,I_index(5));
    Isv2 =    y_mat(:,I_index(6));
    
    h    = dis2.h;
    g1   = dis2.g1;
    g2   = dis2.g2;
    g2_v1 = dis2.g2_v1;
    g2_v2 = dis2.g2_v2;
    h_v1  = dis2.h_v1;
    h_v2  = dis2.h_v2;
    sig1 = dis2.sig1;
    sig2 = dis2.sig2;

    Iadot=     sig1.*E     -g1.*Ia;
    Isdot=     sig2.*E     -(g2+h).*Is;
    Iav1dot=   sig1.*Ev1   -g1.*Iav1;
    Isv1dot=   sig2.*Ev1   -(g2_v1+h_v1).*Isv1;
    Iav2dot=   sig1.*Ev2   -g1.*Iav2;
    Isv2dot=   sig2.*Ev2   -(g2_v2+h_v2).*Isv2;
    
    I = sum(Ia + Is + Iav1 + Isv1 + Iav2 + Isv2);
    
    Idot = sum(Iadot + Isdot + Iav1dot + Isv1dot + Iav2dot + Isv2dot);
    
    IdotoverI = Idot/(I+1e-10);
    
end



