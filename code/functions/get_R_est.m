% function to estimate Rt using growth rate of infections and durations in
% compartments
%
% dis2: struct of pathogen parameters
% compindex: struct of vectors of indices locating compartments
% y_mat: matrix of current disease states
% p3: fraction of asymptomatic infectiousness averted
% p4: fraction of asymptomatic infectiousness averted
%
% R_est: estimated value of Rt

function R_est =  get_R_est(dis2, compindex, y_mat, p3, p4)
    
    I_index = compindex.I_index;
    Ia =    y_mat(:,I_index(1));
    Is =    y_mat(:,I_index(2));
    Iav1 =    y_mat(:,I_index(3));
    Isv1 =    y_mat(:,I_index(4));
    Iav2 =    y_mat(:,I_index(5));
    Isv2 =    y_mat(:,I_index(6));
                
    [~,~,~,~,~,~,growth_rate] = get_Idot(dis2, compindex, y_mat);
    b1 = 1/dis2.Tlat;
    
    Ia_tot = sum(Ia + Iav1 + Iav2)*dis2.red*(1-p3);
    Is_tot = sum(Is + Isv1 + Isv2)*(1-p4);
    Ia_dur = Ia_tot/dis2.Tay;
    Is_dur = sum(Is./dis2.Ts + Isv1./dis2.Ts_v1 + Isv2./dis2.Ts_v2)*(1-p4);
    i_period = (Ia_dur + Is_dur) ;
    b2 = i_period/(Ia_tot+Is_tot+1e-10) + 1e-10;
    
    R_est = (1+growth_rate/b1)*(1+growth_rate/b2);
    
end



