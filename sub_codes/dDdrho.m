function dD_drho = dDdrho(n,mu1,t1,V_s,rho,c,pt,rL,sL,ccXM,all_q1,all_q2,all_q3,all_q4,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) 
%% INPUT
% k = Wave number ; 
% V_s = S-wave velocity ; 
% V_p = P-wave velocity;
% rho = Density ; 
% h = Thickness ; 
% c = Phase velocity ; 
% pt = density changed layer ;
% rL, sL = Eigen functions of half-space
%% OUTPUT
% dD_dVp = Partial Derivative of Dispersion Function w.r.t. density
%%
if pt == n    % if pt is half-space

    [dx2idvi,dx3idvi,dx4idvi,dx5idvi] = dxdrhon(ccXM,pt,V_s,rho,c,all_q1,all_q2,all_q3,all_q4,all_a,all_b,all_y1,all_y2,all_z1,all_z2); 

    dD_drho = dx2idvi + sL.*dx3idvi - rL.*(dx4idvi + sL.*dx5idvi);
else

    [dx2idvi,dx3idvi,dx4idvi,dx5idvi] = dxdrho(n,mu1,t1,pt, V_s, rho, c, ccXM,all_q1,all_q2,all_q3,all_q4,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb); 
    dD_drho = dx2idvi + sL.*dx3idvi - rL.*(dx4idvi + sL.*dx5idvi);
end