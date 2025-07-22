function dD_dVs = dDdVs(n,mu1,t1,k,V_s,rho,h,c,pt,dsL_dVs,rL,sL,ccXM,all_q1,all_q2,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) 
%% INPUT
% k = Wave number ; 
% V_s = S-wave velocity ; 
% V_p = P-wave velocity;
% rho = Density ; 
% h = Thickness ; 
% c = Phase velocity ; 
% pt = V_s changed layer ;
% dsL_dVs = Derivative of sL w.r.t V_s ; 
% rL, sL = Eigen functions of half-space
%% OUTPUT
% dD_dVs = Partial Derivative of Dispersion Function w.r.t. S-wave Velocity
%%
if pt == n  % if pt is half-space

    [dx2idvi,dx3idvi,dx4idvi,dx5idvi] = half_space(ccXM,pt, V_s,c,all_q1,all_q2,all_a,all_b,all_y1,all_y2,all_z1,all_z2,all_E); 
 

    dD_dVs = dx2idvi + sL.*dx3idvi + ccXM(:,end,3).*dsL_dVs - rL.*(dx4idvi + sL.*dx5idvi + ccXM(:,end,5).*dsL_dVs);  % derivative of dispersion function(D) with respect to half-space's S-wave velocity(Vs)

else

    [dx2idvi,dx3idvi,dx4idvi,dx5idvi] = dX_dVs(n,mu1,t1,pt, h, V_s, rho, k, c,ccXM ,all_q1,all_q2,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb); 
    
    dD_dVs = dx2idvi + sL.*dx3idvi - rL.*(dx4idvi + sL.*dx5idvi);   % derivative of dispersion function(D) with respect to a finite layer's S-wave velocity(Vs) 
end
