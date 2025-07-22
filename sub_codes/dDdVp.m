function dD_dVp = dDdVp(n,k,V_p,h,c,pt,drL_dVp,rL,sL,ccXM,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb)  % ,all_p1,all_p2,all_p3,all_p4
%% INPUT
% k = Wave number ; 
% V_s = S-wave velocity ; 
% V_p = P-wave velocity;
% rho = Density ; 
% h = Thickness ; 
% c = Phase velocity ; 
% pt = V_p changed layer ;
% rL, sL = Eigen functions of half-space
% drL_dVp = Derivative of rL w.r.t drL_dVp ; 
%% OUTPUT
% dD_dVp = Partial Derivative of Dispersion Function w.r.t. P-wave Velocity
%%
if pt == n   % if pt is half-space

    dD_dVp = - drL_dVp.*(ccXM(:,end,4) + sL.*ccXM(:,end,5));
else
    [dx2idvi,dx3idvi,dx4idvi,dx5idvi] = dxdvp(n,pt, h, V_p, k, c, all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) ; 
    dD_dVp = dx2idvi + sL.*dx3idvi - rL.*(dx4idvi + sL.*dx5idvi);
end
end
