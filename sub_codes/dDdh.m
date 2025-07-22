function dD_dh = dDdh(n,k,pt,rL,sL,ccXM,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) 
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
% dD_dVp = Partial Derivative of Dispersion Function w.r.t. thickness
%%
[dx2idhi,dx3idhi,dx4idhi,dx5idhi] = dxdh(n,pt,k,ccXM,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) ; 

dD_dh = dx2idhi + sL.*dx3idhi - rL.*(dx4idhi + sL.*dx5idhi);

end
