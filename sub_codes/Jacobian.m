function [f,dVr_dVs,dVr_dVp,dVr_dVrho,dVr_dVh] = Jacobian(Vs, Vp_poi, rho, h, mode, f_min, f_max, df, dc) 
%% INPUT
% Vs = S-wave velocity
% Vp_poi = P-wave velocity or Poisson's ratio
% rho = Density
% h = Layer thickness
% mode = Mode number
% f_min = Minimum frequency
% f_max = Maximum frequency
% df = Frequency resolution
% dc = Velocity resolution
%% OUTPUT
% f = Frequency vector
% dVr_dVs = Jacobian w.r.t. S-wave velocity
% dVr_dVp = Jacobian w.r.t. P-wave velocity
% dVr_dVrho = Jacobian w.r.t. density
% dVr_dVh = Jacobian w.r.t. thickness
%%
f = f_min:df:f_max;             % frequency vector
%% if Vp_poi is Poisson's ratio, convert Vp_poi to P-wave velocity
if Vp_poi(1)<0.5
    multiplier = sqrt((1-Vp_poi)./(0.5 - Vp_poi)) ;
    Vp_poi = Vs.*multiplier ;
else
    Vp = Vp_poi;
end
%% NAVIGATE TO THE FUNCTION CODES
addpath('sub_codes');
%% DISPERSION CURVE
Vr = Theoretical_dispersion(length(Vs)-1,h,Vs,rho,Vp,f',mode,dc);
Vr = Vr(:,mode) ;
%% CALCULATION OF SUPPORTIVE PARAMETERS
tic
nc = length(Vr);                            
n = length(Vs);

rL0 = 1 - (Vr.^2)./(Vp(end)^2);      % rL0 is a parameter relating phase and P-wave velocity 
sL0 = 1 - (Vr.^2)./(Vs(end)^2);      % sL0 is a parameter relating phase and S-wave velocity 
drL0_dVr = -(2*Vr)./(Vp(end)^2);     % derivative of rL0 with respect to phase velocity 
dsL0_dVr = -(2*Vr)./(Vs(end)^2);     % derivative of sL0 with respect to phase velocity

rL = sqrt(rL0);                      % eigen function (refer Eqn 3 on paper)
sL = sqrt(sL0);                      % eigen function (refer Eqn 3 on paper)

drL_dVr = 0.5.*drL0_dVr./sqrt(rL0);  % derivative of rL with respect to phase velocity 
dsL_dVr = 0.5.*dsL0_dVr./sqrt(sL0);  % derivative of sL with respect to phase velocity

dsL0_dVs = 2.*(Vr.^2)./(Vs(end)^3);  % derivative of rL0 with respect to S-wave velocity
dsL_dVs = 0.5.*dsL0_dVs./sL;         % derivative of rL with respect to S-wave velocity

drL0_dVp = 2.*(Vr.^2)./Vp(end)^3;    % derivative of rL0 with respect to P-wave velocity
drL_dVp = 0.5.*drL0_dVp./rL;         % derivative of rL with respect to P-wave velocity

%% JACOBIAN CALCULATION

mu1 = (rho(1)*Vs(1)^2);              % shear modulus of 1st layer 
t1 = 2 - (Vr).^2/Vs(1)^2;            % a parameter relating phase and s-wave velocity of 1st layer
x = zeros (length(t1),5);
x(:,1) = 2.*t1; x(:,2) = -t1.^2; x(:,5) = -4;   %refer Eqn 7 on paper
x = (mu1^2)*x;  

dVr_dVs = zeros(nc,n);
dVr_dVp = zeros(nc,n);
dVr_dVh = zeros(nc,n-1);
dVr_dVrho = zeros(nc,n);

k = 2*pi*f'./Vr;                     % wave number

% State Vector and its Derivative w.r.to Phase Velocity
[X,dx,cXM,all_p1,all_p2,all_p3,all_p4,all_q1,all_q2,all_q3,all_q4,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb] = dDdc(n,nc,x,t1,mu1,Vp, Vs, rho, Vr, k, h);  

% Partial Derivative of Dispersion Function(D) w.r.t Phase Velocity
dD_dVr = dx(:,2) + sL.*dx(:,3) + dsL_dVr.*X(:,3) - (rL.*(dx(:,4) + sL.*dx(:,5) + dsL_dVr.*X(:,5)) + drL_dVr.*(X(:,4) + sL.*X(:,5))); 

for pt = 1:n                         % n = number of finite layer + half-space

        ccXM = cXM(:,1:pt,:) ;
    
        dD_dVs = dDdVs(n,mu1,t1,k,Vs,rho,h,Vr,pt,dsL_dVs,rL,sL,ccXM,all_q1,all_q2,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb);        % derivative of dispersion function w.r.to S-wave velocity
        dVr_dVs(:,pt) = - dD_dVs./dD_dVr;     % Jacobian of phase velocity w.r.to S-wave velocity    
        
        dD_dVp = dDdVp(n,k,Vp,h,Vr,pt,drL_dVp,rL,sL,ccXM,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) ;                                % derivative of dispersion function w.r.to P-wave velocity
        dVr_dVp(:,pt) = - dD_dVp./dD_dVr;     % Jacobian of phase velocity w.r.to P-wave velocity
         
        dD_drho = dDdrho(n,mu1,t1,Vs,rho,Vr,pt,rL,sL,ccXM,all_q1,all_q2,all_q3,all_q4,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) ;   % derivative of dispersion function w.r.to density
        dVr_dVrho(:,pt) = - dD_drho./dD_dVr;  % Jacobian of phase velocity w.r.to density
        
        if pt<n
        dD_dh = dDdh(n,k,pt,rL,sL,ccXM,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) ;                                                  % derivative of dispersion function w.r.to layer thickness
        dVr_dVh(:,pt) = - dD_dh./dD_dVr;      % Jacobian of phase velocity w.r.to layer thickness
        end 
end
toc
end
