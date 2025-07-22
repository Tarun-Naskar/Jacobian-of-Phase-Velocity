function[x,dx,cXM,all_p1,all_p2,all_p3,all_p4,all_q1,all_q2,all_q3,all_q4,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb] = dDdc(n,nc,x,t1,mu1,Vp, Vs, rho, c, k, d)  

% This function calculate partial derivative of dispersion function with
% respect to phase velocity
% inputs : Vp = p-wave velocity ; Vs = s-wave velocity ; rho = densiry ;
% c = phase velocity ; k = wave number ; d = layer thickness
%outputs : x = final state vector ; dx = derivatibe of x; cXM = a 3D matrix 
%containing all x calculated for all Rayleigh wave velocity(c) and all layers

dt1_dc = -2*c./Vs(1)^2;                      % derivative of t1 w.r.to phase velocity  
dx1_dc = 2*dt1_dc*mu1^2;                     % derivative of 1st element of x
dx2_dc = -2*t1.*dt1_dc*mu1^2;                % derivative of 2nd element of x

dx = zeros (nc,5);                            % derivative of state vector
dx(:,1) = dx1_dc; dx(:,2) = dx2_dc;
               

cXM = zeros(nc,n,5) ;
cXM(:,1,:) = x ; 


all_p1 = zeros(nc,n-1) ; all_p2 = zeros(nc,n-1) ; all_p3 = zeros(nc,n-1) ; all_p4 = zeros(nc,n-1) ; 
all_q1 = zeros(nc,n-1) ; all_q2 = zeros(nc,n-1) ; all_q3 = zeros(nc,n-1) ; all_q4 = zeros(nc,n-1) ; 
all_y1 = zeros(nc,n-1) ; all_y2 = zeros(nc,n-1) ;
all_z1 = zeros(nc,n-1) ; all_z2 = zeros(nc,n-1) ;
all_E = zeros(nc,n-1) ; all_a = zeros(nc,n-1) ; all_b = zeros(nc,n-1) ; 
all_r = zeros(nc,n-1) ; all_s = zeros(nc,n-1) ;
all_Ca = zeros(nc,n-1) ; all_Sa = zeros(nc,n-1) ;  all_Cb = zeros(nc,n-1) ; all_Sb = zeros(nc,n-1) ;

for m = 1:n-1
    
    E = rho(m+1)/rho(m);                    % a parameter relating density of two consecutive layers %(see Eqn A7)
    gm1 = Vs(m)^2./c.^2;                    % a parameter relating phase and s-wave velocity of m_th layer     %(see Eqn A2)       
    gm2 = Vs(m+1)^2 ./ c.^2;                % a parameter relating phase and s-wave velocity of (m+1)_th layer %(see Eqn A2)
    et = 2*(gm1 - E*gm2);                   %(see Eqn A7)

    dEdc = 0;                               % derivative of E with respect to phase velocity
    dgm_dc1 = -2*Vs(m)^2./c.^3;             % derivative of gm1 with respect to phase velocity
    dgm_dc2 = -2*Vs(m+1)^2./c.^3;           % derivative of gm2 with respect to phase velocity
    detdc = 2*(dgm_dc1 - E*dgm_dc2);        % derivative of et with respect to phase velocity

    a = E + et;
    a_ = a - 1;
    b = 1 - et;                             % a, a_, b, b_ are some defined parameter to simplify calculations %(see Eqn A9 & A10)
    b_ = b - 1;

    da_dc = dEdc + detdc;
    dadc = da_dc; 
    dbdc = -2*(dgm_dc1 - E*dgm_dc2);        % derivative of a, a_, b, b_ w.r.to phase velocity
    db_dc = dbdc; 

    r0 = 1 - c.^2./Vp(m)^2;                 % r0 is part eigen function(r) of 0_th interface %(see Eqn A4)
    s0 = 1 - c.^2./Vs(m)^2;                 % s0 is part eigen function(s) of 0_th interface %(see Eqn A6)
    dr0_dc = -2*c./Vp(m)^2;                 % derivative of r0 w.r.to phase velocity
    ds0_dc = -2*c./Vs(m)^2;                 % derivative of s0 w.r.to phase velocity

    r = sqrt(r0);                           % eigen function of 0_th interface
    Ca = cosh(k.*r.*d(m));                  %(see Eqn A4)
    Sa = sinh(k.*r.*d(m));                  %(see Eqn A4)

    dr_dc = 0.5*dr0_dc./sqrt(r0);           % derivative of r with respect to phase velocity
    dCa_dc = Sa.*k.*d(m).*(dr_dc-r./c);     % derivative of Ca with respect to phase velocity
    dSa_dc = Ca.*k.*d(m).*(dr_dc-r./c);     % derivative of Sa with respect to phase velocity

    s = sqrt(s0);                           % eigen function of 0_th interface
    Cb = cosh(k.*s.*d(m));                  %(see Eqn A6)
    Sb = sinh(k.*s.*d(m));                  %(see Eqn A6)

    ds_dc = 0.5*ds0_dc./sqrt(s0);           % derivative of s with respect to phase velocity
    dCb_dc = Sb.*k.*d(m).*(ds_dc-s./c);     % derivative of Cb with respect to phase velocity
    dSb_dc = Cb.*k.*d(m).*(ds_dc-s./c);     % derivative of Sb with respect to phase velocity

    p1 = Cb.*x(:,2) + s.*Sb.*x(:,3); p2 = Cb.*x(:,4) + s.*Sb.*x(:,5); p3 = (1./s).*Sb.*x(:,2) + Cb.*x(:,3);  p4 = (1./s).*Sb.*x(:,4) + Cb.*x(:,5); %(see Eqn A13 to A16)

    dp1dc = dCb_dc.*x(:,2) + Cb.*dx(:,2) + s.*Sb.*dx(:,3) + s.*dSb_dc.*x(:,3) + ds_dc.*Sb.*x(:,3); %(see Eqn B11)
    dp2dc = dCb_dc.*x(:,4) + Cb.*dx(:,4) + s.*Sb.*dx(:,5) + s.*dSb_dc.*x(:,5) + ds_dc.*Sb.*x(:,5); %(see Eqn B12)
    dp3dc = (1./s).*Sb.*dx(:,2) + (1./s).*dSb_dc.*x(:,2) - (1./s.^2).*ds_dc.*Sb.*x(:,2) + Cb.*dx(:,3) + dCb_dc.*x(:,3); %(see Eqn B13)
    dp4dc = (1./s).*Sb.*dx(:,4) + (1./s).*dSb_dc.*x(:,4) - (1./s.^2).*ds_dc.*Sb.*x(:,4) + Cb.*dx(:,5) + dCb_dc.*x(:,5); %(see Eqn B14)

    q1 = Ca.*p1 - r.*Sa.*p2; q2 = -(1./r).*Sa.*p3 + Ca.*p4;   q3 = Ca.*p3 - r.*Sa.*p4; q4 = -(1./r).*Sa.*p1 + Ca.*p2;   %(see Eqn A17 to A20)

    dq1dc = (Ca.*dp1dc + dCa_dc.*p1) - (r.*Sa.*dp2dc + r.*dSa_dc.*p2 + dr_dc.*Sa.*p2);                      %(see Eqn B15)
    dq2dc = -((1./r).*Sa.*dp3dc + (1./r).*dSa_dc.*p3 - (1./r.^2).*dr_dc.*Sa.*p3) + Ca.*dp4dc + dCa_dc.*p4;  %(see Eqn B16)
    dq3dc = (Ca.*dp3dc + dCa_dc.*p3) -(r.*Sa.*dp4dc + r.*dSa_dc.*p4 + dr_dc.*Sa.*p4);                       %(see Eqn B17)
    dq4dc = -((1./r).*Sa.*dp1dc + (1./r).*dSa_dc.*p1 - (1./r.^2).*dr_dc.*Sa.*p1) + Ca.*dp2dc + dCa_dc.*p2;  %(see Eqn B18)

    y1 = a_.*x(:,1) + a.*q1;    %(see Eqn A21)
    y2 = a.*x(:,1) + a_.*q2;    %(see Eqn A22)

    dy1dc = (x(:,1).*da_dc + a_.*dx(:,1)) + (a.*dq1dc + da_dc.*q1) ; 
    dy2dc = a.*dx(:,1) + dadc.*x(:,1) + a_.*dq2dc + dadc.*q2; 

    z1 = b.*x(:,1) + b_.*q1; 
    z2 = b_.*x(:,1) + b.*q2; 

    dz1dc = b.*dx(:,1) + dbdc.*x(:,1) + b_.*dq1dc + db_dc.*q1;  %(see Eqn B21)
    dz2dc = b_.*dx(:,1) + db_dc.*x(:,1) + b.*dq2dc + dbdc.*q2;  %(see Eqn B22)

    % update dx
    dx(:,1) = b_.*dy1dc + dbdc.*y1 + b.*dy2dc + dbdc.*y2;       %(see Eqn B23)
    dx(:,2) = (a.*dy1dc + dadc.*y1) + (a_.*dy2dc + da_dc.*y2);  %(see Eqn B24)
    dx(:,3) = E*dq3dc;                                          %(see Eqn B25)
    dx(:,4) = E*dq4dc;                                          %(see Eqn B26)
    dx(:,5) = b_.*dz1dc + db_dc.*z1 + b.*dz2dc + dbdc.*z2;      %(see Eqn B27)

    % update x
    x(:,1) = b_ .* y1 + b .* y2;    %(see Eqn A25)
    x(:,2) = a .* y1 + a_ .* y2;    %(see Eqn A26)
    x(:,3) = E * q3;                %(see Eqn A27)
    x(:,4) = E * q4;                %(see Eqn A28)
    x(:,5) = b_ .* z1 + b .* z2;    %(see Eqn A29)
    
    % store data for future use
    cXM(:,m+1,:) = x ; 
    
    all_p1(:,m) = p1 ;
    all_p2(:,m) = p2 ;
    all_p3(:,m) = p3 ;
    all_p4(:,m) = p4 ;
    
    all_q1(:,m) = q1 ;
    all_q2(:,m) = q2 ;
    all_q3(:,m) = q3 ;
    all_q4(:,m) = q4 ;
    
    all_y1(:,m) = y1 ;
    all_y2(:,m) = y2 ;
    
    all_z1(:,m) = z1 ;
    all_z2(:,m) = z2 ;
    
    all_E(:,m) = E ;
    all_a(:,m) = a ;
    all_b(:,m) = b ;
    
    all_r(:,m) = r ;
    all_s(:,m) = s ;
    
    all_Ca(:,m) = Ca ;
    all_Sa(:,m) = Sa ;
    all_Cb(:,m) = Cb ;
    all_Sb(:,m) = Sb ;
    
end

end