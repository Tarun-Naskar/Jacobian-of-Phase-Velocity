function[dx2pt_1dvpt,dx3pt_1dvpt,dx4pt_1dvpt,dx5pt_1dvpt] = half_space(XM,pt, Vs, c,all_q1,all_q2,all_a,all_b,all_y1,all_y2,all_z1,all_z2,all_E) 
%% INPUT
% Layer Recursion Parameters related to half-space
%% OUTPUT
% Partial derivatives of propagator vector's elements w.r.t S-wave velocity 
% of half-space.
%%
Ept_1 = all_E(:,pt-1) ;

a = all_a(:,pt-1) ;   a_ = a - 1;     b = all_b(:,pt-1);    b_ = b - 1;

daidvi2 = -4*Ept_1*Vs(pt)./c.^2;  da_idvi2 = daidvi2;  dbidvi2 = -daidvi2;  db_idvi2 = dbidvi2;     % refer Eqn 19 to 22 on paper

q1 = all_q1(:,pt-1) ;

q2 = all_q2(:,pt-1) ;

y1 = all_y1(:,pt-1) ;

y2 = all_y2(:,pt-1) ;
dy1pt_1dvpt = XM(:,end-1,1).*da_idvi2 + q1.*daidvi2 ;                           % refer Eqn 54a on paper
dy2pt_1dvpt = XM(:,end-1,1).*daidvi2 + q2.*da_idvi2 ;                           % refer Eqn 54b on paper
dx2pt_1dvpt = y1.*daidvi2 + a.*dy1pt_1dvpt + y2.*da_idvi2 + a_.*dy2pt_1dvpt ;   % refer Eqn 56b on paper

dx3pt_1dvpt = 0 ;       % refer Eqn 58a on paper

dx4pt_1dvpt = 0;        % refer Eqn 58b on paper

z1 = all_z1(:,pt-1) ;

z2 = all_z2(:,pt-1) ;
dz1pt_1dvpt = XM(:,end-1,1).*dbidvi2 + q1.*db_idvi2 ;                           
dz2pt_1dvpt = XM(:,end-1,1).*db_idvi2 + q2.*dbidvi2 ;                           
dx5pt_1dvpt = z1.*db_idvi2 + b_.*dz1pt_1dvpt + z2.*dbidvi2 + b.*dz2pt_1dvpt;    % refer Eqn 56c on paper
end