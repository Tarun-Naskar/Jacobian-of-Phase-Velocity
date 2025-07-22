function[dx2pt_1drpt,dx3pt_1drpt,dx4pt_1drpt,dx5pt_1drpt] = dxdrhon(XM,pt,Vs,rho,c,all_q1,all_q2,all_q3,all_q4,all_a,all_b,all_y1,all_y2,all_z1,all_z2) 
%% INPUT
% Layer Recursion Parameters related to half-space
%% OUTPUT
% Partial derivatives of propagator vector's elements w.r.t density of
% half-space.
%%
gmpt = Vs(pt)^2./c.^2;            
dEidr2 = 1/rho(pt-1) ;        detidr2 = -2*dEidr2*gmpt ;

a = all_a(:,pt-1) ;   a_ = a - 1;     b = all_b(:,pt-1);    b_ = b - 1;

daidri2 = dEidr2 + detidr2 ;   da_idri2 = daidri2;    dbidri2 = - detidr2;   db_idri2 = dbidri2;

q1 = all_q1(:,pt-1) ;
q2 = all_q2(:,pt-1) ;
q3 = all_q3(:,pt-1) ;
q4 = all_q4(:,pt-1) ;

y1 = all_y1(:,pt-1) ;
y2 = all_y2(:,pt-1) ;
dy1dri2 = XM(:,end-1,1).*da_idri2 + q1.*daidri2;
dy2dri2 = XM(:,end-1,1).*daidri2 + q2.*da_idri2;

dx2pt_1drpt = y1.*daidri2 + a.*dy1dri2 + y2.*da_idri2 + a_.*dy2dri2 ; 

dx3pt_1drpt = dEidr2.*q3 ;
dx4pt_1drpt = dEidr2.*q4;

z1 = all_z1(:,pt-1) ;
z2 = all_z2(:,pt-1) ;
dz1dri2 = XM(:,end-1,1).*dbidri2 + q1.*db_idri2;
dz2dri2 = XM(:,end-1,1).*db_idri2 + q2.*dbidri2;
dx5pt_1drpt = z1.*db_idri2 + b_.*dz1dri2 + z2.*dbidri2 + b.*dz2dri2 ; 
end