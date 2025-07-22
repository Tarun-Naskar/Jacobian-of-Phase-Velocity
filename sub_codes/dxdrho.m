function[dx2idri,dx3idri,dx4idri,dx5idri] = dxdrho(n,mu1,t1,pt,Vs, rho, c, XM,all_q1,all_q2,all_q3,all_q4,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb)
%% INPUT
% Layer Recursion Parameters Related to a Finite Layer
%% OUTPUT
% Partial derivatives of propagator vector's elements w.r.t density of a 
% finite layer
   
for i = pt:n-1
    if i == pt
        if pt==1
            dmu1dr1 = Vs(1)^2;
            dx1pt_1drpt = 2*(2*mu1*dmu1dr1*t1);
            dx2pt_1drpt = -(2*mu1*dmu1dr1*t1.^2);
            dx3pt_1drpt = 0;
            dx4pt_1drpt = 0;
            dx5pt_1drpt = -8*mu1*dmu1dr1;
        else
            
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

            dx1pt_1drpt = y1.*db_idri2 + b_.*dy1dri2 + y2.*dbidri2 + b.*dy2dri2 ;
            dx2pt_1drpt = y1.*daidri2 + a.*dy1dri2 + y2.*da_idri2 + a_.*dy2dri2 ; 
            
            dx3pt_1drpt = dEidr2*q3 ;
            dx4pt_1drpt = dEidr2*q4;
            
            z1 = all_z1(:,pt-1) ;
            z2 = all_z2(:,pt-1) ;
            
            dz1dri2 = XM(:,end-1,1).*dbidri2 + q1.*db_idri2;
            dz2dri2 = XM(:,end-1,1).*db_idri2 + q2.*dbidri2;
            dx5pt_1drpt = z1.*db_idri2 + b_.*dz1dri2 + z2.*dbidri2 + b.*dz2dri2 ; 
        end
        
        gmpt_p1 = Vs(pt+1)^2./c.^2; Ept = all_E(:,pt);
        
        dEidr = - rho(pt+1)/rho(pt)^2 ;        detidr = -2*dEidr*gmpt_p1 ;

        r = all_r(:,pt) ; s = all_s(:,pt) ;

        a = all_a(:,pt) ;   a_ = a - 1;     b = all_b(:,pt);    b_ = b - 1;

        daidri = dEidr + detidr ;      da_idri = daidri;        dbidri = - detidr;         db_idri = dbidri;

        Ca = all_Ca(:,pt) ;    Sa = all_Sa(:,pt) ;

        Cb = all_Cb(:,pt) ;    Sb = all_Sb(:,pt) ;

        q1 = all_q1(:,pt) ;
        y1 = all_y1(:,pt) ;
        
        dp1ptdrpt = Cb.*dx2pt_1drpt + s.*Sb.*dx3pt_1drpt ;
        dp2ptdrpt = Cb.*dx4pt_1drpt + s.*Sb.*dx5pt_1drpt ;
        dq1ptdrpt = Ca.*dp1ptdrpt - r.*Sa.*dp2ptdrpt;
        dy1ptdrpt = a_.*dx1pt_1drpt + XM(:,end,1).*da_idri + q1.*daidri + a.*dq1ptdrpt;

        q2 = all_q2(:,pt) ;

        y2 = all_y2(:,pt) ;
        
        dp3ptdrpt = (1./s).*Sb.*dx2pt_1drpt + Cb.*dx3pt_1drpt ;
        dp4ptdrpt = (1./s).*Sb.*dx4pt_1drpt + Cb.*dx5pt_1drpt ;
        dq2ptdrpt = (-1./r).*Sa.*dp3ptdrpt + Ca.*dp4ptdrpt;
        dy2ptdrpt = a.*dx1pt_1drpt + XM(:,end,1).*daidri + q2.*da_idri + a_.*dq2ptdrpt;

        dx1ptdrpt = y1.*db_idri + b_.*dy1ptdrpt  + y2.*dbidri + b.*dy2ptdrpt ;
        dx1idri = dx1ptdrpt; 


        dx2ptdrpt = y1.*daidri + a.*dy1ptdrpt + y2.*da_idri + a_.*dy2ptdrpt ;
        dx2idri = dx2ptdrpt; 

        q3 = all_q3(:,pt) ;
        dq3ptdrpt = Ca.*dp3ptdrpt - r.*Sa.*dp4ptdrpt ;
        dx3ptdrpt = Ept.*dq3ptdrpt + dEidr*q3 ;
        dx3idri = dx3ptdrpt; 
        
        q4 = all_q4(:,pt) ;
        dq4ptdrpt = -(1./r).*Sa.*dp1ptdrpt + Ca.*dp2ptdrpt ;
        dx4ptdrpt = Ept.*dq4ptdrpt + dEidr*q4;
        dx4idri = dx4ptdrpt; 
        
        z1 = all_z1(:,pt) ;
        dz1ptdrpt = b.*dx1pt_1drpt + XM(:,end,1).*dbidri + q1.*db_idri+ b_.*dq1ptdrpt ;

        z2 = all_z2(:,pt) ;
        dz2ptdrpt = b_.*dx1pt_1drpt + XM(:,end,1).*db_idri + q2.*dbidri + b.*dq2ptdrpt ;
        
        dx5ptdrpt = z1.*db_idri + b_.*dz1ptdrpt + z2.*dbidri + b.*dz2ptdrpt ; 
        dx5idri = dx5ptdrpt; 
    else

        E = all_E(:,i);

        r = all_r(:,i) ; s = all_s(:,i) ;

        a = all_a(:,i) ;   a_ = a - 1;     b = all_b(:,i);    b_ = b - 1;

        Ca = all_Ca(:,i) ;    Sa = all_Sa(:,i) ;

        Cb = all_Cb(:,i) ;    Sb = all_Sb(:,i) ;

        
        p1fdrpt = Cb.*dx2idri + s.*Sb.*dx3idri ;
        p2fdrpt = Cb.*dx4idri + s.*Sb.*dx5idri ;
        q1fdrpt = Ca.*p1fdrpt - r.*Sa.*p2fdrpt ;
        y1fdrpt = a_.*dx1idri + a.*q1fdrpt ;
        
        p3fdrpt = (1./s).*Sb.*dx2idri + Cb.*dx3idri ;
        p4fdrpt = (1./s).*Sb.*dx4idri + Cb.*dx5idri ;
        q2fdrpt = (-1./r).*Sa.*p3fdrpt + Ca.*p4fdrpt ;
        y2fdrpt = a.*dx1idri + a_.*q2fdrpt ;
        
        x1fdrpt = b_.*y1fdrpt + b.*y2fdrpt;

        x2fdrpt = a.*y1fdrpt + a_.*y2fdrpt;

        q3fdrpt = Ca.*p3fdrpt - r.*Sa.*p4fdrpt ;
        x3fdrpt = E.*q3fdrpt;

        q4fdrpt = -(1./r).*Sa.*p1fdrpt + Ca.*p2fdrpt ;
        x4fdrpt = E.*q4fdrpt;

        z1fdrpt = b.*dx1idri + b_.*q1fdrpt ;
        z2fdrpt = b_.*dx1idri + b.*q2fdrpt ;
        x5fdrpt = b_.*z1fdrpt + b.*z2fdrpt ;

        dx1idri = x1fdrpt; dx2idri = x2fdrpt; dx3idri = x3fdrpt;  dx4idri = x4fdrpt;  dx5idri = x5fdrpt;
    end

end
end
