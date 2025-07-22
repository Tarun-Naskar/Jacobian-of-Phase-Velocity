function[dx2idvi,dx3idvi,dx4idvi,dx5idvi] = dxdvp(n,pt, d,Vp, k, c,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb)
%% INPUT
% Layer Recursion Parameters Related to a Finite Layer
%% OUTPUT
% Partial derivatives of propagator vector's elements w.r.t P-wave velocity 
% of a Finite Layer
%%
for i = pt:n-1
    if i == pt
        Ept = all_E(:,pt);

        dr0_da = 2*c.^2./Vp(pt)^3;

        r = all_r(:,pt);  

        drdv = 0.5*dr0_da./r;

        a = all_a(:,pt) ;   a_ = a - 1;     b = all_b(:,pt);    b_ = b - 1;

        Ca = all_Ca(:,pt);     Sa = all_Sa(:,pt);
        dCadv = Sa.*k*d(pt).*drdv; dSadv = Ca.*k*d(pt).*drdv;

        p1 = all_p1(:,pt) ;
        p2 = all_p2(:,pt) ;

        dq1ptdvpt = p1.*dCadv - p2.*(r.*dSadv + drdv.*Sa) ;
        dy1ptdvpt = a.*dq1ptdvpt ;

        p3 = all_p3(:,pt) ;

        p4 = all_p4(:,pt) ;

        
        dq2ptdvpt = -p3.*(1./r .* dSadv - (1./r.^2).*drdv.*Sa) + p4.*dCadv;
        dy2ptdvpt = a_.*dq2ptdvpt ;

        dx1ptdvpt = b_.*dy1ptdvpt  + b.*dy2ptdvpt ;
        dx1idvi = dx1ptdvpt; 


        dx2ptdvpt = a.*dy1ptdvpt + a_.*dy2ptdvpt ;
        dx2idvi = dx2ptdvpt; 

        dq3ptdvpt = p3.*dCadv - p4.*(r.*dSadv + drdv.*Sa) ;
        dx3ptdvpt = Ept.*dq3ptdvpt ;
        dx3idvi = dx3ptdvpt; 
        
        dq4ptdvpt = -p1.*(1./r .* dSadv - (1./r.^2).*drdv.*Sa) + p2.*dCadv ;
        dx4ptdvpt = Ept.*dq4ptdvpt ;
        dx4idvi = dx4ptdvpt; 



        dz1ptdvpt = b_.*dq1ptdvpt ;

        dz2ptdvpt = b.*dq2ptdvpt ;
        
        dx5ptdvpt = b_.*dz1ptdvpt + b.*dz2ptdvpt ; 
        dx5idvi = dx5ptdvpt; 
    else
        E = all_E(:,i);

        r = all_r(:,i);  s = all_s(:,i);

        a = all_a(:,i) ;   a_ = a - 1;     b = all_b(:,i);    b_ = b - 1;

        Ca = all_Ca(:,i);     Sa = all_Sa(:,i);

        Cb = all_Cb(:,i);     Sb = all_Sb(:,i);

        
        p1fdvpt = Cb.*dx2idvi + s.*Sb.*dx3idvi ;
        p2fdvpt = Cb.*dx4idvi + s.*Sb.*dx5idvi ;
        q1fdvpt = Ca.*p1fdvpt - r.*Sa.*p2fdvpt ;
        y1fdvpt = a_.*dx1idvi + a.*q1fdvpt ;
        
        p3fdvpt = (1./s).*Sb.*dx2idvi + Cb.*dx3idvi ;
        p4fdvpt = (1./s).*Sb.*dx4idvi + Cb.*dx5idvi ;
        q2fdvpt = (-1./r).*Sa.*p3fdvpt + Ca.*p4fdvpt ;
        y2fdvpt = a.*dx1idvi + a_.*q2fdvpt ;
        
        x1fdvpt = b_.*y1fdvpt + b.*y2fdvpt;

        x2fdvpt = a.*y1fdvpt + a_.*y2fdvpt;

        q3fdvpt = Ca.*p3fdvpt - r.*Sa.*p4fdvpt ;
        x3fdvpt = E.*q3fdvpt;

        q4fdvpt = -(1./r).*Sa.*p1fdvpt + Ca.*p2fdvpt ;
        x4fdvpt = E.*q4fdvpt;

        z1fdvpt = b.*dx1idvi + b_.*q1fdvpt ;
        z2fdvpt = b_.*dx1idvi + b.*q2fdvpt ;
        x5fdvpt = b_.*z1fdvpt + b.*z2fdvpt ;

        dx1idvi = x1fdvpt; dx2idvi = x2fdvpt; dx3idvi = x3fdvpt;  dx4idvi = x4fdvpt;  dx5idvi = x5fdvpt;
    end

end
end
