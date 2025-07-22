function[dx2idhi,dx3idhi,dx4idhi,dx5idhi] = dxdh(n,pt,k,XM,all_p1,all_p2,all_p3,all_p4,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb)
%% INPUT
% Layer Recursion Parameters Related to a Finite Layer
%% OUTPUT
% Partial derivatives of propagator vector's elements w.r.t thickness of a 
% finite layer
   
for i = pt:n-1
    if i == pt
        Ept = all_E(:,pt);

        r = all_r(:,pt);    s = all_s(:,pt);
        
        a = all_a(:,pt) ;   a_ = a - 1;     b = all_b(:,pt);    b_ = b - 1;

        Ca = all_Ca(:,pt);     Sa = all_Sa(:,pt);
        dCadh = Sa.*k.*r;      dSadh = Ca.*k.*r;

        Cb = all_Cb(:,pt);     Sb = all_Sb(:,pt);
        dCbdh = Sb.*k.*s;      dSbdh = Cb.*k.*s;

        p1 = all_p1(:,pt) ;

        p2 = all_p2(:,pt) ;
        
        dp1ptdvpt = dCbdh.*XM(:,end,2) + s.*dSbdh.*XM(:,end,3) ;
        dp2ptdvpt = dCbdh.*XM(:,end,4) + s.*dSbdh.*XM(:,end,5) ;
        dq1ptdvpt = dCadh.*p1 + Ca.*dp1ptdvpt - r.*(Sa.*dp2ptdvpt + dSadh.*p2) ;
        dy1ptdvpt = a.*dq1ptdvpt ;

        p3 = all_p3(:,pt) ;

        p4 = all_p4(:,pt) ;
        
        dp3ptdvpt = (1./s).*dSbdh.*XM(:,end,2) +  dCbdh .* XM(:,end,3) ;
        dp4ptdvpt = (1./s).*dSbdh.*XM(:,end,4) +  dCbdh .* XM(:,end,5) ;
        dq2ptdvpt = (-1./r).*(Sa.*dp3ptdvpt + dSadh.*p3) + Ca.*dp4ptdvpt + dCadh.*p4 ;
        dy2ptdvpt = a_.*dq2ptdvpt ;

        dx1ptdvpt = b_.*dy1ptdvpt  + b.*dy2ptdvpt ;
        dx1idhi = dx1ptdvpt; 


        dx2ptdvpt = a.*dy1ptdvpt + a_.*dy2ptdvpt ;
        dx2idhi = dx2ptdvpt; 

        dq3ptdvpt = Ca.*dp3ptdvpt + dCadh.*p3 - r.*(Sa.*dp4ptdvpt + dSadh.*p4) ;
        dx3ptdvpt = Ept.*dq3ptdvpt ;
        dx3idhi = dx3ptdvpt; 
        
        dq4ptdvpt = (-1./r).*(Sa.*dp1ptdvpt + dSadh.*p1) + Ca.*dp2ptdvpt + dCadh.*p2 ;
        dx4ptdvpt = Ept.*dq4ptdvpt ;
        dx4idhi = dx4ptdvpt; 

        dz1ptdvpt = b_.*dq1ptdvpt ;
        dz2ptdvpt = b.*dq2ptdvpt ;
        
        dx5ptdvpt = b_.*dz1ptdvpt + b.*dz2ptdvpt ; 
        dx5idhi = dx5ptdvpt; 
    else
        E = all_E(:,i);

        r = all_r(:,i);  s = all_s(:,i);

        a = all_a(:,i) ;   a_ = a - 1;     b = all_b(:,i);    b_ = b - 1;

        Ca = all_Ca(:,i);     Sa = all_Sa(:,i);

        Cb = all_Cb(:,i);     Sb = all_Sb(:,i);

        
        p1fdvpt = Cb.*dx2idhi + s.*Sb.*dx3idhi ;
        p2fdvpt = Cb.*dx4idhi + s.*Sb.*dx5idhi ;
        q1fdvpt = Ca.*p1fdvpt - r.*Sa.*p2fdvpt ;
        y1fdvpt = a_.*dx1idhi + a.*q1fdvpt ;
        
        p3fdvpt = (1./s).*Sb.*dx2idhi + Cb.*dx3idhi ;
        p4fdvpt = (1./s).*Sb.*dx4idhi + Cb.*dx5idhi ;
        q2fdvpt = (-1./r).*Sa.*p3fdvpt + Ca.*p4fdvpt ;
        y2fdvpt = a.*dx1idhi + a_.*q2fdvpt ;
        
        x1fdvpt = b_.*y1fdvpt + b.*y2fdvpt;

        x2fdvpt = a.*y1fdvpt + a_.*y2fdvpt;

        q3fdvpt = Ca.*p3fdvpt - r.*Sa.*p4fdvpt ;
        x3fdvpt = E.*q3fdvpt;

        q4fdvpt = -(1./r).*Sa.*p1fdvpt + Ca.*p2fdvpt ;
        x4fdvpt = E.*q4fdvpt;

        z1fdvpt = b.*dx1idhi + b_.*q1fdvpt ;
        z2fdvpt = b_.*dx1idhi + b.*q2fdvpt ;
        x5fdvpt = b_.*z1fdvpt + b.*z2fdvpt ;

        dx1idhi = x1fdvpt; dx2idhi = x2fdvpt; dx3idhi = x3fdvpt;  dx4idhi = x4fdvpt;  dx5idhi = x5fdvpt;
    end

end
end