function[dx2idvi,dx3idvi,dx4idvi,dx5idvi] = dX_dVs(n,mu1,t1,pt, d, Vs, rho, k, c,XM,all_q1,all_q2,all_y1,all_y2,all_z1,all_z2,all_E,all_a,all_b,all_r,all_s,all_Ca,all_Sa,all_Cb,all_Sb) 
%% INPUT
% Layer Recursion Parameters Related to a Finite Layer
%% OUTPUT
% Partial derivatives of propagator vector's elements w.r.t S-wave velocity 
% of a Finite Layer
%%
for i = pt:n-1
    if i == pt
        if pt==1
            dmu1dv1 = 2*rho(1)*Vs(1);                                   % refer Eqn 15 on paper
            dt1dv1 = 2*c.^2./Vs(1)^3;                                   % refer Eqn 8 on paper
            dx1pt_1dvpt = 2*(2*mu1*dmu1dv1*t1 + dt1dv1*mu1^2);          % refer Eqn 29a on paper
            dx2pt_1dvpt = -(2*t1.*dt1dv1*mu1^2 + 2*mu1*dmu1dv1*t1.^2);  % refer Eqn 29b on paper
            dx3pt_1dvpt = 0;                                            % refer Eqn 29c on paper
            dx4pt_1dvpt = 0;                                            % refer Eqn 29d on paper
            dx5pt_1dvpt = -8*mu1*dmu1dv1;                               % refer Eqn 29e on paper
        else
            Ept_1 = all_E(:,pt-1);
           
            a = all_a(:,pt-1) ;   a_ = a - 1;     b = all_b(:,pt-1);    b_ = b - 1;

            daidvi2 = -4*Ept_1*Vs(pt)./c.^2; da_idvi2 = daidvi2; dbidvi2 = -daidvi2; db_idvi2 = dbidvi2; 

            q1 = all_q1(:,pt-1) ;
            q2 = all_q2(:,pt-1) ;

            y1 = all_y1(:,pt-1) ;

            y2 = all_y2(:,pt-1) ;
            dy1dvi2 = (XM(:,end-1,1).*da_idvi2 + q1.*daidvi2);              % refer Eqn 54a on paper
            dy2dvi2 = (XM(:,end-1,1).*daidvi2 + q2.*da_idvi2);              % refer Eqn 54b on paper

            dx1pt_1dvpt = y1.*db_idvi2 + b_.*dy1dvi2 + y2.*dbidvi2 + b.*dy2dvi2 ;   % refer Eqn 56a on paper
            dx2pt_1dvpt = y1.*daidvi2 + a.*dy1dvi2 + y2.*da_idvi2 + a_.*dy2dvi2 ;   % refer Eqn 56b on paper
            dx3pt_1dvpt = 0 ;                                                       
            dx4pt_1dvpt = 0;                                                        

            z1 = all_z1(:,pt-1) ;

            z2 = all_z2(:,pt-1) ;
            dz1dvi2 = (XM(:,end-1,1).*dbidvi2 + q1.*db_idvi2);                      
            dz2dvi2 = (XM(:,end-1,1).*db_idvi2 + q2.*dbidvi2);                      
            dx5pt_1dvpt = z1.*db_idvi2 + b_.*dz1dvi2 + z2.*dbidvi2 + b.*dz2dvi2 ;   % refer Eqn 56c on paper
        end

        Ept = all_E(:,pt);

        ds0_db = 2*c.^2./Vs(pt)^3; 

        r = all_r(:,pt) ; s = all_s(:,pt) ;

        dsdv = 0.5*ds0_db./s;                                                       % refer Eqn 24 on paper 
        
        a = all_a(:,pt) ;   a_ = a - 1;     b = all_b(:,pt);    b_ = b - 1;

        daidvi = 4*Vs(pt)./c.^2; da_idvi = daidvi; dbidvi = -daidvi; db_idvi = dbidvi;  % refer Eqn 19-22 on paper

        Ca = all_Ca(:,pt) ;    Sa = all_Sa(:,pt) ;

        Cb = all_Cb(:,pt) ;    Sb = all_Sb(:,pt) ;

        dCbdv = Sb.*k.*d(pt).*dsdv; dSbdv = Cb.*k.*d(pt).*dsdv;     % refer Eqn 27- and 28 on paper

        q1 = all_q1(:,pt) ;
        y1 = all_y1(:,pt) ;
        
        dp1ptdvpt = Cb.*dx2pt_1dvpt + XM(:,end,2).*dCbdv + s.*Sb.*dx3pt_1dvpt + XM(:,end,3).*s.*dSbdv + XM(:,end,3).*Sb.*dsdv ;     % refer Eqn 50a on paper
        dp2ptdvpt = Cb.*dx4pt_1dvpt + XM(:,end,4).*dCbdv + s.*Sb.*dx5pt_1dvpt + XM(:,end,5).*s.*dSbdv + XM(:,end,5).*Sb.*dsdv ;     % refer Eqn 50b on paper
        dq1ptdvpt = Ca.*dp1ptdvpt - r.*Sa.*dp2ptdvpt;                                                                               % refer Eqn 52a on paper
        dy1ptdvpt = a_.*dx1pt_1dvpt + XM(:,end,1).*da_idvi + q1.*daidvi + a.*dq1ptdvpt;                                             % refer Eqn 54a on paper

        q2 = all_q2(:,pt) ;

        y2 = all_y2(:,pt) ;
        
        dp3ptdvpt = (1./s).*Sb.*dx2pt_1dvpt + XM(:,end,2).*(1./s).*dSbdv - XM(:,end,2).*Sb.*(1./s.^2).*dsdv + Cb.*dx3pt_1dvpt + XM(:,end,3).*dCbdv ;    % refer Eqn 50c on paper
        dp4ptdvpt = (1./s).*Sb.*dx4pt_1dvpt + XM(:,end,4).*(1./s).*dSbdv - XM(:,end,4).*Sb.*(1./s.^2).*dsdv + Cb.*dx5pt_1dvpt + XM(:,end,5).*dCbdv ;    % refer Eqn 50d on paper
        dq2ptdvpt = (-1./r).*Sa.*dp3ptdvpt + Ca.*dp4ptdvpt;                                                                                             % refer Eqn 52b on paper
        dy2ptdvpt = a.*dx1pt_1dvpt + XM(:,end,1).*daidvi + q2.*da_idvi + a_.*dq2ptdvpt;                                                                 % refer Eqn 54b on paper

        dx1ptdvpt = y1.*db_idvi + b_.*dy1ptdvpt  + y2.*dbidvi + b.*dy2ptdvpt ;          % refer Eqn 56a on paper
        dx1idvi = dx1ptdvpt; 


        dx2ptdvpt = y1.*daidvi + a.*dy1ptdvpt + y2.*da_idvi + a_.*dy2ptdvpt ;           % refer Eqn 56b on paper
        dx2idvi = dx2ptdvpt; 

        dq3ptdvpt = Ca.*dp3ptdvpt - r.*Sa.*dp4ptdvpt ;      % refer Eqn 52c on paper
        dx3ptdvpt = Ept.*dq3ptdvpt ;                        % refer Eqn 58a on paper
        dx3idvi = dx3ptdvpt; 
        
        dq4ptdvpt = -1./r.*Sa.*dp1ptdvpt + Ca.*dp2ptdvpt ;      % refer Eqn 52a on paper
        dx4ptdvpt = Ept.*dq4ptdvpt ;                            % refer Eqn 58b on paper
        dx4idvi = dx4ptdvpt; 

        z1 = all_z1(:,pt) ;
        dz1ptdvpt = b.*dx1pt_1dvpt + XM(:,end,1).*dbidvi + q1.*db_idvi+ b_.*dq1ptdvpt ;    

        z2 = all_z2(:,pt) ;
        dz2ptdvpt = b_.*dx1pt_1dvpt + XM(:,end,1).*db_idvi + q2.*dbidvi + b.*dq2ptdvpt ;    
        
        dx5ptdvpt = z1.*db_idvi + b_.*dz1ptdvpt + z2.*dbidvi + b.*dz2ptdvpt ;               % refer Eqn 56c on paper
        dx5idvi = dx5ptdvpt; 
    else
        % This 'else' portion is to complete the recursive process to reach
        % to half-space
        E = all_E(:,i);

        r = all_r(:,i) ; s = all_s(:,i) ;

        a = all_a(:,i) ;   a_ = a - 1;     b = all_b(:,i);    b_ = b - 1;

        Ca = all_Ca(:,i) ;    Sa = all_Sa(:,i) ;

        Cb = all_Cb(:,i) ;    Sb = all_Sb(:,i) ;
        
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

