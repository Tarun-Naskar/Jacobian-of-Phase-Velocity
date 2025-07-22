function [F_R] = fast_delta(n, h, Vs, Vp, rho, f, c)
%% INPUT 

% n = Number of finite layer
% Vs = S-wave velocity
% Vp = P-wave velocity
% rho = Density
% h = Layer thickness
% f = Frequency vector
% c = Velocity vector
%% OUTPUT

% F_R = Dispersion function (Buchen & Ben-Hador 1996)
%% CALCULATION STARTING
if length(Vs) ==1
    Vs = [Vs Vs];
    rho = [rho rho];
end
c = c';
f1 = repmat(f',length(c),1);
k = 2*pi*f1./c;

% Initial Value
t1 = (2-c.^2./Vs(1)^2);
mu1 = (rho(1)*Vs(1)^2);
X = zeros (length(t1),5);
X(:,1) = 2.*t1; X(:,2) = -t1.^2; X(:,5) = -4;
X = (mu1.^2) .* X;

X1 = repmat(X(:,1),1,length(f));
X2 = repmat(X(:,2),1,length(f));
X3 = repmat(X(:,3),1,length(f));
X4 = repmat(X(:,4),1,length(f));
X5 = repmat(X(:,5),1,length(f));

% Half-space Properties
r_l = sqrt(1-c.^2./Vp(end)^2);
s_l = sqrt(1-c.^2./Vs(end)^2);

% Layer Recursion
for m = 1 : n
    Vs1 = Vs(m);
    Vs2 = Vs(m+1);
    rho1 = rho(m);
    rho2 = rho(m+1);
    alpha1 = Vp(m);
    d = h(m);
    gamma1 = Vs1^2 ./ c.^2;
    gamma2 = Vs2^2 ./ c.^2;
    e = rho2 / rho1;
    et = 2 * (gamma1 - e * gamma2);

    a = e + et;
    a_ = a - 1;
    b = 1 - et;
    b_ = b - 1;

    r = sqrt(1 - (c.^2 ./ alpha1^2));
    C_a = cosh(k.*r.*d);
    S_a = sinh(k.*r.*d);

    s = sqrt(1 - (c.^2 ./ Vs1^2));
    C_B = cosh(k.*s.*d);
    S_B = sinh(k.*s.*d);

    p1 = C_B .* X2 + s .* S_B .* X3;
    p2 = C_B .* X4 + s .* S_B .* X5;

    m1 = (1./s) .* S_B; n1 = isnan(m1);
    m1(n1) = k(n1)*d;    

    p3 =  m1.* X2 + C_B .* X3; 
    p4 = m1 .* X4 + C_B .* X5; 

    m2 = (1./r) .* S_a; n2 = isnan(m2);
    m2(n2) = k(n2)*d;   

    q1 = C_a .* p1 - r .* S_a .* p2;
    q2 = -m2 .* p3 + C_a .* p4; 
    q3 = C_a .* p3 - r .* S_a .* p4;
    q4 = -m2 .* p1 + C_a .* p2; 

    y1 = a_ .* X1 + a .* q1;
    y2 = a .* X1 + a_ .* q2;

    z1 = b .* X1 + b_ .* q1;
    z2 = b_ .* X1 + b .* q2;

    X1 = b_ .* y1 + b .* y2;
    X2 = a .* y1 + a_ .* y2;
    X3 = e .* q3;
    X4 = e .* q4;
    X5 = b_ .* z1 + b .* z2;
end
% Dispersion Function
F_R = X2 + s_l .* X3 - r_l .* (X4 + s_l .* X5);
F_R = real(F_R);
end


