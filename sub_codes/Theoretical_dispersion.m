function [D] = Theoretical_dispersion(n,h,Vs,rho,Vp,f,mode,dc)
%% INPUT

% n = Number of finite layer
% Vs = S-wave velocity
% Vp = P-wave velocity
% rho = Density
% h = Layer thickness
% f = Frequency vector
% mode = Mode number
% dc = Velocity resolution
%% OUTPUT

% D = Dispersion curve
%% CREATE VELOCITY VECTOR
c_min = round(0.8*min(Vs));     % minimum Rayleigh wave velocity
c_max = 1.05*max (Vs);          % maximum Rayleigh wave velocity
c = c_min : dc : c_max;         % velocity vector
%% DISPERSION CURVE CALCULATION
F_R = zeros(length(c), length(f));
D = zeros (length(f), 15);
F_R  = fast_delta(n, h, Vs, Vp, rho, f, c);  % dispersion function
% Root Refinement
for ii = 1 : length(f)
    freq = f(ii);
    x1 = find(F_R(2:end,ii) < 0 & F_R(1:end-1,ii)>0);
    x2 = find(F_R(1:end-1,ii) < 0 & F_R(2:end,ii)>0);
    x3 = [x1;x2];
    x = dc.* x3;
    x = sort(x); x = x + c_min;
    [x] = refine(x, n, h, Vs, Vp, rho, freq, dc);  
    D(ii,1:length(x))=x;
end
D(D==0)=NaN; D = real(D);
D = D(:,1:mode);
end