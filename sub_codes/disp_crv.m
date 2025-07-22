function [D] = disp_crv(f,mode,h,beta,rho,alpha) % disp_crv(f_min,df,f_max,mode,h,beta,rho,alpha)

% this function calculate dispersion curve using fast delta matrix algorithm
% f_min = minimum frequency ; df = frequency interval ; f_max = maximum frequency
% mode = mode of Raykeigh wave ; h = layer thickness of model; 
% beta = s-wave velocity ; alpha = p-wave velocity
%  rho = density

n = length(beta);

%%
c_min = 0.9*min(beta); % minimum Rayleigh wave velocity
c_max = max(beta);     % maximum Rayleigh wave velocity
dc = 0.5;              % velocity step
c = c_min : dc : c_max;
%% Dispersion function calculation
F_R = zeros(length(c), length(f));

for ii = 1 : length(f)
    freq = f(ii);
    F_R (:,ii)  = fast_delta_v2(n, h, beta, alpha, rho, freq, c); % fastest
    F_R(:,ii) = fillmissing(F_R(:,ii),'spline'); % fillmissing replace the 'NaN' by interpolating the other value
    x1 = find(F_R(2:end,ii) < 0 & F_R(1:end-1,ii)>0);
    x2 = find(F_R(1:end-1,ii) < 0 & F_R(2:end,ii)>0);
    x3 = [x1;x2];
    x = dc.* x3;
    x = sort(x); x = x + c_min;
    
    [x] = refine2(x, n, h, beta, alpha, rho, freq,dc);  
    
    D(ii,1:length(x))=x;
end

%% Ploting
D(D==0)=NaN; D = real(D);
D = D(:,1:mode);

end