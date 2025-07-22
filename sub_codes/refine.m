%% Root Refinement using False-position Method
function [x5]=refine(x3, n, h, Vs, Vp, rho, freq, dc) 
x5 = zeros(1,length(x3));
for ii=1:length(x3) 
    x4 = x3(ii);
    a = x4 - dc;
    b = x4 + dc;

    for jj=1:1:8
        f_l = fast_delta(n, h, Vs, Vp, rho, freq, a);
        f_u = fast_delta(n, h, Vs, Vp, rho, freq, b);
        
        new = b - f_u * (b - a)/(f_u - f_l);
        f_new = fast_delta(n, h, Vs, Vp, rho, freq, new);
        
        if (f_new*f_l)<0
            b=new;
        else
            a=new;
        end 
    end 
    x5(ii) = new;
end
end
    