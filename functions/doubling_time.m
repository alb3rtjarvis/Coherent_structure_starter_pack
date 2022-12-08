function [t_d] = doubling_time(t0,x0,vfield)
% Function used to compute doubling time from a point x0 at t0 in vfield
sigma_1 = 0;
nk = 1000;
t = linspace(1e-5,100,nk);
k = 0;
while sigma_1 < 2 && k < nk
    k = k+1;
    dphi = Dphi(t0,x0,vfield,t(k));
    [~,S,~] = svd(dphi);
    sigma_1 = S(1,1);
end

t_d = t(k);

