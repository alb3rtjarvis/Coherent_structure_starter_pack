function CG = CG_tensor(t0,x0,vfield,T,delta)
% - Compute Cauchy-Green Tensor - %
if ~exist('delta','var')
    delta = 1e-7;
end
dphi = Dphi(t0,x0,vfield,T,delta);
CG = dphi'*dphi;
end