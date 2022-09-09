function [ev] = CG_eigenvector(t0,x0,vfield,T,ev_flag)
% - Compute max/min eigenvectors of Cauchy-Green Tensor - %
if ~exist('ev_flag','var')
    ev_flag = 0;
end
delta = 1e-10;
CG = CG_tensor(t0,x0,vfield,T,delta); % compute rotation-independent deformation tensor
[V,D] = eig(CG);
if ev_flag == 0
    [~, idx_max] = max(diag(D));
    ev = V(:,idx_max);
else
    [~, idx_min] = min(diag(D));
    ev = V(:,idx_min);
end
end

