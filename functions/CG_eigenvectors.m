function [ev1,ev2] = CG_eigenvectors(t0,x0,vfield,T)
% - Compute max/min eigenvectors of Cauchy-Green Tensor - %
CG = CG_tensor(t0,x0,vfield,T); % compute rotation-independent deformation tensor
[V,D] = eig(CG);
[~, idx_max] = max(diag(D));
ev1 = V(:,idx_max);
[~, idx_min] = min(diag(D));
ev2 = V(:,idx_min);
end