function [delta] = CG_tensorlines(t0,x0,vfield,T)
% Function used in minimizing shear components of CG-tensor and minimizing
% difference of the diagonal entries, used to find degenerate points of CG
d0 = zeros(2,5);
eps = 0.001;
X0 = aux_grid(x0,eps);
for k = 2:5
    CG = CG_tensor(t0,X0(k,:),vfield,T);
    d0(1,k) = CG(1,1) - CG(2,2);
    d0(2,k) = CG(1,2);
end
a = 1/2*(d0(1,3)-d0(1,2))/(2*eps);
b = 1/2*(d0(1,5)-d0(1,4))/(2*eps);
c = (d0(2,3)-d0(2,2))/(2*eps);
d = (d0(2,5)-d0(2,4))/(2*eps);
delta = a*d-b*c;
end