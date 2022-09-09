function [F] = CG_opt(t0,x,vfield,T)
% Function used in minimizing shear components of CG-tensor and minimizing
% difference of the diagonal entries, used to find degenerate points of CG

CG = CG_tensor(t0,x,vfield,T);
F(1) = abs(CG(1,1) - CG(2,2));
F(2) = abs(CG(1,2));
end

