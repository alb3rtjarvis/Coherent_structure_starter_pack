function F = CG_search(t0,x,vfield,T)
% Function used in minimizing shear components of CG-tensor and minimizing
% difference of the diagonal entries, used to find degenerate points of CG

CG = CG_tensor(t0,x,vfield,T);
a = abs(CG(1,1) - CG(2,2));
b = abs(CG(1,2));
F = a+b;
end

