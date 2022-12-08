function [ftle] = FTLE(t0,x0,vfield,T)
% Compute FTLE at a point x0 at time t0 in vfield for integration time T
    CG = CG_tensor(t0,x0,vfield,T);
    lambda_max = max(eig(CG));
    ftle = 1/(2*abs(T))*log(lambda_max);
end
