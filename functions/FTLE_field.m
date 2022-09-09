function [ftle] = FTLE_field(vfield,t0,x,y,T)
% - Compute FTLE field from equations given by vfield - %
if size(x,1) == 1 || size(x,2) == 1
    nx = length(x);
    ny = length(y);
    [X0,Y0] = meshgrid(x,y);
else
    nx = size(x,2);
    ny = size(x,1);
    X0 = x;
    clear x;
    Y0 = y;
    clear y;
end
ftle = zeros(ny,nx);
for i = 1:ny
    for j = 1:nx
        x0 = [X0(i,j),Y0(i,j)];
        CG = CG_tensor(t0,x0,vfield,T,1e-2);
        lambda_max = max(eig(CG));
        ftle(i,j) = 1/(2*abs(T))*log(lambda_max);        
    end
end
end

