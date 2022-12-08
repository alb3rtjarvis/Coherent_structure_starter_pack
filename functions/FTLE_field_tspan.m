function [ftle] = FTLE_field_tspan(vfield,t,x,y,T)
% - Compute FTLE field over time tspan from equations given by vfield - %
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
dt = t(2)-t(1);
tspan = t(1):dt:t(end)-T;
nt = length(tspan);
ftle = zeros(ny,nx,nt);
for k = 1:nt
    starttime = tic;     
    for i = 1:ny
        for j = 1:nx
            x0 = [X0(i,j),Y0(i,j)];
            CG = CG_tensor(tspan(k),x0,vfield,T);
            lambda_max = max(eig(CG));
            ftle(i,j,k) = 1/(2*abs(T))*log(lambda_max);       
        end
    end
    endtime = toc(starttime);
    timeleft = (nt-k)*endtime;
    if timeleft > 3600
        disptext = [num2str(timeleft/3600),' hours remaining'];
    elseif timeleft < 60
        disptext = [num2str(timeleft),' seconds remaining'];
    else
        disptext = [num2str(timeleft/60),' minutes remaining'];
    end
    disp(disptext);
    disp(k);
end
end

