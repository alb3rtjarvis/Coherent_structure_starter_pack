function [lambda_max,lambda_min,ev_max,ev_min] = CG_ev_field(vfield,t,x,y,T)
% - Compute max/min eigenvalues and eigenvectors of Cauchy-Green - %
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
lambda_max = zeros(ny,nx,nt);
lambda_min = zeros(ny,nx,nt);
ev_max = zeros(ny,nx,2,nt);
ev_min = zeros(ny,nx,2,nt);
k = 1;
for t1 = tspan
    starttime = tic;     
    for i = 2:(ny-1)
        for j = 2:(nx-1)
            CG = CG_tensor(t1,[X0(i,j),Y0(i,j)],vfield,T);
            [V,D] = eig(CG);
            [lambda_max(i,j,k), idx_max] = max(diag(D)); % Find the maximum eigenvalue
            [lambda_min(i,j,k), idx_min] = min(diag(D));
            ev_max(i,j,:,k) = V(:,idx_max);
            ev_min(i,j,:,k) = V(:,idx_min);
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
    k = k+1;
end
end

