function ftle = FTLE_field_fast(vfield,t,x,y,T)
% - Compute FTLE from equations given by vfield on grid given by x and y - %
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
ftlespan = @(t0,T) [t0 t0+T];
X_adv = zeros(ny,nx);
Y_adv = zeros(ny,nx);
k = 1;
for t0 = tspan
    starttime = tic;     
    for i = 1:ny
        for j = 1:nx
            FTLEspan = ftlespan(t0,T);
            [~, X_out] = ode45(vfield,FTLEspan,[X0(i,j),Y0(i,j)]);
            X_adv(i,j) = X_out(end,1); 
            Y_adv(i,j) = X_out(end,2);         
        end
    end
    phi = zeros(2);
    ftle_t = zeros(ny,nx);
    for i = 2:(ny-1)
        for j = 2:(nx-1)
            phi(1,1) = (X_adv(i+1,j)-X_adv(i-1,j))/(X0(i,j+1)-X0(i,j-1));
            phi(1,2) = (X_adv(i,j+1)-X_adv(i,j-1))/(Y0(i+1,j)-Y0(i-1,j));
            phi(2,1) = (Y_adv(i+1,j)-Y_adv(i-1,j))/(X0(i,j+1)-X0(i,j-1));
            phi(2,2) = (Y_adv(i,j+1)-Y_adv(i,j-1))/(Y0(i+1,j)-Y0(i-1,j));

            delta = phi'*phi; % compute rotation-independent deformation tensor
            [~,D] = eig(delta);
            lambda_max = max(diag(D)); % Find the maximum eigenvalue
            ftle_t(i,j)= 1/abs(T)*log(sqrt(lambda_max)); % the FTLE formula
        end
    end
    ftle(:,:,k) = ftle_t;
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

