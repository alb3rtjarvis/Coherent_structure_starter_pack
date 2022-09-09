function [ftle] = FTLE_data(vel,t,x,y,T,fixed_timestep_flag,plotflag)
% - Compute FTLE from stacked velocity data vel - %
if ~exist('fixed_timestep_flag','var')
    fixed_timestep_flag = 0;
end
if ~exist('plotflag','var')
    plotflag = 0;
end
if size(x,1) == 1 || size(x,2) == 1
    nx = length(x);
    ny = length(y);
    dx = abs(x(2)-x(1));
    dy = abs(y(2)-y(1));
    [X0,Y0] = meshgrid(x,y);
else
    nx = size(x,2);
    ny = size(x,1);
    dx = abs(x(1,2) - x(1,1));
    dy = abs(y(2,1) - y(1,1));
    X0 = x;
    clear x;
    Y0 = y;
    clear y;
end
for k = 1:size(vel,2)
    XX(:,:,k) = X0;
    YY(:,:,k) = Y0;
    u(:,:,k) = reshape(vel(1:ny*nx,k),[ny nx]);
    v(:,:,k) = reshape(vel(ny*nx+1:end,k),[ny nx]);
    TT(:,:,k) = t(k)*ones(ny,nx);
end
P = [2 1 3];    
Xprime = permute(XX,P);
Yprime = permute(YY,P);
Tprime = permute(TT,P);
Uprime = permute(u,P);
Vprime = permute(v,P);

u_interp = griddedInterpolant(Xprime,Yprime,Tprime,Uprime,'linear');
v_interp = griddedInterpolant(Xprime,Yprime,Tprime,Vprime,'linear');   

function vel_interp = velocity_interp_func(tt,X_in)
    xx = X_in(1);
    yy = X_in(2);

    if xx < X0(1,1) 
        xx = X0(1,1);
    end
    if xx > X0(1,end)
        xx = X0(1,end);
    end
    if yy < Y0(1,1) 
        yy = Y0(1,1);
    end
    if yy > Y0(end,1)
        yy = Y0(end,1);
    end
    if tt < t(1)
        tt = t(1);
    end
    if tt > t(end)
        tt = t(end);
    end
    uu = u_interp(xx,yy,tt);
    vv = v_interp(xx,yy,tt);
    vel_interp = [uu;vv];
end
dt = t(2)-t(1);
tspan = t(1):dt:t(end)-T;
nt = length(tspan);
ftle = zeros(ny,nx,nt);
ftle_t = zeros(ny,nx);
if fixed_timestep_flag == 0
    intspan = @(t0,T) [t0 t0+T];
else
    intspan = @(t0,T) linspace(t0,t0+T,51);
end
for k = 1:nt
    starttime = tic;    
    for i = 1:ny
        for j = 1:nx
            [t, X_out] = ode45(@velocity_interp_func,intspan(tspan(k),T),[X0(i,j),Y0(i,j)]);
            X_adv(i,j) = X_out(length(FTLEspan),1); 
            Y_adv(i,j) = X_out(length(FTLEspan),2);         
        end
        disp(i);
    end
    phi = zeros(2);
    ftle_t = zeros(ny,nx);
    for i = 2:(ny-1)
        for j = 2:(nx-1)
            phi(1,1) = (X_adv(i+1,j)-X_adv(i-1,j))/(X0(i,j+1)-X0(i,j-1));
            phi(1,2) = (X_adv(i,j+1)-X_adv(i,j-1))/(Y0(i+1,j)-Y0(i-1,j));
            phi(2,1) = (Y_adv(i+1,j)-Y_adv(i-1,j))/(X0(i,j+1)-X0(i,j-1));
            phi(2,2) = (Y_adv(i,j+1)-Y_adv(i,j-1))/(Y0(i+1,j)-Y0(i-1,j));
 
            C = phi'*phi; % compute rotation-independent deformation tensor 
            lambda_max = max(eig(C)); % Find the maximum eigenvalue  
            ftle_t(i,j)= 1/abs(T)*log(sqrt(lambda_max)); % the FTLE formula
        end
    end
    if plotflag == 1
        [~,ch] = contourf(X0,Y0,ftle_t,35);
        set(ch,'edgecolor','none'); % remove the lines between contour levels
        daspect([1 1 1]); % set the aspect ratio to 1:1
        drawnow
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
end
end

