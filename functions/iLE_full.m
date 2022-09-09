function [s1,s2,s1_v,s2_v] = iLE_full(vfield,t,x,y)
% - Function to compute iLE from vector field vfield over tspan t - %
if size(x,1) == 1 || size(x,2) == 1
    nx = length(x);
    ny = length(y);
    dx = abs(x(2)-x(1));
    dy = abs(y(2)-y(1));
    [X,Y] = meshgrid(x,y);
else
    nx = size(x,2);
    ny = size(x,1);
    dx = abs(x(1,2) - x(1,1));
    dy = abs(y(2,1) - y(1,1));
    X = x;
    Y = y;
    clear x y
end
nt = length(t);
vel = zeros(ny,nx,2);
s1 = zeros(ny,nx,nt);
s2 = zeros(ny,nx,nt);
s1_v = zeros(ny,nx,2,nt);
s2_v = zeros(ny,nx,2,nt);
for k = 1:nt
    for i = 1:ny
        for j = 1:nx
            [vel(i,j,:)] = vfield(t(k),[X(i,j) Y(i,j)]);
        end
    end
    u = vel(:,:,1); v = vel(:,:,2);
    [dudx, dudy] = gradient(u,dx,dy);
    [dvdx, dvdy] = gradient(v,dx,dy);
    for i = 1:ny
        for j = 1:nx
            grad = [dudx(i,j) dudy(i,j); dvdx(i,j) dvdy(i,j)];
            S = 0.5*(grad + grad');
            [V,D] = eig(S);
            [s1(i,j,k), idx_max] = max(diag(D)); % Find the maximum eigenvalue
            [s2(i,j,k), idx_min] = min(diag(D));
            s1_v(i,j,:,k) = V(:,idx_max);
            s2_v(i,j,:,k) = V(:,idx_min);                  
        end
    end
end
end
