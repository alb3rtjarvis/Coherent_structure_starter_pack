function [s1,s2,e1,e2] = iLE_vel_full(vel,t,x,y)
% - Function to compute iLE from velocity array with x,y components stacked
% - in vectors
if size(x,1) == 1 || size(x,2) == 1
    nx = length(x);
    ny = length(y);
    dx = abs(x(2)-x(1));
    dy = abs(y(2)-y(1));
else
    nx = size(x,2);
    ny = size(x,1);
    dx = abs(x(1,2) - x(1,1));
    dy = abs(y(2,1) - y(1,1));
end
nt = length(t);
s1 = zeros(ny,nx);
s2 = zeros(ny,nx);
for k = 1:nt
    u = reshape(vel(1:nx*ny,k),[ny nx]);
    v = reshape(vel(nx*ny+1:end,k),[ny nx]);
    [dudx, dudy] = gradient(u,dx,dy);
    [dvdx, dvdy] = gradient(v,dx,dy);
    for i = 1:ny
        for j = 1:nx
            grad = [dudx(i,j) dudy(i,j); dvdx(i,j) dvdy(i,j)];
            S = 0.5*(grad + grad');
            [V,D] = eig(S);
            [s1(i,j,k), idx_min] = min(diag(D));
            [s2(i,j,k), idx_max] = max(diag(D));
            e1(i,j,:,k) = V(:,idx_min);
            e2(i,j,:,k) = V(:,idx_max);
        end
    end
end
end
