function X0 = aux_grid(x0,delta)
% Create 4-point auxilary grid around x0 point with spacing delta
    X0(1,:) = x0;
    X0(2,:) = [x0(1)-delta x0(2)];
    X0(3,:) = [x0(1)+delta x0(2)];
    X0(4,:) = [x0(1) x0(2)-delta];
    X0(5,:) = [x0(1) x0(2)+delta];
end