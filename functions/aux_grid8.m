function X0 = aux_grid8(x0,delta)
% Create 8-point auxilary grid around X0 where distance between any two
% adjacent grid points is given by delta
X0(1,:) = x0;
X0(2,:) = [x0(1)-delta, x0(2)+delta];
X0(3,:) = [x0(1), x0(2)+delta];
X0(4,:) = [x0(1)+delta, x0(2)+delta];
X0(5,:) = [x0(1)+delta, x0(2)];
X0(6,:) = [x0(1)+delta, x0(2)-delta];
X0(7,:) = [x0(1), x0(2)-delta];
X0(8,:) = [x0(1)-delta, x0(2)-delta];
X0(9,:) = [x0(1)-delta, x0(2)];    
end