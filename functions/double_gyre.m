function [dxdt] = double_gyre(t,x)
% - Double Gyre system of equations - %
A = 0.1;
eps = 0.1;
omega = 2*pi/10;

f = (eps*sin(omega*t))*(x(1))^2 + (1 - 2*eps*sin(omega*t))*x(1);
dfdx = 2*eps*sin(omega*t)*(x(1)) + (1 - 2*eps*sin(omega*t));
dxdt = zeros(2,1);
dxdt(1) = -pi*A*sin(pi*f)*cos(pi*x(2));
dxdt(2) = pi*A*cos(pi*f)*sin(pi*x(2))*dfdx;
end

