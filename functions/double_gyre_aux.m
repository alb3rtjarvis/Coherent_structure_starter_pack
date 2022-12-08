function [dxdt] = double_gyre_aux(t,x)
% - Double Gyre system of equations with phase shift and perturbation- %
A = 0.1;
eps = 0.1;
omega = 2*pi/10;
alpha = 0;
phi = 0;

f = (eps*sin(omega*t + phi))*(x(1))^2 + (1 - 2*eps*sin(omega*t + phi))*x(1);
dfdx = 2*eps*sin(omega*t + phi)*(x(1)) + (1 - 2*eps*sin(omega*t + phi));
dxdt = zeros(2,1);
dxdt(1) = -pi*A*sin(pi*f)*cos(pi*x(2)) - alpha*x(1);
dxdt(2) = pi*A*cos(pi*f)*sin(pi*x(2))*dfdx - alpha*x(2);
end

