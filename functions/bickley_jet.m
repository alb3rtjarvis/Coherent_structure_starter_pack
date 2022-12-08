function [dxdt] = bickley_jet(t,x)
%%% Function defining system of equations representing bickley jet
U0=62.66; %m/s
%U0 = 3.6*U0; %km/hr
U0 = U0/24; %km/day
L = 1770.0; %kms
% L = L/4000; % no dim
A1 = 0.075;
A2 = 0.4;
A3 = 0.3;
c2 = 0.205*U0;
c3 = 0.461*U0;
re = 6371.0; %earths radius, kms
% re = re/4000; %no dim
k1 = 2/re;
k2 = 4/re;
k3 = 6/re;
dxdt = zeros(4,1);

dxdt(1) = x(3);
dxdt(2) = x(4);
dxdt(3) = U0*((1/cosh(x(2))/L)^2) + 2*A3*U0*tanh(x(2)/L)*((1/cosh(x(2)/L))^2)*cos(k3*(x(1))-c3*t) +...
          2*A2*U0*tanh(x(2)/L)*((1/cosh(x(2)/L))^2)*cos(k2*(x(1)-c2*t));
dxdt(4) = -k3*A3*L*U0*((1/cosh(x(2)/L))^2)*sin(k3*(x(1)-c3*t)) - k2*A2*L*U0*((1/cosh(x(2)/L))^2)*sin(k2*(x(1)-c2*t));  
end

