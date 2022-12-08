function dphi = Dphi(t0,x0,vfield,T,delta)
% - Compute derivative of flow map - %
if ~exist('delta','var')
    delta = 1e-3;
end
X0 = aux_grid(x0,delta);
ftlespan = @(t0,T) [t0 t0+T];   % Integration span
x_adv = zeros(5,1);
y_adv = zeros(5,1);
opts = odeset('RelTol',1e-7,'AbsTol',1e-8);    
FTLEspan = ftlespan(t0,T);
for k = 2:5
    [~, X_out] = ode45(vfield,FTLEspan,[X0(k,1),X0(k,2)],opts);
    x_adv(k) = X_out(end,1); 
    y_adv(k) = X_out(end,2);    
end

% Creation of dphi using center differencing
dphi = zeros(2);
dphi(1,1) = (x_adv(3)-x_adv(2))/(2*delta);
dphi(1,2) = (x_adv(5)-x_adv(4))/(2*delta);
dphi(2,1) = (y_adv(3)-y_adv(2))/(2*delta);
dphi(2,2) = (y_adv(5)-y_adv(4))/(2*delta);
end