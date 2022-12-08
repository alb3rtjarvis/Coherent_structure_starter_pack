function [LV_s,LV_u,LE_s,LE_u] = Lyapunov_vector(t0,x0,vfield,T,n)
% Function for finite time approximation of leading Lyapunov vector in  
% forward and backwards time using gram schmidt orthogonalization in 2D
% - t0 -- initial time
% - x0 -- initial point
% - vfield -- vector field defined by a function
% - T -- finite integration time (must be positive)
% - n -- number of steps in re-orthogonaliztion process

% Defines full integration span
intspan_s = linspace(t0+T,t0,n+1); % stable integration span
intspan_u = linspace(t0-T,t0,n+1);  % unstable integration span

delta = 1e-11;                       % spacing between auxillary grid points
E_s(:,:,1) = eye(2);               % initialize vectors for LV search
E_u(:,:,1) = eye(2);
opts = odeset('RelTol',1e-6,'AbsTol',1e-7);
lambda_s = zeros(2,n);
lambda_u = zeros(2,n);
[~, X_s] = ode45(vfield,[t0 t0+T],[x0(1),x0(2)],opts); % define point for pullback to act on to find stable LV 
[~, X_u] = ode45(vfield,[t0 t0-T],[x0(1),x0(2)],opts);  % define point for pushforward to act on to find unstable LV 
x0_s = X_s(end,:);           
x0_u = X_u(end,:);
for i = 1:n
    X0_s = aux_grid(x0_s,delta);  % define auxillary grid around point
    x_adv_s = zeros(5,1);
    y_adv_s = zeros(5,1);    
    INTspan_s = [intspan_s(i), intspan_s(i+1)]; % define intspan for each step
    
    X0_u = aux_grid(x0_u,delta);        % same as above for stable
    x_adv_u = zeros(5,1);
    y_adv_u = zeros(5,1);    
    INTspan_u = [intspan_u(i), intspan_u(i+1)];    
    for k = 1:5
        % Integration used to define Dphi for stable LV
        [~, X_out_s] = ode45(vfield,INTspan_s,[X0_s(k,1),X0_s(k,2)],opts);
        x_adv_s(k) = X_out_s(end,1); 
        y_adv_s(k) = X_out_s(end,2);
        
        % Integration used to define Dphi for unstable LV
        [~, X_out_u] = ode45(vfield,INTspan_u,[X0_u(k,1),X0_u(k,2)],opts);
        x_adv_u(k) = X_out_u(end,1); 
        y_adv_u(k) = X_out_u(end,2);        
    end
    % Pullback
    Dphi_s = zeros(2);
    Dphi_s(1,1) = (x_adv_s(3)-x_adv_s(2))/(2*delta);
    Dphi_s(1,2) = (x_adv_s(5)-x_adv_s(4))/(2*delta);
    Dphi_s(2,1) = (y_adv_s(3)-y_adv_s(2))/(2*delta);
    Dphi_s(2,2) = (y_adv_s(5)-y_adv_s(4))/(2*delta);
    
    % Pushforward
    Dphi_u = zeros(2);
    Dphi_u(1,1) = (x_adv_u(3)-x_adv_u(2))/(2*delta);
    Dphi_u(1,2) = (x_adv_u(5)-x_adv_u(4))/(2*delta);
    Dphi_u(2,1) = (y_adv_u(3)-y_adv_u(2))/(2*delta);
    Dphi_u(2,2) = (y_adv_u(5)-y_adv_u(4))/(2*delta);    
    
    % Modified Gram-Schmidt to approximate stable LV
    [E_s(:,:,i+1),R_s] = mgs(Dphi_s*E_s(:,:,i));
    lambda_s(:,i) = log(diag(R_s));   % log of norm of vector at each step 
    x0_s = [x_adv_s(1) y_adv_s(1)];  % defines new point to be used in next step
    
    % Modified Gram-Schmidt to approximate unstable LV
    [E_u(:,:,i+1),R_u] = mgs(Dphi_u*E_u(:,:,i));
    lambda_u(:,i) = log(diag(R_u));     % log of norm of vector at each step 
    x0_u = [x_adv_u(1) y_adv_u(1)];     % defines new point to be used in next step
    
    clear X0_s X0_u
end
%[~,idx_s] = max(sum(lambda_s,2)/n);   % approximation of max LE from pullback starting at x0_s
LV_s = E_s(:,1,end); % corresponding LV
LE_s = sum(lambda_s(1,:),2)/n;

%[~,idx_u] = max(sum(lambda_u,2)/n);     % approximation of max LE from pushforward starting at x0_u
LV_u = E_u(:,1,end); % corresponding LV
LE_u = sum(lambda_u(1,:),2)/n; 
end


