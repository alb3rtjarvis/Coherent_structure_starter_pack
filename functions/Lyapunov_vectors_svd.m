function [LV_s,LV_u,LV_s_v,LV_u_v] = Lyapunov_vectors_svd(t0,x0,vfield,T,n)
% Script for finite time Lyapunov vectors using SVD %
intspan_s = linspace(t0+T,t0,n+1); % stable integration span
intspan_u = linspace(t0-T,t0,n+1);  % unstable integration span
delta = 1e-11;
opts = odeset('RelTol',1e-13,'AbsTol',1e-14);
[~, X_f] = ode45(vfield,[t0 t0+T],[x0(1),x0(2)],opts); % define point for pullback to act on to find stable LV 
[~, X_b] = ode45(vfield,[t0 t0-T],[x0(1),x0(2)],opts);  % define point for pushforward to act on to find unstable LV 
x0_f = X_f(end,:);           
x0_b = X_b(end,:);
Dphi_si = zeros(2,2,n+1);
Dphi_ui = zeros(2,2,n+1);
Dphi_si(:,:,1) = eye(2);
Dphi_ui(:,:,1) = eye(2);
for i = 1:n
    X0_s = aux_grid(x0_f,delta);
    x_adv_s = zeros(5,1);
    y_adv_s = zeros(5,1);    
    INTspan_s = [intspan_s(i), intspan_s(i+1)];
    
    X0_u = aux_grid(x0_b,delta);
    x_adv_u = zeros(5,1);
    y_adv_u = zeros(5,1);    
    INTspan_u = [intspan_u(i), intspan_u(i+1)];
    for k = 1:5
        [~, X_out_s] = ode45(vfield,INTspan_s,[X0_s(k,1),X0_s(k,2)],opts);
        x_adv_s(k) = X_out_s(end,1); 
        y_adv_s(k) = X_out_s(end,2);
        
        [~, X_out_u] = ode45(vfield,INTspan_u,[X0_u(k,1),X0_u(k,2)],opts);
        x_adv_u(k) = X_out_u(end,1); 
        y_adv_u(k) = X_out_u(end,2);        
    end
    
    Dphi_b = zeros(2);
    Dphi_b(1,1) = (x_adv_s(3)-x_adv_s(2))/(2*delta);
    Dphi_b(1,2) = (x_adv_s(5)-x_adv_s(4))/(2*delta);
    Dphi_b(2,1) = (y_adv_s(3)-y_adv_s(2))/(2*delta);
    Dphi_b(2,2) = (y_adv_s(5)-y_adv_s(4))/(2*delta);
    Dphi_si(:,:,i+1) = Dphi_b*Dphi_si(:,:,i);
    x0_f = [x_adv_s(1) y_adv_s(1)];
    
    Dphi_f = zeros(2);
    Dphi_f(1,1) = (x_adv_u(3)-x_adv_u(2))/(2*delta);
    Dphi_f(1,2) = (x_adv_u(5)-x_adv_u(4))/(2*delta);
    Dphi_f(2,1) = (y_adv_u(3)-y_adv_u(2))/(2*delta);
    Dphi_f(2,2) = (y_adv_u(5)-y_adv_u(4))/(2*delta);
    Dphi_ui(:,:,i+1) = Dphi_f*Dphi_ui(:,:,i);
    x0_b = [x_adv_u(1) y_adv_u(1)];
    
%     clear X0_f X0_b Dphi_f Dphi_b
end

Dphi_s = Dphi_si(:,:,end);
[Us,~,Vs] = svd(Dphi_s);

Dphi_u = Dphi_ui(:,:,end);
[Uu,~,Vu] = svd(Dphi_u);

LV_s = Us(:,1);
LV_u = Uu(:,1);

LV_s_v = Vs(:,1);
LV_u_v = Vu(:,1);
end