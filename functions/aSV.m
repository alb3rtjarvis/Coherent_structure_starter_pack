function [aSV_b,aSV_f] = aSV(t0,x0,vfield)
% - Compute approximation of asymptotic singular vectors - %
delta = 0.5;
tol = 1e-5;
U_b(:,:,1) = zeros(2);
V_f(:,:,1) = zeros(2);
svb_diff = 1;
svf_diff = 1;
k = 0;
while svb_diff > tol
    k = k+1;
    t1 = t0 - delta*k;
    dphi_b = Dphi(t1,x0,vfield,delta*k);
    [U_b(:,:,k+1),~,~] = svd(dphi_b);
    svb_diff = norm(U_b(:,:,k+1)-U_b(:,:,k));
    if k == 100
        fprintf('Timeout\n');        
        break
    end
end
aSV_b = U_b(:,:,k);
k = 0;
while svf_diff > tol
    k = k+1;
    t2 = t0 + delta*k;
    dphi_f = Dphi(t0,x0,vfield,delta*k);
    [~,~,V_f(:,:,k+1)] = svd(dphi_f);
    svf_diff = norm(V_f(:,:,k+1)-V_f(:,:,k));
    if k == 100
        fprintf('Timeout\n');        
        break
    end    
end
aSV_f = V_f(:,:,k);


