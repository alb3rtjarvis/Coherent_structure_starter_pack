function [df] = fivepoint2(f,steps,h)
%%%---5 point differentiation formula for 2nd derivative of f---%%%
df = zeros(size(f));
for i = 1:steps
    if i == 1 || i == 2 
        df(:,i) = (-f(:,i+3) + 4*f(:,i+2) - 5*f(:,i+1) + 2*f(:,i))/(h^2);
    elseif i == steps || i == steps - 1
        df(:,i) = (2*f(:,i) - 5*f(:,i-1) + 4*f(:,i-2) - f(:,i-3))/(h^2);
    else
        df(:,i) = (-f(:,i+2) + 16*f(:,i+1) - 30*f(:,i) + 16*f(:,i-1) - f(:,i-2))/(12*h^2);
    end
end
end

