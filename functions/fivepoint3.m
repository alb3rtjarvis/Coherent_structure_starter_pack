function [df] = fivepoint3(f,steps,h)
%%%---5 point differentiation formula for 3rd derivative of f---%%%
df = zeros(size(f));
for i = 1:steps
    if i == 1 || i == 2 || i == 3
        df(:,i) = (-5*f(:,i) + 18*f(:,i+1) - 24*f(:,i+2) + 14*f(:,i+3) - 3*f(:,i+4))/(2*h^3);
    elseif i == steps || i == steps - 1 || i == steps - 2
        df(:,i) = (3*f(:,i-4) - 14*f(:,i-3) + 24*f(:,i-2) - 18*f(:,i-1) + 5*f(:,i))/(2*h^3);
    else
        df(:,i) = (f(:,i-3) - 8*f(:,i-2) + 13*f(:,i-1) - 13*f(:,i+1) + 8*f(:,i+2) - f(:,i+3))/(8*h^3);
    end
end
end

