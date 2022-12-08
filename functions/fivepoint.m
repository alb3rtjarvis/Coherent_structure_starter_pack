function [df] = fivepoint(f,steps,h)
%%%---5 point differentiation formula for vectorized f---%%%
df = zeros(size(f));
for i = 1:steps
    if i == 1 || i == 2
        df(:,i) = (-3*f(:,i+4) + 16*f(:,i+3) - 36*f(:,i+2) + 48*f(:,i+1)-25*f(:,i))/(12*h);
    elseif i == steps || i == steps - 1
        df(:,i) = (3*f(:,i) - 4*f(:,i-1) + f(:,i-2))/(2*h);
    else
        df(:,i) = (-f(:,i+2) + 8*f(:,i+1) - 8*f(:,i-1) + f(:,i-2))/(12*h);
    end
end
end

