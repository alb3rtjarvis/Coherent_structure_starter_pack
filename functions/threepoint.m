function [df] = threepoint(f,steps,h)
%%%---3 point differentiation formula for vectorized f---%%%
df = zeros(size(f));
for i = 1:steps
    if i == 1
        df(:,i) = (-f(:,i+2) + 4*f(:,i+1) - 3*f(:,i))/(2*h);
    elseif i == steps
        df(:,i) = (3*f(:,i) - 4*f(:,i-1) + f(:,i-2))/(2*h);
    else
        df(:,i) = (f(:,i+1) - f(:,i-1))/(2*h);
    end
end
end

