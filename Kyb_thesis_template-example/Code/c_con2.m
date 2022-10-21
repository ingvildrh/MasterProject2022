function [c, ceq] = c_con2(z)
global alpha beta lambda_t mx N lambda_dot_t
    j = 1;
    for i = 1:mx:(N-1)*mx
        c(j) = alpha*exp(-beta*(z(i)-lambda_t)^2) - z(i+4);
        j = j + 1;
    end
    ceq = [];
    
    % Add constraint on travel rate
    for i = 1:mx:(N-1)*mx
        c(j) = abs(z(i+1)) - lambda_dot_t;
        j = j + 1;
    end
end