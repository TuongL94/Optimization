function [lambda, No_of_iterations] = armijo(func,x,d)

F_zero = func(x); % F(0) = f(x+0*d)
No_of_iterations = 1;
lambda = 1;
tol = 10^-50;
satisfied = false;

while satisfied ~= true
    F_lambda = func(x+lambda*d);
    if F_lambda == Inf
        lambda = lambda/10;
    elseif F_lambda > F_zero
        lambda = lambda/10;
    elseif F_lambda == F_zero 
        lambda = lambda*10;
    elseif isnan(F_lambda)
        lambda = lambda/2;
    else
        satisfied = true;
    end
end

