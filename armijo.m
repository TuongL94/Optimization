function [lambda, No_of_iterations] = armijo(func,x,d)

F_zero = func(x);
No_of_iterations = 1;
lambda = 1;
alpha = 2;
epsilon = 0.7;
h = 1/realmax;

fprim0 = (func(x+h*d) - func(x))/h;
T = @(t) F_zero + epsilon*t*fprim0;

while func(x+lambda*d) >= T(lambda) || func(x+alpha*lambda*d) <= T(alpha*lambda)
    if func(x+lambda*d) > T(lambda)
        lambda = lambda/alpha;
    else
        lambda = lambda*alpha;
    end
end

