function [lambda, No_of_iterations] = armijo2(func,x,d)
No_of_iterations = 1;
epsilon = 0.7;
alpha = 2;
F_0 = func(x);
F = @(t) func(x+t*d);
h = 1/realmax;
F_prim = (F(h)-F(0))/h;
if F_prim == 0 
    F_prim = h;
end
T = @(t) F_0 + epsilon*t*F_prim;


lambda = 1;

while (F(lambda) > T(lambda)) || (F(alpha*lambda) < T(alpha*lambda))
    if F(lambda) > T(lambda)
        lambda = lambda/alpha;
    else
        lambda = lambda*alpha;
    end
end

