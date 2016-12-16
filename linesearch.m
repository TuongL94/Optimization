function [lambda, No_of_iterations] = linesearch(func,x,d)
% Finds an approximate minimum of the function func along the line d using
% inexact line search with Armijo's rule.
% Input:   func - the objective function
%          x - the search point
%          d - the search direction
% Output:  lambda - an approximation of the minimum
%          No_of_iterations - number of iterations made by the line search
%          algorithm (always 1)

No_of_iterations = 1;
F_zero = func(x);
lambda = 1;
alpha = 2;
epsilon = 0.7;
h = 1/realmax; % Step size used for derivative approximation
fprim0 = (func(x+h*d) - func(x))/h; % Approximation of the derivative of func with forward difference
T = @(t) F_zero + epsilon*t*fprim0;

while func(x+lambda*d) >= T(lambda) || func(x+alpha*lambda*d) <= T(alpha*lambda)
    if func(x+lambda*d) > T(lambda)
        lambda = lambda/alpha;
    else
        lambda = lambda*alpha;
    end
end