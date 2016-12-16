function [lambda,No_of_iterations]=linesearch(func,x,d)
% Finds an interval of uncertainty using the Golden Section search method
% and returns the left end point of this interval as an approximation of
% the minimum of func along the line d.
% Input:   func - the objective function
%          x - the search point
%          d - the search direction
% Output:  lambda - left end point of the final interval of uncertainty
%          No_of_iterations - number of iterations made by the line search
%          algorithm

alpha = 0.618; % Golden section number
No_of_iterations = 0;
max_no_of_iterations = 200;
tol = 10^-6; % The tolerance which determines when the interval of uncertainty is small enough
a = 0; % The initial left end point of the interval
b = 10^6; % The initial right end point of the interval (should be set to a high number)
lambda = a + (1-alpha)*(b-a);
flambda = func(x+lambda*d);
fa = func(x+a*d);
fb = func(x+b*d);

% Handles the case where a function varies very slowly.
while fa-fb > tol || norm(fa-flambda) < tol
    a = lambda;
    btemp = b;
    b = b + (1+alpha)*(b-lambda);
    lambda = btemp;
    fa = func(x+a*d);
    fb = func(x+b*d);
    flambda = func(x+lambda*d);
end

mu = a + alpha*(b-a);
fmu = func(x+mu*d);
while b-a > tol && No_of_iterations <= max_no_of_iterations
    if flambda > fmu
        a = lambda;
        lambda = mu;
        mu = a + alpha*(b-a);
        flambda = fmu;
         if b-a > tol
            fmu = func(x+mu*d);
         end
      else
          b = mu;
          mu = lambda;
          lambda = a + (1-alpha)*(b-a);
          fmu = flambda;
           if b-a > tol
            flambda = func(x+lambda*d);
           end
    end
    No_of_iterations = No_of_iterations + 1;
end


% if isnan(func(x+lambda*d)) || func(x+lambda*d)>func(x)
%     error('Bad job of the line search!')
% end
