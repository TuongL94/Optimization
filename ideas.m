fa = func(x+a*d);
fb = func(x+b*d);
n = 1;
while fa > fb 
    a = 1.5^n*a+1;
    fa = func(x+a*d);
    n = n+1;
end

% while fa < fb
%     b = b/(1.5^n);
%     fb = func(x+b*d);
%     n = n+1;
% end

while fa < flambda && flambda < func(x+b*d)
    b = b/(2^(alpha^n));
    lambda = a + (1-alpha)*(b-a);
    flambda = func(x+lambda*d);
    n = n + 1;
end

% if func(x) < Inf
%     while flambda == Inf
%         b = b/2;
%         lambda = a + (1-alpha)*(b-a);
%         flambda = func(x+lambda*d);
%     end
% end
% 
% if func(a+(b-a)*tol*d) < func(b-(b-a)*tol*d)
%     lambda = (b-a)*tol;
% else
%     lambda = b-(b-a)*tol;
% end