% project work
tol = 10^-3;
x = [1 2 3 4]';
[t_data, y_data] = data2;

m = length(t_data);

r_i = @(x,t,y) (phi2(x,t) - y);
%r_i = @(x,t,y) (phi1(x,t) - y);
grad_ri = @(x,t) [exp(-x(2)*t) -x(1)*t.*exp(-x(2)*t) exp(-x(4)*t) -x(3)*t.*exp(-x(4)*t)];
%grad_ri = @(x,t) [exp(-x(2)*t) -x(1)*t.*exp(-x(2)*t)]; 
grad_r = @(x) grad_ri(x, t_data);
r_vec = @(x) r_i(x,t_data, y_data);
f = @(x) sum(r_vec(x).^2);

while true
    J = grad_r(x);
    r = r_vec(x);
    d = (J'*J)\(-J'*r);
    lambda = armijo2(f,x,d);
    xnew = x+lambda*d;
    rel_diff = abs(f(xnew)-f(x))/abs(f(x));
    grad_f = 2*J'*r;
    if (rel_diff < tol) && (norm(grad_f) < tol)
        break;
    end
    x = xnew;
end
f(x)
grid = linspace(0,2,100);
figure(1)
plot(grid,phi2(x,grid))
figure(2)
plot(t_data,y_data,'o')