x1 = linspace(0,10,100);
x2 = linspace(0,10,100);
[X,Y] = meshgrid(x1,x2);
F = zeros(100,100);
F2 = zeros(100,100);

% phi1
[t1,y1] = data1;
r_i = @(x,t,y) (phi1(x,t) - y);
f = @(x) 0; % Initialize the objective function as the zero function
% Constructs the objective function with fitting function phi1
for i = 1:length(t1)
    f = @(x) f(x) + (r_i(x,t1(i),y1(i)))^2;
end
% Computes the value of the objective function at the grid points
for i = 1:length(x1)
    for j = 1:length(x2)
        F(i,j) = f([x1(i) x2(j)]);
    end
end
figure
surf(X,Y,F)
