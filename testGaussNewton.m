[t1 y1] = data1;
[t2 y2] = data2;
tol = 10^-6;

%%
dim = 4;
r_1 = @(x) (phi1(x,t2) - y2);
f_1 = @(x) sum(r_1(x).^2);
r_2 = @(x) (phi2(x,t2) - y2);
f_2 = @(x) sum(r_2(x).^2);
% Upper bounds for initial point for Gauss Newton method 
% wihtout linesearch. Lower bound is 0 for all variables.
if dim == 2
    F = zeros(100,100);
    x1 = y2(1);
    x2 = 2;
    x1_grid = linspace(0,x1,100);
    x2_grid = linspace(0,x2,100);
    % Computes the value of the objective function at the grid points
    for i = 1:length(x1_grid)
        for j = 1:length(x2_grid)
            F(i,j) = f_1([x1_grid(i) x2_grid(j)]);
        end
    end
    [fmin,ind] = min(F(:));
    [m,n] = ind2sub(size(F),ind);
    xstart = [x1_grid(m) ; x2_grid(n)];
    min = gaussnewton(@phi1,t2,y2,xstart,tol,0,1,0);
elseif dim == 4
    x1 = y2(1);
    x2 = 10;
    x3 = y2(1);
    x4 = 10;
    x1_grid = linspace(0,x1,20);
    x2_grid = linspace(0,x2,20);
    x3_grid = linspace(0,x3,20);
    x4_grid = linspace(0,x4,20);
    fmin = realmax;
    index_vec = zeros(4,1);
    for i = 1:length(x1_grid)
        for j = 1:length(x2_grid)
            for k = 1:length(x3_grid)
                for l = 1:length(x4_grid)
                    fval = f_2([x1_grid(i);x2_grid(j);x3_grid(k);x4_grid(l)]);
                    if fval < fmin
                        fmin = fval;
                        index_vec(1) = i;
                        index_vec(2) = j;
                        index_vec(3) = k;
                        index_vec(4) = l;
                    end
                end
            end
        end
    end
    xstart = [x1_grid(index1);x2_grid(index2);x3_grid(index3);x4_grid(index4)];
    min = gaussnewton(@phi2,t2,y2,xstart,tol,0,1,0);
end

%%

xstart = [2 5 4 6]';
min = gaussnewton(@phi2,t2,y2,xstart,tol,1,1,0);

%%
% Test case: phi1 and data1
n_points = 20; % number of random points
xmin_real = [10.7959;2.4790]; % real minimum
epsilon = 1; % length of the  2-dim cube centered around the real mimimum

x1_interval = (xmin_real(1) - 1) + (2*epsilon).*rand(n_points,1);
x2_interval = (xmin_real(2) - 1) + (2*epsilon).*rand(n_points,1);

converge = zeros(1,n_points);
for i = 1:n_points
    x_start = [x1_interval(i); x2_interval(i)];
    xmin =  gaussnewton(@phi1,t1,y1,x_start,tol,0,0,0);
    diff = norm(xmin_real-xmin);
    if diff < tol*10^3
        converge(i) = 1;
    end
end
converge

%%
% Test case: phi1 and data2
n_points = 20; % number of random points
xmin_real = [12.9826;1.7864]; % real minimum
epsilon = 1; % length of the  2-dim cube centered around the real mimimum

x1_interval = (xmin_real(1) - 1) + (2*epsilon).*rand(n_points,1);
x2_interval = (xmin_real(2) - 1) + (2*epsilon).*rand(n_points,1);

converge = zeros(1,n_points);
for i = 1:n_points
    x_start = [x1_interval(i); x2_interval(i)];
    xmin =  gaussnewton(@phi1,t2,y2,x_start,tol,0,0,0);
    diff = norm(xmin_real-xmin);
    if diff < tol*10^3
        converge(i) = 1;
    end
end
converge

%%
% Test case: phi2 and data1
n_points = 10; % number of random points
xmin_real = [6.3446; 10.5867; 6.0958; 1.4003]; % real minimum
epsilon = 1; % length of the 4-dim cube centered around the real minimum

x1_interval = (xmin_real(1)-epsilon) + (2*epsilon).*rand(n_points,1);
x2_interval = (xmin_real(2)-epsilon) + (2*epsilon).*rand(n_points,1);
x3_interval = (xmin_real(3)-epsilon) + (2*epsilon).*rand(n_points,1);
x4_interval = (xmin_real(4)-epsilon) + (2*epsilon).*rand(n_points,1);

converge = zeros(1,n_points);
for i = 1:n_points
    x_start = [x1_interval(i); x2_interval(i); x3_interval(i); x4_interval(i)];
    xmin =  gaussnewton(@phi2,t1,y1,x_start,tol,1,0,0);
    diff = norm(xmin_real - xmin);
    if diff < tol*10^5
        converge(i) = 1;
    end
end
converge

%%
% Test case: phi2 and data2
n_points = 20;
xmin_real = [9.7389; 2.9208; 4.1742; 0.8748];
epsilon = 1;

x1_interval = (xmin_real(1)-epsilon) + 2*epsilon.*rand(n_points,1);
x2_interval = (xmin_real(2)-epsilon) + 2*epsilon.*rand(n_points,1);
x3_interval = (xmin_real(3)-epsilon) + 2*epsilon.*rand(n_points,1);
x4_interval = (xmin_real(4)-epsilon) + 2*epsilon.*rand(n_points,1);

converge = zeros(1,n_points);

for i = 1:n_points
    x_start = [x1_interval(i); x2_interval(i); x3_interval(i); x4_interval(i)];
    xmin =  gaussnewton(@phi2,t2,y2,x_start,tol,1,0,0);
    diff = norm(xmin_real - xmin);
    if diff < tol*10^4
        converge(i) = 1;
    end
end
sum(converge)/n_points
     
        
