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
% if dim == 2
%     F = zeros(100,100);
%     x1 = y2(1);
%     x2 = 2;
%     x1_grid = linspace(0,x1,100);
%     x2_grid = linspace(0,x2,100);
%     % Computes the value of the objective function at the grid points
%     for i = 1:length(x1_grid)
%         for j = 1:length(x2_grid)
%             F(i,j) = f_1([x1_grid(i) x2_grid(j)]);
%         end
%     end
%     [fmin,ind] = min(F(:));
%     [m,n] = ind2sub(size(F),ind);
%     xstart = [x1_grid(m) ; x2_grid(n)];
%     min = gaussnewton(@phi1,t2,y2,xstart,tol,0,1,0);
% elseif dim == 4
%     x1 = y2(1);
%     x2 = 10;
%     x3 = y2(1);
%     x4 = 10;
%     x1_grid = linspace(0,x1,20);
%     x2_grid = linspace(0,x2,20);
%     x3_grid = linspace(0,x3,20);
%     x4_grid = linspace(0,x4,20);
%     fmin = realmax;
%     index_vec = zeros(4,1);
%     for i = 1:length(x1_grid)
%         for j = 1:length(x2_grid)
%             for k = 1:length(x3_grid)
%                 for l = 1:length(x4_grid)
%                     fval = f_2([x1_grid(i);x2_grid(j);x3_grid(k);x4_grid(l)]);
%                     if fval < fmin
%                         fmin = fval;
%                         index_vec(1) = i;
%                         index_vec(2) = j;
%                         index_vec(3) = k;
%                         index_vec(4) = l;
%                     end
%                 end
%             end
%         end
%     end
%     xstart = [x1_grid(index1);x2_grid(index2);x3_grid(index3);x4_grid(index4)];
%     min = gaussnewton(@phi2,t2,y2,xstart,tol,0,1,0);
% end
%%

xstart = [3;8;5;1];
min = gaussnewton(@phi2,t1,y1,xstart,tol,1,1,0);

%%
% Test case: phi1 and data1
n_points = 20;

x_low = 4;
x_high = 16;

y_low = -2;
y_high = 6;

x_interval = x_low + (x_high-x_low).*rand(n_points,1);

y_interval = y_low + (y_high-y_low).*rand(n_points,1);

converge = zeros(1,n_points);

for i = 1:n_points
    x_start = [x_interval(i); y_interval(i)];
    xmin =  gaussnewton(@phi1,t1,y1,x_start,tol,0,0,0);
    diff = norm(gaussnewton(@phi1,t1,y1,[10; 2],tol,1,0,0)-xmin);
    if diff < tol*10^3
        converge(i) = 1;
    end
end
converge
%%
% Test case: phi1 and data2
n_points = 20;

x_low = 8;
x_high = 17;
y_low = -3;
y_high = 6;

x_interval = x_low + (x_high-x_low).*rand(n_points,1);
y_interval = y_low + (y_high-y_low).*rand(n_points,1);

converge = zeros(1,n_points);
for i = 1:n_points
    x_start = [x_interval(i); y_interval(i)];
    xmin =  gaussnewton(@phi1,t2,y2,x_start,tol,0,0,0);
    diff = norm(gaussnewton(@phi1,t2,y2,[10; 2],tol,1,0,0)-xmin);
    if diff < tol*10^3
        converge(i) = 1;
    end
end
converge
%%
% Test case: phi2 and data1
n_points = 10;
xmin_real = [6.3446; 10.5867; 6.0958; 1.4003];
x1 = y1(1);
x2 = 10;
x3 = y1(1);
x4 = 10;

x13_interval = x1*rand(n_points,1);
x24_interval = x2.*rand(n_points,1);

converge = zeros(1,n_points);
for i = 1:n_points
    x_start = [x13_interval(i); x24_interval(i); x13_interval(i); x24_interval(i)];
    xmin =  gaussnewton(@phi2,t1,y1,x_start,tol,1,0,0);
    diff = norm(xmin_real - xmin);
    if diff < tol*10^3
        converge(i) = 1;
    end
end
converge
     
        
