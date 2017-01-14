[t1 y1] = data1;
[t2 y2] = data2;
tol = 10^-5;
xstart = [5; 5; 5; 5];

min = gaussnewton(@phi2,t1,y1,xstart,tol,0,1,1);

% n_points = 20;
% 
% x_low = 4;
% x_high = 16;
% 
% y_low = -2;
% y_high = 6;
% 
% x_interval = x_low + (x_high-x_low).*rand(n_points,1);
% 
% y_interval = y_low + (y_high-y_low).*rand(n_points,1);
% 
% converge = zeros(1,n_points);
% 
% for i = 1:n_points
%     x_start = [x_interval(i); y_interval(i)];
%     xmin =  gaussnewton(@phi1,t1,y1,x_start,tol,0,0,0);
%     diff = norm(gaussnewton(@phi1,t1,y1,[10; 2],tol,1,0,0)-xmin);
%     if diff < tol*10^3
%         converge(i) = 1;
%     end
% end
% converge
%         


% n_points = 20;
% 
% x_low = 8;
% x_high = 17;
% 
% y_low = -3;
% y_high = 6;
% 
% x_interval = x_low + (x_high-x_low).*rand(n_points,1);
% 
% y_interval = y_low + (y_high-y_low).*rand(n_points,1);
% 
% converge = zeros(1,n_points);
% 
% for i = 1:n_points
%     x_start = [x_interval(i); y_interval(i)];
%     xmin =  gaussnewton(@phi1,t2,y2,x_start,tol,0,0,0);
%     diff = norm(gaussnewton(@phi1,t2,y2,[10; 2],tol,1,0,0)-xmin);
%     if diff < tol*10^3
%         converge(i) = 1;
%     end
% end
% converge
     
        