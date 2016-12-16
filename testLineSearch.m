a = 10;
d = 1;
xstart = -0.5;
func = @(x) (1-10^a*x).^2;
% xgrid = linspace(-200,200,1000);
% plot(xgrid,func(xgrid))
[lambda,No_of_iterations] = linesearch(func,xstart,d);
func(xstart+lambda*d)

%lambda_1=linesearch(@test_func,[0;0],[1;0]);
%[lambda_2,No_of_iterations]=linesearch(@test_func,[0;0],[0;1]);

