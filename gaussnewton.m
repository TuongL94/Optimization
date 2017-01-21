function [xmin] = gaussnewton(phi,t,y,start,tol,use_linesearch,printout,plotout)
% Input:   phi - fitting function
%            t - independent data point values
%            y - dependent data point values
%        start - initial point, given as a column vector
%          tol - tolerance used for stopping criterion
%        use_linesearch - determines if linesearch is to be used. 0 means
%        that linesearch will not be used, 1 implies that it will be used.
%        printout - Setting this to 1 generates printout of the
%        optimization process
%        plotout - Setting this to 1 generates plots of the data points and
%        the fitted function
%
% Output:   xmin - the optimal parameters for the fitting function

if length(start) ~= 2 && length(start) ~= 4
    error('Please choose a starting point of 2 or 4 dimensions.')
end
    
% Parameters
max_no_of_iterations = 100;
no_of_iterations = 0;
norm_tol = 0.5; % tolerance for the norm of the gradient at local minimum 

% Define the residual functions and their gradients. If the start point is
% of 2 variables the residual functions and the gradients are defined for
% the fitting function phi1. If start is of 4 variables the procedure is
% done for phi2.
if length(start) == 2
    r = @(x) (phi(x,t) - y);
    J = @(x) [exp(-x(2)*t) -x(1)*t.*exp(-x(2)*t)]; 
    f = @(x) sum(r(x).^2);
elseif length(start) == 4
    r = @(x) (phi(x,t) - y);
    J = @(x) [exp(-x(2)*t) -x(1)*t.*exp(-x(2)*t) exp(-x(4)*t) -x(3)*t.*exp(-x(4)*t)];
    f = @(x) sum(r(x).^2);    
end

xcurrent = start;
for i = 1:max_no_of_iterations
    epsilon = 0.01;
    H = J(xcurrent)'*J(xcurrent);
    H_mod = H;
    [~,p] = chol(H_mod);
    while (p > 0)
        H_mod = H + epsilon*eye(size(H,1));
        [~,p] = chol(H_mod);
        epsilon = 4*epsilon;
    end
    H = H_mod;
    epsilon = 0.01;
    while cond(H_mod) > 10000000
        H_mod = H + epsilon*eye(size(H,1));
        epsilon = 4*epsilon;
    end
    H = H_mod;
    d = H\(-J(xcurrent)'*r(xcurrent));
    if use_linesearch == 1
        lambda = linesearch(f,xcurrent,d);
        xnew = xcurrent+lambda*d;
    elseif use_linesearch == 0
        xnew = xcurrent+d;
    else
        error('Please enter the number 0 or 1 in the sixth argument of the function');
    end
    fcurrent = f(xcurrent);
    rel_diff = abs(f(xnew)-f(xcurrent))/abs(f(xcurrent));
   
    grad_f = 2*J(xnew)'*r(xnew);
    norm_grad_f = norm(grad_f);
    xcurrent = xnew;
    
    if no_of_iterations == 0 && printout == 1
            fprintf('%s %7s %15s %15s\n','iter', 'x', 'f(x)','norm(grad)');
            fprintf('%d %13.4f %13.4f %12.4f\n',no_of_iterations,start(1),fcurrent,norm_grad_f);
            for i = 2:length(start)
                fprintf('%s %13.4f %10.2s\n',' ',start(i),' ');
            end
    end
    
    no_of_iterations = no_of_iterations + 1;
    
    % Print out
    if printout == 1
       if no_of_iterations >= 10
            fprintf('%s %7s %15s %15s %15s\n','iter', 'x', 'f(x)','norm(grad)','rel.diff(f)');
            fprintf('%d %12.4f %13.4f %12.4f %16.7f\n',no_of_iterations,xcurrent(1),...
            fcurrent,norm_grad_f,rel_diff);
            for i = 2:length(start)
               fprintf('%s %13.4f %13.2s %10s\n',' ',xcurrent(i),' ',' ');
            end
        else   
            fprintf('%s %7s %15s %15s %15s\n','iter', 'x', 'f(x)','norm(grad)','rel.diff(f)');
            fprintf('%d %13.4f %13.4f %12.4f %16.7f\n',no_of_iterations,xcurrent(1),...
            fcurrent,norm_grad_f,rel_diff);
            for i = 2:length(start)
                fprintf('%s %13.4f %10.2s %10s\n',' ',xcurrent(i),' ',' ');
            end
        end
    end 
    
    % Termination criterion
    if (rel_diff < tol) && norm_grad_f < norm_tol
        break;
    end
    
    if no_of_iterations == max_no_of_iterations && ((rel_diff > tol) || (norm_grad_f > norm_tol))
        disp('The algorithm failed to converge to a local minimum within the maximum number of iterations')
    end
end
xmin = xcurrent;

if plotout == 1
    grid = linspace(floor(min(t)),ceil(max(t)),100);
    plot(grid,phi(xmin,grid))
    hold on
    plot(t,y,'ro')
end




