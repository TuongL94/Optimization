function [xmin] = gaussnewtonwithprint(phi,t,y,start,tol,use_linesearch,printout,plotout)
%GAUSSNEWTON Summary of this function goes here
%   Detailed explanation goes here

if length(start) ~= 2 && length(start) ~= 4
    error('Please choose a starting point of 2 or 4 dimensions.')
end
    
% Parameters
max_no_of_iterations = 101;
no_of_iterations = 0;

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
    d = (J(xcurrent)'*J(xcurrent))\(-J(xcurrent)'*r(xcurrent));
    if use_linesearch == 1
        lambda = linesearch(f,xcurrent,d);
        xnew = xcurrent+lambda*d;
    elseif use_linesearch == 0
        xnew = xcurrent+d;
    else
        error('Please enter the number 0 or 1 in the sixth argument of the function');
    end
    no_of_iterations = no_of_iterations + 1;
    fcurrent = f(xcurrent);
    rel_diff = abs(f(xnew)-fcurrent)/abs(fcurrent);
    grad_f = 2*J(xnew)'*r(xnew);
    norm_grad_f = norm(grad_f);
    xcurrent = xnew;
    % Print out
    if printout == 1
        if no_of_iterations >= 100
            fprintf('%s %7s %12s %13s %13s\n','iter', 'x', 'f(x)','norm(grad)','rel.diff');
            fprintf('%d %10.2f %10.2f %10.2f %14.2f\n',no_of_iterations,xcurrent(1),fcurrent,norm_grad_f,rel_diff);
            for i = 2:length(start)
                fprintf('%s %12.2f %10.2s %10s\n',' ',xcurrent(i),' ',' ');
            end
        elseif no_of_iterations >= 10
            fprintf('%s %6s %12s %13s %13s\n','iter', 'x', 'f(x)','norm(grad)','rel.diff');
            fprintf('%d %10.2f %10.2f %10.2f %14.2f\n',no_of_iterations,xcurrent(1),fcurrent,norm_grad_f,rel_diff);
            for i = 2:length(start)
                fprintf('%s %11.2f %10.2s %10s\n',' ',xcurrent(i),' ',' ');
            end
        else   
            fprintf('%s %5s %12s %13s %13s\n','iter', 'x', 'f(x)','norm(grad)', 'rel.diff');
            fprintf('%d %10.2f %10.2f %10.2f %14.2f\n', no_of_iterations,xcurrent(1),fcurrent,norm_grad_f,rel_diff);
            for i = 2:length(start)
                fprintf('%s %10.2f %10.2s %10s\n',' ', xcurrent(i),' ',' ');
            end
        end
    end 
    % Termination criterion
    if (rel_diff < tol) && (norm_grad_f < tol)
        break;
    end
    
end

xmin = xcurrent;

if plotout == 1
    grid = linspace(floor(min(t)),ceil(max(t)),100);
    plot(grid,phi(xmin,grid))
    hold on
    plot(t,y,'ro')
end



