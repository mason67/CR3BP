function [Lx,Fx] = getLagrangePts(mu, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GETLAGRANGEPTS Calculates the Lagrange/Equlibrium/Libration Points for 
%   a specific system mu with a tolerance tol
%   Inputs:
%       mu -    [fl]  Mass ratio of the system 
%       tol -   [fl]  F(x) tolerance for L1-L3
%   Outputs:
%       Lx -   [5x3] Matrix of L<1-5> x,y,z points
%       Fx -   [5x3] Matrix of L<1-5> F(x) outputs 
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Lx = zeros(5,3); 
    Fx = Lx;
    % Calculate L1-L3 x coordinates using fsolve
    xIG = [1 - mu - 0.1, 1 - mu + 0.05, -mu - 0.05];

    dudx = @(x) x - (1 - mu)*(x + mu)/(abs(x + mu)^3) - mu*(x - 1 + mu)/(abs(x - 1 + mu)^3);
    options = optimoptions('fsolve', 'FunctionTolerance', tol, 'OptimalityTolerance', tol, 'Display', 'off');
    k = 1;
    for xI = xIG
        [x, fval] = fsolve(dudx,xI,options);
        if fval < tol
            Lx(k,1) = x;
            Fx(k,1) = fval;
            k = k + 1;
        else
            error('fsolve could not solve to the specified tolerance for L%d',k)
        end
    end
    
    % Next, find the triangular equilibrium points L4, L5:
    Lx(4,:) = [1/2 - mu,  sqrt(3)/2, 0];
    Lx(5,:) = [1/2 - mu, -sqrt(3)/2, 0];

end

