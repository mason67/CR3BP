function dydt = CR3BP(t,yi,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CR3BP Equations of Motion for the CR3BP defined for integration schemes in
%the rotating frame.
%   Inputs:
%       t - [--]  Npn dimensional time, t
%       y - [6x1] Non-dimensional state vector in the rotating frame,
%           expecting coordinates defined as [X,Y,Z,dX,dY,dZ]
%       mu - [--] Mass ratio of the CR3BP system
%   Outputs:
%       dydt - [6x1] Equations of motion definining the CR3BP system
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % local copies to help with understanding
    x  = yi(1); y  = yi(2); z  = yi(3);
    xd = yi(4); yd = yi(5); zd = yi(6);
    
    % define rotating frame vectors r1, r2:
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    dydt(1) = xd;
    dydt(2) = yd;
    dydt(3) = zd;
    dydt(4) =  2*yd + x - (1 - mu)*(x + mu)/(r1^3) - mu*(x - 1 + mu)/(r2^3);
    dydt(5) = -2*xd + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    dydt(6) = -(1 - mu)*z/(r1^3) - mu*z/(r2^3);
    
    dydt = dydt';

end

