function C = calcJacobiConst(R, V, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCJACOBICONST Calculates Jacobi constant at X,Y,Z values
%   Inputs:
%       R -     [fl]  [mx3] Position data in X,Y,Z format
%       V -     [fl]  [mx3] Velocity data in X,Y,Z format
%       mu -    [fl]  Mass ratio of the system 
%   Outputs:
%       C - [mx1] vector of Jacobi constants
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % separate out the state vars
    x = R(:,1);
    y = R(:,2);
    z = R(:,3);
    xd = V(:,1);
    yd = V(:,2);
    zd = V(:,3);

    % Calculate the Jacobi Const, C:
    r1 = sqrt((x + mu).^2 + y.^2 + z.^2);
    r2 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);
    C  = x.^2 + y.^2 + 2*(1 - mu)./r1 + 2*mu./r2 - xd.^2 - yd.^2 - zd.^2;
end

