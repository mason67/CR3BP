function dydt = CR3BP_STM(t,yi,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CR3BP Equations of Motion for the CR3BP defined for integration schemes in
%the rotating frame.
%   Inputs:
%       t  - [--]  Npn dimensional time, t
%       y  - [42x1] Non-dimensional state vector in the rotating frame,
%            expecting coordinates defined as [X,Y,Z,dX,dY,dZ, STMij]. NOTE:
%            STM Components are expected to be in ROW MAJOR format, e.g.
%            [..., STM(1,1), STM(1,2), STM(1,3), ..., STM(6,6)]
%       mu - [--] Mass ratio of the CR3BP system
%   Outputs:
%       dydt - [42x1] Equations of motion definining the CR3BP system with
%              STM integration
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % local copies to help with understanding
    x  = yi(1); y  = yi(2); z  = yi(3);
    xd = yi(4); yd = yi(5); zd = yi(6);
    
    % define rotating frame vectors r1, r2:
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    % Evaluate Variational Equations at input position
    U = evalVariationalEqns([x, y, z], mu);
    Uxx = U(1); Uxy = U(2); Uyx = U(4); Uyy = U(5); Uzz = U(9);

    % Construct 6x6 variational matrix, A(t):
    A = [0   0   0    1 0 0;
            0   0   0    0 1 0;
            0   0   0    0 0 1;
            Uxx Uxy 0    0 2 0;
            Uyx Uyy 0   -2 0 0;
            0   0   Uzz  0 0 0];
    % Flatten into a 1x36 vector with ROW MAJOR ORDER:
%     A = reshape(Amat',1,[]); % default behavior is column major, so transpose to get "row major"
    
    % Construct integration output
    dydt = zeros(42,1);

    % State Vector Derivative
    dydt(1) = xd;
    dydt(2) = yd;
    dydt(3) = zd;
    dydt(4) =  2*yd + x - (1 - mu)*(x + mu)/(r1^3) - mu*(x - 1 + mu)/(r2^3);
    dydt(5) = -2*xd + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    dydt(6) = -(1 - mu)*z/(r1^3) - mu*z/(r2^3);

    % STM Derivative: Phid(tn, tn-1) = A(t)*Phi(tn, tn-1)
    PHI = [yi(7:12)';
           yi(13:18)';
           yi(19:24)';
           yi(25:30)';
           yi(31:36)';
           yi(37:42)'];
    PHId = A*PHI;
    dydt(7:42) = reshape(PHId',1,[]);

end

