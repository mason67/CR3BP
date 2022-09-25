function dydt = linIPEOM(t,y,U)
%LINIPEOM Linearized equations of motion 
%   Inputs:
%       t - [--]  Npn dimensional time, t
%       y - [6x1] Non-dimensional state vector in the rotating frame,
%           expecting coordinates defined as [X,Y,Z,dX,dY,dZ]
%       U - [1x9] Vector of Uxx,Uxy,...,Uzz values
%   Outputs:
%       dydt - [6x1] Equations of motion definining the linearized system
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Uxx = U(1);
    Uxy = U(2);
    Uyx = U(4);
    Uyy = U(5);
    [m,~] = size(y);
    if m ~= 4
        y = y';
    end

    % Use A2d
    dydt = [0   0    1 0;
            0   0    0 1; 
            Uxx Uxy  0 2;
            Uyx Uyy -2 0]*y;

end

