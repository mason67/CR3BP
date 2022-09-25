function [A, alpha, lambda] = getLinearizedIPVals(xEq, mu, xI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GETLINEARIZEDAVALS Returns Ai, alpha<1,3> and lambda<1,3> values for
%in-plane linearized dynamics around an equilibrium point.
%   Detailed explanation goes here
%   Inputs:
%       xEq -   [1x3] X,Y,Z data of equilibrium point to study
%       mu -    [fl]  Mass ratio of the system 
%       xI -    [4x1] X,Y,Xd,Yd initial conditions of particle about eq. pt
%   Outputs:
%       A  -     [4x1] Ai components of system
%       alpha -  [1x4] Calculated alpha<1-4> values
%       lambda - [1x4] Calculated lambda<1-4> values
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % First, evaluate the variational eqns at the specific equilibrium pt:
    U = evalVariationalEqns(xEq, mu);

    % Get modes:
    [~, LIPLeq] = getEquilModes(xEq, mu);

    % find the real and imaginary components...
    foundL1 = false;
    foundL3 = false;
    for i =1:length(LIPLeq)
        if isreal(LIPLeq(i)) && ~foundL1
            L1 = LIPLeq(i);
            L2 = -L1;
            foundL1 = true;
        elseif imag(LIPLeq(i)) && ~foundL3
            L3 = LIPLeq(i);
            L4 = conj(L3);
            foundL3 = true;
        end
    end

    % Now calculate alphas:
    alpha1 = (L1^2 - U(1,1))/(2*L1);
    alpha2 = (L2^2 - U(1,1))/(2*L2);
    alpha3 = (L3^2 - U(1,1))/(2*L3);
    alpha4 = (L4^2 - U(1,1))/(2*L4);

    % Now, use initial conditions to get the A values:
    th  = 1/(2*(L1*alpha1 - L3*alpha3));
    phi = 1/(2*(L1*alpha3 - L3*alpha1));
    A = [-th*alpha3*L3  phi*alpha3 -phi*L3  th;
         -th*alpha3*L3 -phi*alpha3  phi*L3  th;
          th*alpha1*L1 -phi*alpha1  phi*L1 -th;
          th*alpha1*L1  phi*alpha1 -phi*L1 -th]*xI;
    
    % return alpha and lambdas too
    alpha = [alpha1 alpha2 alpha3 alpha4];
    lambda = [L1 L2 L3 L4];
end

