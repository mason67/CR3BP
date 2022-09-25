function [lambda_oop, lambda_ip] = getEquilModes(R, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GETEQUILMODES Returns the in plane and out of plane modes for equilibrium
%points in the CR3BP
%   Inputs:
%       R -     [fl]  [mx3] Position data in X,Y,Z format
%       mu -    [fl]  Mass ratio of the system 
%   Outputs:
%       lambda_oop - [mx2] Out-of-plane modes
%       lambda_ip  - [mx4] In-plane modes
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Evaluate variational equations for given position data
    U = evalVariationalEqns(R, mu);
    
    % Define local vars for ease of use
    Uxx = U(:,1); Uxy = U(:,2); Uxz = U(:,3);
    Uyx = U(:,4); Uyy = U(:,5); Uyz = U(:,6);
    Uzx = U(:,7); Uzy = U(:,8); Uzz = U(:,9);

    % Calculate the modes
    % out of plane
    lambda_oop = 1i.*sqrt(abs(Uzz))*[1 -1]; 
    
    % In plane
    Lambda_ip  = (-4 + Uxx + Uyy)./2 + (sqrt((4 - Uxx - Uyy).^2 - 4.*(Uxx.*Uyy - Uxy.^2))./2)*[1 -1];
    lambda_ip_pos = sqrt(Lambda_ip);
    lambda_ip_neg = -lambda_ip_pos;
    lambda_ip = [lambda_ip_pos lambda_ip_neg];
    
end

