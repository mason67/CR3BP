function U = evalVariationalEqns(R, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EVALVARIATIONALEQNS Calculates the Uxx, Uxy, Uxz... Variational Equations
%in the CR3BP
%   Inputs:
%       R -     [fl]  [mx3] Position data in X,Y,Z format
%       mu -    [fl]  Mass ratio of the system 
%   Outputs:
%       U - [Mx9] Matrix of variational equations, ROW MAJOR ORDER
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % separate out the state vars
    x = R(:,1);
    y = R(:,2);
    z = R(:,3);
    
    % Calculate r1, r2:
    r1 = sqrt((x + mu).^2 + y.^2 + z.^2);
    r2 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);

    % Calculate variational equations
    Uxx = 1 - (1 - mu)./(r1.^3) - mu./(r2.^3)+ (3.*(1 - mu).*(x + mu).^2)./(r1.^5) + (3.*mu.*(x - 1 + mu).^2)./(r2.^5);
    Uxy = (3.*(1 - mu).*(x + mu).*y)./(r1.^5) + (3.*mu.*(x - 1 + mu).*y)./(r2.^5);
    Uxz = (3.*(1 - mu).*(x + mu).*z)./(r1.^5) + (3.*mu.*(x - 1 + mu).*z)./(r2.^5);
    Uyx = Uxy;
    Uyy = 1 - (1 - mu)./(r1.^3) - mu./(r2.^3) + (3.*(1 - mu).*(y.^2))./(r1.^5) + (3.*mu.*(y.^2))./(r2.^5);
    Uyz = (3.*(1 - mu).*y.*z)./(r1.^5) + (3.*mu.*y.*z)./(r2.^5);
    Uzx = Uxz;
    Uzy = Uyz;
    Uzz = -(1 - mu)./(r1.^3) - mu./(r2.^3) + (3.*(1 - mu).*(z.^2))./(r1.^5) + (3.*mu.*(z.^2))./(r2.^5);

    U = [Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz];
end

