function varargout = calcZVC(C, mu, step, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCZVC Calculates pseudopotential of a CR3BP system and plots the ZVCs at
%a given Jacobi Constant, C.
%   Inputs:
%       C -     [fl]  Desired Jacobi Constant value for the CR3BP system
%       mu -    [fl]  Mass ratio of the system 
%       step -  [fl]  Step size to use for X,Y meshgrid
%       varargin - matlab variable argument in, supported params:
%           -# 'plot' [bool] flag to plot ZVCs, default: true
%           -# 'sys'  [str]  System identifier, i.e. 'Earth-Moon'
%   Outputs:
%       X -    [mxm] meshgrid matrix of X coordinates
%       Y -    [mxm] meshgrid matrix of Y coordinates
%       CMAT - [mxm] meshgrid matrix of jacobi constants
%       h -    [--]  figure handle for the plot
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse variable inputs
    P = inputParser;
    addParameter(P, 'plot', true, @(x) assert(isbool(x),'plot input must be boolean.'))
    addParameter(P, 'sys', 'System', @(x) assert(isstring(x) || ischar(x), 'sys input must be string.'))
    parse(P,varargin{:});

    % Create X,Y meshgrid for pseudopotential calculation
    x = -2.0:step:2.0;
    [X,Y] = meshgrid(x);

    % Calculate the pseudo potential for all values X and Y
    U = 0.5*(X.^2 + Y.^2) + (1 - mu)./sqrt((X + mu).^2 + Y.^2) + mu./sqrt((X - 1 + mu).^2 + Y.^2);
    CMAT = 2*U;

    % Plot the ZVC at the given Jacobi Constant
    if P.Results.plot
        h = figure;
        contourf(X, Y, CMAT, [C C])
        
        % Add primaries
        hold on
        plot(-mu,  0, 'kx')
        plot(1-mu, 0, 'rx')
        hold off
        axis equal
    
        % Modify color map 
        map = [0.6 0.6 0.6; 1 1 1];
        colormap(map)
        set(gca, 'color', [0.6 0.6 0.6])
        legend('','P1','P2','color','none')
        xlabel('X [--]')
        ylabel('Y [--]')
        title(sprintf('%s ZVC at C = %0.3f', P.Results.sys, C))
    end
    varargout = {X,Y,CMAT,h};

end

