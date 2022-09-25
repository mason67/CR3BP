function varargout = plotOrbit(y, ttitle, sys, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTORBIT Plots a trajectory in the CR3BP rotating frame with options 
% to plot primaries.
%   Inputs:
%       y -      [mx3] matrix of X,Y,Z coordinates in the rotating frame.
%       ttitle - [str] Title for the plot
%       sys -    [struct] structure that contains at the very least the 
%                system mu value.
%       varargin - MATLAB variable argument in, supports:
%           -# pPrim <both|P1|P2> tells it which primary(ies) to plot
%           -# p1Size <0.0-1.0> scale factor of the larger primary
%           -# p2Size <0.0-1.0> scale factor of the smaller primary
%   Outputs:
%       h -    [--]  figure handle for the plot
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %parse inputs
    P = inputParser;
    addParameter(P, 'pPrim', 'both', @(x) assert(strcmp(x,'both') || strcmp(x,'P1') || strcmp(x,'P2') , 'both input must be: P1, P2 or both.'))
    addParameter(P, 'p1Size', sys.P1Size, @(x) assert(x<1.0 && x>0.0, 'p1Size must be between 0.0 and 1.0.')) 
    addParameter(P, 'p2Size', sys.P2Size, @(x) assert(x<1.0 && x>0.0, 'p1Size must be between 0.0 and 1.0.')) 
    parse(P,varargin{:});

    h = figure;
    plot3(y(:,1),y(:,2),y(:,3))
    xlabel('X [nondim]')
    ylabel('Y [nondim]')
    zlabel('Z [nondim]')
    title(ttitle)
    grid on
    hold on
    
    % Plot initial, final points
    plot3(y(1,1), y(1,2), y(1,3), 'gx')
    plot3(y(end,1), y(end,2), y(end,3), 'rx')

    % process primaries
    P1.r = [-sys.mu, 0, 0];
    P1.rad = P.Results.p1Size;
    P2.r = [1 - sys.mu, 0, 0];
    P2.rad = P.Results.p2Size;
    [X, Y, Z] = sphere();
    
    % Plot primaries 
    if strcmp(P.Results.pPrim, 'both')
        surf(X.*P1.rad + P1.r(1), Y.*P1.rad + P1.r(2), Z.*P1.rad + P1.r(3))
        surf(X.*P2.rad + P2.r(1), Y.*P2.rad + P1.r(2), Z.*P2.rad + P2.r(3))
        legend('Traj','init','final','P1','P2')
    elseif strcmp(P.Results.pPrim, 'P1')
        surf(X.*P1.rad + P1.r(1), Y.*P1.rad + P1.r(2), Z.*P1.rad + P1.r(3))
        legend('Traj','init','final','P1')
    elseif strcmp(P.Results.pPrim, 'P2')
        surf(X.*P2.rad + P2.r(1), Y.*P2.rad + P1.r(2), Z.*P2.rad + P2.r(3))
        legend('Traj','init','final','P2')
    end
    axis equal

    varargout = {h};
end

