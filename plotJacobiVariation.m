function varargout = plotJacobiVariation(t,X,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTJACOBIVARIATION Helpful function to identify inaccuracies in numerical
%intergration via the Jacobi Constant
%   Inputs:
%       t -     [fl]  [mx1] Time vector
%       X -     [fl]  [mx6] State data in X,Y,Z,Xd,Yd,Zd format
%       mu -    [fl]  Mass ratio of the system 
%   Outputs (via varargout):
%       C - [mx1] vector of Jacobi constants
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C = calcJacobiConst(X(:,1:3), X(:,4:6), mu);

    figure
    subplot(2,1,1)
    plot(t, C)
    ylabel('C_j [--]')
    ylim([min(C) max(C)])
    xlabel('Nondim. Time [--]')
    title('C_j vs. Time')
    subplot(2,1,2)
    plot(t, C-C(1))
    ylabel('C_j [--]')
    xlabel('Nondim. Time [--]')
    title('Variation from init C_j')
    yline(0,'k--')

    varargout = {C};
end

