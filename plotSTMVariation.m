function varargout = plotSTMVariation(t,PHI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTSTMVARIATION Plots State Transition Matrix (STM) determinant variation
%from unity to assist with debugging and understanding of integration
%tolerances.
%   Inputs:
%       t -     [fl]  [mx1] Time vector
%       PHI -   [fl]  [mx36] ROW MAJOR form of the 6x6 STM Matrix
%   Outputs (via varargout):
%       determinant - [mx1] vector of STM determinants
%   Author:
%       Tyler Mason, tyma4153@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [nDataPts, ~] = size(PHI);
    determinant = zeros(nDataPts,1);
    for i = 1:nDataPts
        % Reconstruct STM from row major to 6x6
        STM = [PHI(i,1:6);
               PHI(i,7:12);
               PHI(i,13:18);
               PHI(i,19:24);
               PHI(i,25:30);
               PHI(i,31:36)];
        % Calculate difference between unity
        determinant(i) = det(STM);
    end
    
    % Plot determinant and difference
    figure
    subplot(2,1,1)
    plot(t, determinant)
    ylabel('Det(STM) [--]')
    xlabel('Nondim. Time [--]')
    title('STM Determinant vs. Time')
    subplot(2,1,2)
    plot(t, determinant - 1)
    ylabel('Variation [--]')
    xlabel('Nondim. Time [--]')
    title('Variation from Unity (det(STM) - 1)')
    yline(0,'k--')
    
    varargout = {determinant};
end

