function [yVal,isterminal,direction] = CR3BPYIntDYPosEv(t,y)
%CR3BPYINTDYPOSEV Event function to capture the Y-Plane crossover (y = 0)
%   of the CR3BP with derivitave ydot > 0.
    yVal       = y(2); % Y-Plane crossover point that must be zero 
    isterminal = 1;    % Halt the integration when this condition is reached
    direction  = 1;    % This zero must be approached while ydot is increasing
end

