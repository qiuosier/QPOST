function [ SP_input, SP_output, PM_output, PM_time ] =...
    QP_PrepareData( ES_input, ES_output, SamplePerCycle, ...
    PhasorPerCycle, WindowSize)
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function resample the raw sinusoid data and 
%   provide sampled data at a fixed interval.
% Phasor data is also calculated based on the resampled data.
%
% Note that the data generated by simulation may have variable steps. 
%
% The phasor data is calculated based on the nonrecursive method 
%   described in Chapter 2 (pp. 30-32) of:
%   Synchronized Phasor Measurements and Their Applications
%   A.G. Phadke and J.S. Thorp
%
% INPUT VARIABLES:
% ES_input: Data to be resampled, but for which no Phasor is calculated
% ES_output: Data to be resampled, and for which Phasors are calculated
% SamplePerCycle: Number of samples per cycle in the resampled data
% 	Frequency of 60 Hz is used in this function, i.e. 60 cycles/second
% PhasorPerCycle: Number of Phasors per cycle in the Phasor data
% WindowSize: The number of samples used to calculate each Phasor
%
% OUTPUT VARIABLES:
% SP_input: Resampled sinusoid inputs corresponding to ES_input
% SP_output: Resampled sinusoid outputs corresponding to ES_output
% PM_output: Phasor outputs estimated based on SP_output
% PM_time: the time stamps corresponding the values in PM_output
% Note that SP_input and SP_output are timeseries. PM_output is a matrix.
% PM_output is tall matrix with each row as the Phaosrs at a time.


% Set default values for phasor per cycle and window size, if not supplied
if nargin == 3
    PhasorPerCycle = 4;
    WindowSize = 24;
end

%% Resample Data
% Data Length (seconds)
EndTime = ES_output.Time(end);
% 
ResampleTime = 0:1/60/SamplePerCycle:EndTime;
%
SP_output = resample(ES_output, ResampleTime);
SP_input = resample(ES_input, ResampleTime);

%% Phasor Estimation
PM_time = 0:1/60/PhasorPerCycle:EndTime;
% Window size
N_Window = WindowSize;
% Total number of phasors
N_Phasor = floor(EndTime * 60 * PhasorPerCycle) + 1;
PM_time = PM_time(1:N_Phasor);
% Number of outputs
N_Output = size(ES_output.Data,2);
% Phasor estimates
PM_output = zeros(N_Phasor,N_Output);
for i = 1:N_Output
    XN = SP_output.Data(:,i);
    Theta = 2 * pi / SamplePerCycle;
    % Calculate S
    S = zeros(N_Window, 2);
    for p = 1:N_Window
        S(p,1) = cos((p-1)*Theta);
        S(p,2) = -sin((p-1)*Theta);
    end
    %S = sqrt(2) * S;
    X_Estimated = ones(N_Phasor,1);
    for k = 2:N_Phasor
        j = floor((k-1) * SamplePerCycle/PhasorPerCycle);
        if j < N_Window
            continue;
        end
        Window = XN(j - N_Window + 1 : j);
        V = (S' * S)^-1 * S' * Window;
        % Calculate Phasor, angle shift is added
        X_Estimated(k) = (V(1) + 1i*V(2)) * ...
            exp(-1i * 2 * pi * k / PhasorPerCycle);
    end
    PM_output(:,i) = X_Estimated;
end


end

