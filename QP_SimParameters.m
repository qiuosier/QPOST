function [ E, Pm, Y, delta0 ] = ... 
    QP_SimParameters( busData, lineData, generator, ... 
    faultyLine, faultDistance, removeLine )
% Version $\Delta$. Qiu Qin, December 10, 2014. All Rights Reserved.
%
% This funtion calculates the parameters for simulation.
% The parameters are based on classical generator model.
% For configurations without a power flow solution,
%   the pre-fault E, Pm and delta0 are returned.
%
% OPTIMIZED FOR 9 BUS SYSTEM
%   9 bus system power flow has been solved and store in data files
%   this function will load the power flow solution from the file 
%   if input data match the 9 bus system data stored in QP_9B_Data.mat
%   The filename for the pre-fault power flow solution is PFS9B_0L.mat
% 
% This function requires the following files:
%   QP_ClassicalEquilibrium: 
%       Calculates the equilibrium point of a power system.
%   QP_Constants.m, script defining constant variables.
%   QP_GenVoltage, the function to calculate generator interal voltage
%   QP_ReducedYMatrix, calculate the reduced network admittance matrix
%   QP_RegularYMatrix, calculate the regular admittance matrix
%   QP_RemoveLine, Removes a line from the system data
%   QP_ShortToGroundFault, 
%       Creates new data matrices by breaking the faulted line
%   QP_SolvePowerFlow, solves the power flow
%
% INPUT ARGUMENTS:
% busData, lineData and generator are matrices decribing the power system
%   The data format is the same as PST by Prof. Chow at RPI
%   Please refer to the data file for details.
% faultyLine: the line that is shorted to ground
% faultDistance: the location of the fault
%   It is the portion (between 0 and 1) of line length from the "FROM BUS" 
%   (which is defined in the data file).
% removeLine: the line that is removed from the system,
%   which can be considered as mis-operated.
% Use 0 if no fault or no line removal.
%
% OUTPUT ARGUMENTS:
% delta0 : initial rotor angle in rad
% Pm     : mechanical power in pu
% E      : generator internal voltage magnitude in pu 
% Y      : reduced network admittance matrix

% Check if the input data is the same as the 9 bus system
is9BusSystem = 0;
if exist('QP_9B_Data.mat','file')
    load('QP_9B_Data.mat');
    if isequal(busData, busData9B) && ...
            isequal(lineData, lineData9B) && isequal(generator, macData9B)
        is9BusSystem = 1;
    end
end

% Pre-fault equilibrium point
if is9BusSystem && exist('QP_9B_PowerFlow/PFS9B_0L.mat','file')
    load('QP_9B_PowerFlow/PFS9B_0L.mat');
else
    [delta0, Pm, E0, busData, ~] = ...
        QP_ClassicalEquilibrium(busData, lineData, generator, 0);
end

%% Remove a line
% Remove one line from the system
rLine = QP_RemoveLine( lineData, removeLine );
% Adjust the index of the faulty line
if faultyLine == removeLine
    idx = -removeLine;
elseif faultyLine > removeLine && removeLine ~= 0
    idx = -1;
else
    idx = 0;
end

%% Create a short circuit fault
% Generate a faulty point
% The faulty point is to be used as input for other functions
[ fBus, fLine, faultyPoint ] = ...
    QP_ShortToGroundFault(busData, rLine, faultyLine + idx, faultDistance);

% Calculate Ym matrix, considering the generator impedance
Y = QP_ReducedYMatrix(fBus, fLine, generator, faultyPoint);

%% Post-fault equilibrium point, if exist
if removeLine ~= 0 && faultyLine + idx == 0
    % No fault in the system
    % OPTIMIZED FOR 9 BUS SYSTEM
    %   The filename for post-fault parameters solution is PFS9B_XL.mat,
    %   where X is the removed line.
    filename = ...
        strcat('QP_9B_PowerFlow/PFS9B_',num2str(removeLine),'L.mat');
    if is9BusSystem && exist(filename,'file')
        load(filename);
    else
        [delta0, Pm, E0, ~, ~] = ...
            QP_ClassicalEquilibrium(fBus, fLine, generator, 0);
    end
end

E = abs(E0);

end

