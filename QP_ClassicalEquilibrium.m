function [ delta, Pm, E0, busData, Y ] = ...
    QP_ClassicalEquilibrium( busData, lineData, generator, outputResult )
% Version $\Epsilon$. Qiu Qin, March 22, 2015. All Rights Reserved.
%
% WARNING: The symbolic toolbox in MATLAB 2014 becomes VERY SLOW
%   The reason is unknown.
%
% This function calculate the steady state values for rotor angle,
%   mechanical power, and internal voltage of generators in power 
%   system using classical generator models.
%
% This function requires the following files:
%   QP_Constants, script defining constant variables
%   QP_GenVoltage, the function to calculate generator interal voltage
%   QP_ReducedYMatrix, calculate the reduced network admittance matrix
%   QP_RegularYMatrix, calculate the regular admittance matrix
%   QP_SolvePowerFlow, the function to solve the power flow
%
%
% INPUT ARGUMENTS:
% busData, lineData and generator are matrices decribing the power system
% The data format is the same as PST by Prof. Chow at RPI
% Please refer to the data file for details.
% outputResult indicate whether the results are output to MATLAB command
%   window, 1 = show results, 0 = do not results
%
% OUTPUT ARGUMENTS:
% At steady state, rotor speed is 1 pu.
% delta  : rotor angle in rad
% Pm     : mechanical power in pu
% E0     : generator internal voltage (complex number) in pu 
% busData: solved power flow
% Y      : reduced network admittance matrix
%
% EXAMPLE:
% Assuming QP_9B_data3m9b is loaded.
%
% [delta, Pm, E0, busData, Y ] = ...
%    QP_ClassicalEquilibrium( bus, line, mac_con, 1 )
%
% calculates the equilibrium point and show the results.
% (including power flow solutions).
%
% FUTURE DEVELOPMENTS:
%   Check if power flow has a solution.
%   Check post fault case and compare with Transient EQ function



% Check the number of input arguments
if nargin < 4
   outputResult = 0;
end

% Initialization
% Load Constants
QP_Constants;
% Number of Generators
N_gen = size(generator,1);

% Solve power flow
busData = QP_SolvePowerFlow(busData, lineData, outputResult);
% Generator Internal Voltage
E0 = QP_GenVoltage(generator,busData);
% Rotor angle = Voltage angle in classical model
delta = angle(E0);

% The following Code use reduced admittance matrix
% Reduced Admittance Matrix
Y = QP_ReducedYMatrix(busData, lineData, generator);
% Calculate initial value of Pm 
%   Ref: Anderson's book pp. 37, eq. 2.57
EM = abs(E0);
Pm = zeros(N_gen,1);
for i = 1:N_gen
    Pm(i) = EM(i)^2 * real(Y(i,i));
    for j = 1:N_gen
        if j == i 
            continue;
        end
        Pm(i) = Pm(i) + EM(i) * EM(j) * abs(Y(i,j)) * ...
            cos(angle(Y(i,j)) - delta(i) + delta(j));
    end
end


end

