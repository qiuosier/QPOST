function [ t, delta, omega] = ...
    QP_SimulateOperation( operationData, busData, lineData, generator)
% Version $\Delta$. Qiu Qin, December 10, 2014. All Rights Reserved.
%
% This function simulate the system using classical generator model
% The system will operate according to the sequence defined in the 
%   operationData.
% 
% This function requires the following files:
%   QP_ClassicalEquilibrium: 
%       Calculates the equilibrium point of a power system.
%   QP_ClassicalModel.m, clculate the rotor angle and rotor speed 
%       deviations in  classical generator model
%   QP_Constants.m, script defining constant variables.
%   QP_GenVoltage, the function to calculate generator interal voltage
%   QP_ReducedYMatrix, calculate the reduced network admittance matrix
%   QP_RegularYMatrix, calculate the regular admittance matrix
%   QP_RemoveLine, Removes a line from the system data
%   QP_ShortToGroundFault, 
%       Creates new data matrices by breaking the faulted line
%   QP_SimParameters, calculates the parameters for simulation
%   QP_SimulateSystem, 
%       simulates the system using classical generator model
%   QP_SolvePowerFlow, solves the power flow
%
% INPUT ARGUMENTS:
% busData, lineData and generator are matrices decribing the power system
%   The data format is the same as PST by Prof. Chow at RPI
%   Please refer to the data file for details.
% operationData is a matrix defined as follows:
%   Each row define a period for a particular operating mode
%   The matrix should consist of four columns:
%   Column 1: (simulationTime) Duration of the period.
%   Column 2: (faultyLine) Specify a line fault, if any. Use 0 if no fault.
%   Column 3: (faultDistance) Specify the location of the fault.
%               This is the distance from a bus, which is the "FROM BUS"
%               defined in the busData. This parameter will be no effect
%               when the corresponding column 2 is set to 0.
%   Column 4: (removeLine) Specify a line to be removed during the peroid,
%               if any. Use 0 if no line to be removed.
%
% OUTPUT ARGUMENTS:
% t: Time stamp for each simulation step, a column vector
% delta: rotor angle for each generator (rad)
% omega: rotor speed for each generator (pu)
%   Norminal rotor speed is 1 pu
% t, delta and omega are tall matrices.
% Each row of t, delta, and omega corresponds to each other.
% NOTE: the time step size is NOT fixed.
%
% EXAMPLE:
% Please see QP_EXP_SimulateOperation for an example


% Number of rows (periods) in operationData
N_row = size(operationData, 1);
% Initialize variables
t = [];
delta = [];
omega = [];
% Simulation
for i = 1:N_row
    if isempty(t)
        [tX, deltaX, omegaX] = ...
            QP_SimulateSystem( busData, lineData, generator, ...
            operationData(i,2), operationData(i,3), ...
            operationData(i,4), operationData(i,1), ...
            [], []);
    else
        [tX, deltaX, omegaX] = ...
            QP_SimulateSystem( busData, lineData, generator, ...
            operationData(i,2), operationData(i,3), ...
            operationData(i,4), operationData(i,1), ...
            delta(end,:)', omega(end,:)');
        tX = tX + t(end);
    end
    t = [t;tX]; %#ok<*AGROW>
    delta = [delta;deltaX];
    omega = [omega;omegaX];
end
end

