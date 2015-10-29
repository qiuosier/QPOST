function [ t, delta, omega, varargout ] = ...
    QP_SimulateSystem( busData, lineData, generator, ...
    faultyLine, faultDistance, removeLine, simulateTime, ...
    delta0, omega0)
% Version $\Delta$. Qiu Qin, December 10, 2014. All Rights Reserved.
%
% This function simulates the system using classical generator model
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
% NOTE: line fault and line removal occurs at the beginning of simulation.
% simulationTime: Length of the simulation in seconds
% delta0: Initial rotor angle for simulation (rad)
% omega0: Initial rotor speed for simulation (pu)
%   Norminal omega is 1 pu.
% If delta0 and omega0 are not specified, the pre-fault equilibrium point
%   will be used.
%
% For the simulation parameters,
%   if an stable equilibrium point does not exist for the configuration, 
%   the pre-fault equilibrium point will be used for P and E
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
% OPTIMIZED FOR 9 BUS SYSTEM
%   9 bus system power flow has been solved and store in data files
%   this function will load the power flow solution from the file 
%   if input data match the 9 bus system data stored in QP_9B_Data.mat


if nargin == 7 || isempty(delta0)
    [ ~, ~, ~, delta0 ] = ... 
        QP_SimParameters( busData, lineData, generator, ...
        0, 0, 0 );
    omega0 = ones(size(delta0));
end

% Convert delta and omega to column vector
if size(delta0,1) == 1
    delta0 = delta0';
end
if size(omega0,1) == 1
    omega0 = omega0';
end

% QP_Constants defines constant variables used in this program.
QP_Constants;
% Inertia of Generator
H = generator(:,GEN_H);
% Damping coefficient
D = generator(:,GEN_D);
% Number of Generators
N_gen = size(H,1);
% Parameters for simulation
[ E, P, Y, ~ ] = ... 
    QP_SimParameters( busData, lineData, generator, ... 
    faultyLine, faultDistance, removeLine );

odefun = @(t,X)QP_ClassicalModel(P, E, Y, D, H, X);
tspan = [0, simulateTime];
X0 = [delta0;omega0];
options = odeset('RelTol',1e-3);
[t, X] = ode45(odefun,tspan,X0,options);

delta = X(:,1:N_gen);
omega = X(:,N_gen+1:end);

if nargout == 3
    varargout = [];
end
end

