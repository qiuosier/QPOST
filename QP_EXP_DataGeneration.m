%% EXAMPLE: Simulate the 9-Bus System
% This example shows how to simulate the 9-bus power system and obtained
% sampled data. Both generator dynamics and electrical dynamics are
% simulated.
% The QP_9B_RunModel function is specifically designed for the 9-bus
% system.
%
% Version $\Delta$. Qiu Qin, December 11, 2014. All Rights Reserved.

%% The Operation Sequence of the Simulation
% The system will be simulated as 5 sequential segments:
% 
% # Warm Up, system operates at normal condition but data are not recorded
% # Pre-Fault, system operates at normal condition
% # Fault-On, system operates with a line fault
% # Misoperated, system operates with a fault and a line removed
% # Post-Fault, system operates with the faulted line removed
%
% Load constant variables. Theses constant variables are used to refer to a
% line in the 9 bus system.
clear;
QP_9B_LineDefinition;
%%
% Specify the simulation configurations:

% The line short to ground
faultyLine = LINE57;
% Fault distance from the bus
faultDistance = 0.5;
% The mis-operated line
removeLine = LINE78;

%%
% Specify the time for each segment in seconds
WarmUpTime = rand(1) + 1; % This time can be a random number
PreFaultTime = 18/60;
FaultOnTime = 6/60;
MisOptTime = 6/60;
PostFaultTime = 30/60;

%% Run the Simulink Model and Record Data
% Run the simulink model and generate raw data.
% Please see the function file for details of the variables
[ MS_delta, MS_w, ES_input, ES_output, ES_state ] = ...
    QP_9B_RunModel( faultyLine, faultDistance, removeLine,...
    WarmUpTime, PreFaultTime, FaultOnTime, MisOptTime, PostFaultTime);

%% Resample Data
% The data recorded above will generally have variable step size. 
% QP_PrepareData can be used to resample data with a fixed sampling
% interval.
% Specify the number of samples per cycle
SamplePerCycle = 120;
[ SP_input, SP_output, ~, ~ ] =...
    QP_PrepareData( ES_input, ES_output, SamplePerCycle );
%%
% Note that SP_input and SP_output are input/output data for the 
%   state space model.
% See the QP_9B_StateSpaceEquation.pdf file for state space equations.