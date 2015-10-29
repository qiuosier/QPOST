%% EXAMPLE: Simulate the System with a Series of Operations
% This example shows how to use the QP_SimulateOperations function to 
%   simulate a power system with a series of Operations.
%
% Version $\Delta$. Qiu Qin, December 12, 2014. All Rights Reserved.

%% Load Data File
% Run the data file to load power system data. Please refer to the data
% file for detail description of the data format.
clear;
run('QP_9B_data3m9b.m');
%%
% Load constant variables. Theses constant variables are used to refer to a
% line in the 9 bus system.
QP_9B_LineDefinition;

%%
% Rename the data matrices.
generator = mac_con;
busData = bus(:,:);
lineData = line(:,:);

%% Specify the Operations
% The operations of the system during the simulation are specified in a
%   matrix. Each row of the matrix represent a system configuration.
% The columns are defined as:
%
% * Column 1: Simulation Time
% * Column 2: Line Fault
% * Column 3: Line Fault Location
% * Column 4: Line Removal
%
% Using the following specifications, the system will be simulated for 18
% cycles at normal operating condition, then fault-on configuration for 6
% cycles and a misoperated configuration for another 6 cycles, following by
% a post-fault operating condition for 90 cycles.
operationData = [
    12/60,  NO_FAULT,     0,    NO_REMOVAL;
    30/60,  LINE57,     0.5,    NO_REMOVAL;
    90/60,  NO_FAULT,     0,    LINE57;
    ];
%%
% The function will simulate the system for each configuration
% sequencially.
[ t, delta, omega ] = ...
    QP_SimulateOperation( operationData, busData, lineData, generator);

%%
% The simulation results are as follows:
figure
plot(t,delta);
xlabel('t');
ylabel('\delta');
title('Rotor Angle');
figure
plot(t,omega);
xlabel('t');
ylabel('\omega');
title('Rotor Speed');