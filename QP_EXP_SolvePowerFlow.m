%% EXAMPLE: Solve Power Flow
% This example shows how to use the QP_SolvePowerFlow function to solve the
% power flow.
%
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.

%% Load Data File
% Run the data file to load power system data. Please refer to the data
% file for detail description of the data format.
clear;
run('QP_9B_data3m9b');

%%
% Three matrix describing the power system are stored in the data file.
% In this file, they are named as bus, line and mac_con.
% To solve the power flow, only bus data and line data are needed.
%
% The bus data matrix consist of all information about the load flow.
% In the process of solving power flow, part of the bus matrix will be
% re-calculated. 
% For swing bus, only the bus voltage and angle are used. 
%   Generator output active and reactive power are solved by power flow,
%   regardless what value are specified in the data matrix
% For PV bus, only active power and voltage are used. 
%   Bus voltage angle and reactive power are solved by power flow,
%   regardless what values are specified in the data matrix
% For PQ bus, only active power and reactive power are used. 
%   Bus voltage and angle are solved by power flow,
%   regardless what values are specified in the data matrix

%% Solve the power flow and display results
% The following shows an example of solving power flow.
disp('Power Flow Solution for Normal Operating Condition');
solvedBus = QP_SolvePowerFlow( bus, line, 1 );

%% Modify the data to remove a line and re-solve the power flow
% The following shows an example of removing a line from the system and
% re-solving the power flow. Note that removing a line from the system
% creates a post-fault operating condtion.
modifiedLine = QP_RemoveLine( line, 5 );
disp('Power Flow Solution for Post-fault Operating Condition');
modifiedBus = QP_SolvePowerFlow( bus, modifiedLine, 1 );