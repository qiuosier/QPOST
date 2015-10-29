function [ MS_delta, MS_w, ES_input, ES_output, varargout] = ...
    QP_9B_RunModel( faultyLine, faultDistance, removeLine,...
    WarmUpTime, PreFaultTime, FaultOnTime, MisOptTime, PostFaultTime, ...
    PWR)
% Version $\Epsilon$. Qiu Qin, Feburary 5, 2015. All Rights Reserved.
%
% Feb 5: Added process noise
% 
% This function runs simulink model for the 9 bus system for the following
% operations:
%   1. At the beginning,
%       the system starts a warm-up period around the equilibrium point.
%   2. After a period of WarmUpTime,
%       the system runs at normal operating condition.
%   3. After a period of PreFaultTime, 
%       a short to ground fault occurs.
%   4. After a period of FaultOnTime,
%       removeLine is removed while the fault still exist
%   5. After a period of MisOptTime,
%       the fault is cleared and only faultyLine is removed.
% 
% This function requires the following files:
%   QP_9B_SimModel.slx, the Simulink model.
%   QP_9B_data3m9b.m, the data file.
%   QP_9B_SSModel.m, the function to calculate the state space model.
%   QP_9B_SimParameter Folder, the parameters for Simulink model.
%       This should be eliminated in the future.
%   QP_9B_StoreState.m, script used to record the data after running the
%       Simulink model.
%   QP_9B_LineDefinition.m, script defining constant variables for
%       transmission lines in the 9-bus system.
%   QP_Constants.m, script defining constant variables.
% 
% INPUT ARGUMENTS:
% 
% faultyLine is the line that shorted to ground
% For the 9-bus system, if QP_9B_LineDefinition is loaded,
%   its value can be set as LINEXX, 
%   where XX is the bus numbers at two ends of the line.
% 
% faultDistance specifies the location of the fault.
%   It is the percentage of line length from the "FROM BUS" (which is
%   defined in the data file).
% 
% removeLine is the line that is removed from the system,
%   which can be considered as mis-operated.
% 
% System Operation Times:
% 
% The following times are specified in SECONDS.
% Please note that the times are DURATION of the period.
% 
% WarmUpTime is a preiod to run the simulation without recording data.
% The purpose of using this is to let the system settle to steady state.
% Warm up preiod can be a relatively large random number
% 
% PreFaultTime is the time duration that the system is running normally.
% FaultOnTime is the time duration that the system is running
%   with _faultyLine_ shorted
% 
% MisOptTime is the time duration that the system is running with
%   both _faultyLine_ shorted and _removeLine_ removed,
%   which can be considered as the mis-operated preiod.
% 
% PostFaultTime is the time duration that the system is at post-fault
%
% PWR (Optional) is the noise power. This value is the ratio of noise power 
%   and signal power when the system is operating in normal mode.
% 
% OUTPUT ARGUMENTS:
% The following data are obtained by Simulink
% MS_delta: Rotor angle
% MS_w:     Rotor speed
% Input/Output/State variables for the electrical network model:
% ES_input: Raw sinusoid input samples for voltages at the 3 generator, 
%               i.e. e
% ES_output: Raw sinusoid output samples for bus voltages and currents at 
%               bus 4, 7, 9
% varargout is an optional output.
% When varargout is specified, it will be ES_state.
% ES_state: Raw sinusoid waveforms for states
% CAUTION:
% ES_state will always return values for 18 states, as defined in the 
%   state space equation of the normal operating condition.
%   When there is a fault on a line, the state corresponds to the current
%   at the faulted line will be the current from the generator to the
%   faulted point. When a line is removed, the corresponding state will
%   be set to 0.
% Please refer to QP_9B_StateSpaceEquation.pdf for details of the
%   state space model
% 
% EXAMPLE:
% Assuming QP_9B_LineDefinition is loaded, LINEXX can be used to refer a
%   transmission line in the system.
% 
% [ MS_delta, MS_w, ES_input, ES_output] = ...
%    QP_9B_RunModel( LINE57, 0.5, LINE78, 1.0, 0.3, 0.1, 0.1, 0.5);
% 
% Runs a simulation of 1 second (excluding the warm-up time).
% At t = 0.3, LINE57 shorted to ground at midpoint.
% At t = 0.4, LINE78 is removed while LINE57 is still shorted.
% At t = 0.5, LINE57 is removed and LINE78 is reconnected.
% 
% HOW IT WORKS:
% This function loads the simulink model named QP_9B_SimModel and run it
% for a few times with different configurations. The simulated results for
% each run is recoreded and concatenated as the output of this function.
% Parameters for the generator models are stored in 
%   QP_9B_SimParameter folder
% Parameters for the transmission network are calculated 
%   using QP_9B_SSModel


% PROGRAM START FROM HERE
% No process noise, if noise level is not specified
if nargin == 8
    PWR = 0;
end
% The signal power when system is normal
% State variables noise
statePWR = [
    1.0258
    0.9955
    1.0127
    1.0257
    1.0158
    1.0325
    1.3510
    0.9363
    1.0422
    0.7378
    0.5082
    0.3137
    1.5883
    0.8458
    0.7465
    0.8368
    0.5886
    0.2734] .^ 2 / 2;
% Output variables noise
outputPWR = [
    1.0258
    0.7378
    0.4579
    0.3006
    1.0257
    0.8463
    1.5883
    0.7434
    1.0325
    0.6141
    0.2376
    0.8368] .^ 2 / 2;

% State index
% This is the order index for the 6 line currents in the state space
%   equations.
StateIndex = [11, 14, 12, 17, 15, 18];
% Three zeros are added to the beginning of StateIndex to make it aligned
%   with the line number in the data file.
%   So, while current of line 57 is the 14th state in the state space 
%   equations, we can have: StateIndex(LINE57) = 14
StateIndex = [0, 0, 0, StateIndex];

% Load simulink model
load_system('QP_9B_SimModel');
SimOption = simset('SrcWorkspace','current');

%% Define variables
% Create time series
MS_delta = timeseries('delta');
MS_w = timeseries('omega');
ES_input = timeseries;
ES_output = timeseries;
ES_state = timeseries;
%% Calculated state space model for all cases
%EST_SSModel;
sys_Normal = QP_9B_SSModel( 0, 0, 0 );
A_Normal = sys_Normal.A;
B_Normal = sys_Normal.B;
C_Normal = sys_Normal.C;
sys_FaultOn = QP_9B_SSModel( faultyLine, faultDistance, 0 );
A_FaultOn = sys_FaultOn.A;
B_FaultOn = sys_FaultOn.B;
C_FaultOn = sys_FaultOn.C;
sys_MisOpt = QP_9B_SSModel( faultyLine, faultDistance, removeLine );
A_MisOpt = sys_MisOpt.A;
B_MisOpt = sys_MisOpt.B;
C_MisOpt = sys_MisOpt.C;
sys_PostFault = QP_9B_SSModel( 0, 0, faultyLine );
A_PostFault = sys_PostFault.A;
B_PostFault = sys_PostFault.B;
C_PostFault = sys_PostFault.C;
%% Load Data File
% Using the same data file format as Power System Toolbox by Chow
run('QP_9B_data3m9b');
generator = mac_con;
QP_Constants;
%% The following parameters are for simulink model
% Inertia of Generator
H = generator(:,GEN_H);
% Damping Coefficient
D = generator(:,GEN_D);
%% Initial states
w0 = zeros(3,1)+0.5/60; % Will be overrided?
RC_State0 = 0;
%% System Warm-up
% Use normal state space model
SysA = A_Normal;
SysB = B_Normal;
SysC = C_Normal;
XNoisePWR = statePWR * PWR;
YNoisePWR = outputPWR * PWR;

% The following data file contains Y matrix and delta angles
% IMPORTANT: w0 and delta0 will be overwrited
load(strcat('QP_9B_SimParameter\SimParaF0R0.mat'))
% Set initial rotor angle
E0_delta = delta0; %#ok<NODEF> delta0 is in the SimPara file
deltaState = delta0;
omegaState = w0;

%% Warm up Model
if WarmUpTime > 0
    % Run the model for a warm-up preiod, 
    %   so that the states will reach steady-state
    set_param('QP_9B_SimModel', 'StopTime', num2str(WarmUpTime));
    sim('QP_9B_SimModel',[],SimOption)
    % Save the final states
    % These two states are recorded as arrays
    deltaState = deltaValues;
    omegaState = omegaValues.Data(end,:);
    % States of the electric network are recorded as timeseries
    RC_State0 = RC_State.Data(end,:); %#ok<*NASGU,NODEF>
end
%% Prefault Model
if PreFaultTime > 0
    % Set system states
    delta0 = deltaState;
    w0 = omegaState;
    % Start Simulation
    set_param('QP_9B_SimModel', 'StopTime', num2str(PreFaultTime));
    sim('QP_9B_SimModel',[],SimOption);
    QP_9B_StoreState;
end

%% Fault On Model
if FaultOnTime > 0
    % Use fault on state space model
    SysA = A_FaultOn;
    SysB = B_FaultOn;
    SysC = C_FaultOn;
    XNoisePWR = [statePWR;statePWR(StateIndex(faultyLine))] * PWR;
    % Add a new state in fault-on model, the new state's intial value is
    %   the current of the line where the fault occured
    RC_State0 = [RC_State.Data(end,:), ...
        RC_State.Data(end,StateIndex(faultyLine))];
    % Load system parameters
    load(strcat('QP_9B_SimParameter\SimParaF', ...
        num2str(faultyLine),'R0.mat'))
    % Set system states
    delta0 = deltaState;
    w0 = omegaState;
    % Run fault-on system
    set_param('QP_9B_SimModel', 'StopTime', num2str(FaultOnTime));
    sim('QP_9B_SimModel',[],SimOption);
    % Separate the original 18 states and the addional state
    FT_State = RC_State.Data(end, 19);
    RC_State.Data = RC_State.Data(:,1:18);
    QP_9B_StoreState;
end

%% Misoperated Model
if MisOptTime > 0
    % Use fault on state space model
    SysA = A_MisOpt;
    SysB = B_MisOpt;
    SysC = C_MisOpt;
    % Remove a state in misoperated model
    RC_State0 = RC_State.Data(end,[1:StateIndex(removeLine)-1, ...
            StateIndex(removeLine)+1:18]);
    XNoisePWR = statePWR([1:StateIndex(removeLine)-1, ...
            StateIndex(removeLine)+1:end]) * PWR;
    % Add additional state if a fault exists
    if faultyLine > 0
        if FaultOnTime > 0
            RC_State0 = [RC_State0, FT_State];
        else
            RC_State0 = [RC_State0, ...
                RC_State.Data(end,StateIndex(faultyLine))];
        end
        XNoisePWR = [XNoisePWR;statePWR(StateIndex(faultyLine))];
    end
    % Load system parameters
    load(strcat('QP_9B_SimParameter\SimPara', ...
        'F',num2str(faultyLine),'R',num2str(removeLine),'.mat'))
    % Set system states
    delta0 = deltaState;
    w0 = omegaState;
    % Run fault-on system
    set_param('QP_9B_SimModel', 'StopTime', ...
        num2str(MisOptTime));
    sim('QP_9B_SimModel',[],SimOption);
    % Separate the original 18 states and the addional state
    if faultyLine > 0
        FT_State = RC_State.Data(end,18);
    end
    % Insert 0 values for the removed state
    RC_State.Data = [RC_State.Data(:,1:StateIndex(removeLine)-1), ...
            0 * RC_State.Data(:,1), ...
            RC_State.Data(:,StateIndex(removeLine):17)];
    QP_9B_StoreState;
end

%% Post Fault Model
if PostFaultTime > 0
    % Use post fault state space model
    SysA = A_PostFault;
    SysB = B_PostFault;
    SysC = C_PostFault;
    % Remove the state of faulted line
    RC_State0 = [RC_State.Data(end,1:StateIndex(faultyLine)-1),...
        RC_State.Data(end,StateIndex(faultyLine)+1:end)];
    XNoisePWR = statePWR([1:StateIndex(faultyLine)-1, ...
            StateIndex(faultyLine)+1:end]) * PWR;
    % Load system parameters
    load(strcat('QP_9B_SimParameter\SimParaF0R', ...
        num2str(faultyLine),'.mat'))
    % Set system states
    delta0 = deltaState;
    w0 = omegaState;
    % Run fault-on system
    set_param('QP_9B_SimModel', 'StopTime', ...
        num2str(PostFaultTime));
    sim('QP_9B_SimModel',[],SimOption);
    % Insert 0 values for the removed state
    RC_State.Data = [RC_State.Data(:,1:StateIndex(faultyLine)-1), ...
            0 * RC_State.Data(:,1), ...
            RC_State.Data(:,StateIndex(faultyLine):17)];
    QP_9B_StoreState;
end
if nargout == 4
    varargout = [];
elseif nargout == 5
    varargout{1} = ES_state;
end
close_system('QP_9B_SimModel',0);
end

