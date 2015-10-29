%% EXAMPLE: Multiple Model Filtering for 9-Bus System
% This example shows how to implement a simple multiple model filtering
% algorithm for the 9-bus system.
%
% Version $\Delta$. Qiu Qin, December 11, 2014. All Rights Reserved.

%% Data Generation
% Please refer to the example of data generation for details.
%
% Specify the configurations:

clear;
% Load constant variables
QP_9B_LineDefinition;
% Specify the simulation configuration
% The line short to ground
faultyLine = LINE57;
% Distance from the bus
faultDistance = 0.5;
% The mis-operated line
removeLine = LINE78;
% Specify the time in seconds
WarmUpTime = 5.111; % This time can be a random number
PreFaultTime = 18/60;
FaultOnTime = 6/60;
MisOptTime = 6/60;
PostFaultTime = 30/60;

%%
% Run the simulink model and generate raw data:
[ MS_delta, MS_w, ES_input, ES_output, ES_state ] = ...
    QP_9B_RunModel( faultyLine, faultDistance, removeLine,...
    WarmUpTime, PreFaultTime, FaultOnTime, MisOptTime, PostFaultTime);
%%
% Resample data with a fixed sampling time
SamplePerCycle = 24;
[ SP_input, SP_output, ~, ~ ] =...
    QP_9B_PrepareData( ES_input, ES_output, SamplePerCycle );

% SP_input and SP_output are input/output data for the state space model.
% See the QP_9B_StateSpaceEquation.pdf file for state space equations.
% Note that it may not be necessary to generate data for each run.

%% Calculate State Space Matrices for filters
% Here for the 9-bus system, the state space matrices are calculated based
% on equations given in QP_9B_StateSpaceEquation.pdf.
%
% Load 9-bus system data:
run('QP_9B_data3m9b.m');
generator = mac_con;
QP_Constants;
QP_9B_LineDefinition;
%%
% System frequency is assumed to be 60Hz.
f = 60;
w = 2 * pi * f;
%%
% Obtain Line Parameters, i.e. R, L, C
SC = line(:,LINE_B) / w; %#ok<*NODEF>
SR = line(:,LINE_R);
SL = line(:,LINE_X) / w;
%%
% Loads are modeled as constant impedances. 
% The following are Load impedances
LR = [0.6836;1.0255;0.9194];
LL = [0.2735;0.3418;0.3218] / w;
%%
% Calculate state space model for normal operating condition
CN5 = (SC(LINE75) + SC(LINE75)) / 2;
A_Normal = [...
    -SR(LINE75)/SL(LINE45), 0, 0, -1/SL(LINE75);
    0, -SR(LINE45)/SL(LINE45), 0, -1/SL(LINE45);
    0, 0, -LR(1)/LL(1), 1/LL(1);
    1/CN5, 1/CN5, -1/CN5, 0;
    ];
B_Normal = [0, 1/SL(LINE75); 1/SL(LINE45), 0; 0, 0; 0, 0];
C_Normal = [1, 0, 0, 0; 0, 1, 0, 0];
% Calculate state space model for fault on Line 5-7
CF5A = SC(LINE75)/4 + SC(LINE45)/2;
A_F75 = [...
    -SR(LINE75)/SL(LINE75), 0, 0, 0, 0;
    0, -SR(LINE75)/SL(LINE75), 0, 0, -1/SL(LINE75)/2;
    0, 0, -SR(LINE45)/SL(LINE45), 0, -1/SL(LINE45);
    0, 0, 0, -LR(1)/LL(1), 1/LL(1);
    0, -1/CF5A, 1/CF5A, -1/CF5A, 0;
    ];
B_F75 = [0, 1/SL(LINE75)/2; 0, 0; 1/SL(LINE45), 0; 0, 0; 0, 0];
C_F75 = [1, 0, 0, 0, 0; 0, 0, 1, 0, 0];
% Fault on Line 4-5
CF5B = SC(LINE75)/2 + SC(LINE45)/4;
A_F45 = [...
    -SR(LINE75)/SL(LINE75), 0, 0, 0, -1/SL(LINE75);
    0, -SR(LINE45)/SL(LINE45), 0, 0, 0;
    0, 0, -SR(LINE45)/SL(LINE45), 0, -1/SL(LINE45)/2;
    0, 0, 0, -LR(1)/LL(1), 1/LL(1);
    1/CF5B, 0, -1/CF5B, -1/CF5B, 0;
    ];
B_F45 = [0, 1/SL(LINE75); 1/SL(LINE45)/2, 0; 0, 0; 0, 0; 0, 0];
C_F45 = [1, 0, 0, 0, 0; 0, 1, 0, 0, 0];
% Discretize
sysd = c2d(ss(A_Normal,B_Normal,C_Normal,0), 1/60/SamplePerCycle);
A_Normal = sysd.A;
B_Normal = sysd.B;
C_Normal = sysd.C;
sysd = c2d(ss(A_F75,B_F75,C_F75,0), 1/60/SamplePerCycle);
A_F75 = sysd.A;
B_F75 = sysd.B;
C_F75 = sysd.C;
sysd = c2d(ss(A_F45,B_F45,C_F45,0), 1/60/SamplePerCycle);
A_F45 = sysd.A;
B_F45 = sysd.B;
C_F45 = sysd.C;
%% Multiple Filtering
% Determine Q and R for kalman filter
Q = 0.05 * eye(5);% * diag([MAG_state,MAG_state(14)]);
R = 0.02 * eye(2);% * diag(MAG_output);
% No interacting between filters
T = eye(3);
% CAUTION: The states for different filters in this example are not in the
% same order.

% Specify input/output data, 
% Output data for the three filters are different
u = SP_output.Data(:,[1,5])';
dt = SP_output.Time(2)-SP_output.Time(1); 
yN = [...
    SP_output.Data(:,6)'- SC(LINE75)/2 * ...
    (SP_output.Data(:,5)' - [0, SP_output.Data(1:end-1,5)'])/dt;
    SP_output.Data(:,3)'- SC(LINE45)/2 * ...
    (SP_output.Data(:,1)' - [0, SP_output.Data(1:end-1,1)'])/dt;
    ];
yF75 = [...
    SP_output.Data(:,6)'- SC(LINE75)/4 * ...
    (SP_output.Data(:,5)' - [0, SP_output.Data(1:end-1,5)'])/dt;
    SP_output.Data(:,3)'- SC(LINE45)/2 * ...
    (SP_output.Data(:,1)' - [0, SP_output.Data(1:end-1,1)'])/dt;
    ];
yF45 = [...
    SP_output.Data(:,6)'- SC(LINE75)/2 * ...
    (SP_output.Data(:,5)' - [0, SP_output.Data(1:end-1,5)'])/dt;
    SP_output.Data(:,3)'- SC(LINE45)/4 * ...
    (SP_output.Data(:,1)' - [0, SP_output.Data(1:end-1,1)'])/dt;
    ];

%% Initialization
N_output = 2;           % Total number of outputs
N_sample = size(yN,2);  % Total number of samples
N_state = 5;            % Number of states, including states during fault
N_filter = 3;           % Number of filters
% SystemState store the diagnosis output
SystemState = zeros(1,N_sample);

x0 = rand(18,1);
u0 = ones(3,1);
mu0 = [0.94, 0.03, 0.03];
mu = zeros(N_filter,N_sample);
mu(:,1) = mu0;
xk = ones(N_state, N_filter);
Pk = zeros(N_state, N_state, N_filter);
for i = 1:N_filter
    Pk(:,:,i) = 1 * eye(N_state);
end

% Px is not used, forget it
Px = Pk;

for k = 1:N_sample-1
    %% Interaction / Re-initialization
    x0 = xk; %zeros(N_state, N_filter);
    P0 = Pk; %zeros(N_state, N_state, N_filter);
    for j = 1:N_filter
        % For states in all filters
        % Predicted mode probability
        for i = 1:N_filter
            mu(j,k+1) = mu(j,k+1) + T(i,j) * mu(i,k);
        end
        % NO Mixing estimate
        % NO Mixing convariance
    end
    
    %% Model conditional filtering
    L = zeros(1,N_filter);
    uk = u(:,k);
    % j is used to indicate different filters
    % Normal
    j = 1;
    yk = yN(:,k+1);
    [xe, ~, Pk(1:N_state-1,1:N_state-1,j), L(j)] = ...
        QP_END_Kalman(A_Normal, B_Normal, C_Normal, ...
        P0(1:N_state-1,1:N_state-1,j), ...
        Q(1:N_state-1,1:N_state-1), R, ...
        x0(1:N_state-1,j), uk, yk, [] );
    xk(1:4,j) = xe;
    % Faulted Line 75
    j = 2;
    yk = yF75(:,k+1);
    SSF = 1:5;
    [xe, ~, Pk(SSF,SSF,j), L(j)] = ...
        QP_END_Kalman(A_F75, B_F75, C_F75, ...
        P0(SSF,SSF,j), Q(SSF,SSF), R, x0(SSF,j), uk, yk, [] );
    xk(SSF,j) = xe;
    % Faulted Line 45
    j = 3;
    yk = yF45(:,k+1);
    SSF = 1:5;
    [xe, ~, Pk(SSF,SSF,j), L(j)] = ...
        QP_END_Kalman(A_F45, B_F45, C_F45, ...
        P0(SSF,SSF,j), Q(SSF,SSF), R, x0(SSF,j), uk, yk, [] );
    xk(SSF,j) = xe;
    % Removed Line 75
    % Removed Line 45
    %% Mode probability update
    mu_sum = 0;
    for i = 1:N_filter
        mu_sum = mu_sum + mu(i,k+1) * L(i);
    end
    for j = 1:N_filter
        mu(j,k+1) = mu(j,k+1) * L(j) / mu_sum;
    end
    
    %% Fault Decision
    H = 0.1;
    idx = find(mu(:,k+1) == max(mu(:,k+1)),1);
    if idx > 1 && mu(idx,k+1) < H
        idx = 1;
    end
    SystemState(k+1) = idx;
    
    %% set minimum mu, this is not necessary
    for j = 1:N_filter
        if mu(j,k+1) < 1e-5
            mu(j,k+1) = 1e-5;
        end
    end
end
plot(SP_output.Time,mu');
