function [ solvedBus ] = QP_SolvePowerFlow(busData, lineData, printResult)
% Version $\Epsilon$. Qiu Qin, March 22, 2015. All Rights Reserved.
%
% WARNING: The symbolic toolbox in MATLAB 2014 becomes VERY SLOW
%   The reason is unknown.
%
% This program solve the power flow using Newton-Raphson method.
%
% This function requires the following files:
%   QP_Constants.m, script defining constant variables.
%   QP_RegularYMatrix, calculate the regular admittance matrix
%
% INPUT ARGUMENTS:
% busData and lineData are matrices decribing the power system
%   The data format is the same as PST by Prof. Chow at RPI
%   Please refer to the data file for details.
% printResult is an optional argument. When set to nonzero value,
%   this function will display the power flow solution in the MATLAB
%   command window, otherwise nothing will be displayed.
%
% OUTPUT ARGUMENTS:
% solvedBus: the solved power flow results
% solvedBus has the same data format as bus
% 
% $\Epsilon$ Revision: The tolerance (stopping criteria) will based on the
%   number of buses in the system, smaller tolerance is set for smaller
%   systems.

if nargin < 3
   printResult = 0;
end
% The following file define the constant variables used in this function
QP_Constants;

%%
% Use vectors to store P, Q, voltage magnitude (V) and voltage angle (a)
P_v = busData(:,P_GEN) - busData(:,P_LOAD);
Q_v = busData(:,Q_GEN) - busData(:,Q_LOAD);
V_v = busData(:,VOLTAGE);
a_v = busData(:,ANGLE) * pi / 180;
%% 
% Number of bus
N_bus = size(busData,1);
% Find swing bus
[BusSW, ~] = find(busData(:,BUS_TYPE) == REF_BUS);
% Find PV bus
[BusPV, ~] = find(busData(:,BUS_TYPE) == PV_BUS);
% Find PQ bus
[BusPQ, ~] = find(busData(:,BUS_TYPE) == PQ_BUS);
%% Define symbolic variables
% Here _a_ is used to denote the angle $\theta$.
V = sym('V', [N_bus,1]);
a = sym('a', [N_bus,1]);
%%
% _EP_ and _EQ_ are used to denote the right hand side of equations (1) and
% (2).
EP = sym('EP',[N_bus,1]);
EQ = sym('EQ',[N_bus,1]);
%%
% Initialize _EP_ and _EQ_.
for i = 1:N_bus
    EP(i) = 0;
    EQ(i) = 0;
end
%% Construct equations and Jacobian for Newton-Raphson method
% Admittance Matrix
Y = QP_RegularYMatrix(busData, lineData);
G = real(Y);
B = imag(Y);
% Construct the right hand side of equations with symbolic variables
for i = 1:N_bus
    for j = 1:N_bus
        EP(i) = EP(i) + V(i) * V(j) * ...
            (G(i,j) * cos(a(i) - a(j)) + B(i,j) * sin(a(i) - a(j)));
        EQ(i) = EQ(i) + V(i) * V(j) * ...
            (G(i,j) * sin(a(i) - a(j)) - B(i,j) * cos(a(i) - a(j)));
    end
end
%%
% Take the equations with known $P_i$ and $Q_i$.
Equations = [EP([BusPV;BusPQ]);EQ(BusPQ)];
%%
% Define Jacobian
J = jacobian(Equations,[a([BusPV;BusPQ]); V(BusPQ)]);
%%
% Take the known values of $P_i$ and $Q_i$ 
PQ = [P_v([BusPV;BusPQ]);Q_v(BusPQ)];
%%
% Intial guesses for Newton-Raphson method
guess = [a_v([BusPV;BusPQ]);V_v(BusPQ)];
%% Newton-Raphson Method
% Initialize the change of angles and voltages, _daV_ as arbitary numbers.
daV = ones(size(Equations,1),1);
%%
% Set the tolerance for stopping condition. 
% A common stopping condition is to terminate if the norm of the 
% mismatch equations is below a specified tolerance.
if N_bus < 10
    tolerance = 1e-3;
elseif N_bus < 20
    tolerance = 1e-2;
elseif N_bus < 50
    tolerance = 1e-1;
elseif N_bus < 80
    tolerance = 1;
else
    tolerance = 1;
end
error = norm(daV);
%%
% The following is the main body of Newton-Raphson Method
% Count the iterations
counter = 0;
while (error > tolerance)
    % Stop solving, if the iterations excess certain limit
    counter = counter + 1;
    if counter > 30
        solvedBus = 0;
        disp('Iteration limit reached. Power Flow NOT solved');
        return;
    end
    % Calculate mismatches
    dPQ = double(subs(Equations,[a;V],[a_v;V_v])) - PQ;
    % Calculate Jacobian
    J_v = double(subs(J,[a;V],[a_v;V_v]));
    J_i = J_v^-1;
    % Calculate the changes
    daV = -J_i * dPQ;
    % Update the guess values
    guess = guess + daV;
    a_v([BusPV;BusPQ]) = guess(1:size([BusPV;BusPQ],1));
    V_v(BusPQ) = guess(size([BusPV;BusPQ],1)+1:end);
    % Calculate the norm of mismatch
    error = norm(daV);
end
%% Calculate active and reactive power
P_v = double(subs(EP,[a;V],[a_v;V_v]));
Q_v = double(subs(EQ,[a;V],[a_v;V_v]));

%% Return values
solvedBus = busData;
solvedBus(:,[VOLTAGE,ANGLE]) = [V_v, a_v * 180 / pi];
solvedBus(BusPV, Q_GEN) = Q_v(BusPV);
solvedBus(BusSW, [P_GEN, Q_GEN]) = [P_v(BusSW), Q_v(BusSW)];

%% Display Power Flow Results
% To be implemented
if printResult
    disp('    Bus       Voltage   Voltage   Generator Generator Load      Load');
    disp('    No.       Magnitude Angle     P         Q         P         Q');
    disp(solvedBus(:,1:7));
end
end

