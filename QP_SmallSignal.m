%% Small Signal Model for Power Systems with Classical Generator Model
% This function calculate the small signal model at steady state
%   using classical model for generators.
% A 3-machine 9-bus system is used as example in this document.
%
% Version $\Epsilon$. Qiu Qin, March 22, 2015. All Rights Reserved.
%
% $\Epsilon$ Revision: 
%   Generators may have different MVA Base, which will affect H, D and the
%   generator impedance. For system data with different generator MVA
%   base, the results calculated from this function is different from the
%   results of power system toolbox (PST). It is not sure which one is
%   correct. To be safe, please convert the data to 100MVA in the data file
%   before using this function.
%
%
%%
function [ SysA, SysB_Pm, SysB_V, SysB_X, SysC_V, SysC_A, busData ] = ...
    QP_SmallSignal( busData, lineData, generator, outputResult, ...
    faultyLine, faultDistance, removeLine, delta, omega)
%% Function Description
% This function calculate the small signal model at steady state
%   using classical model for generators. This function considered only
%   generators, transmission lines and static loads.
% Only the dynamics of generators are considered. Second order classical
%   model with rotor anlge and rotor speed as states is used for each 
%   generator. The rotor angle of the first generator is selected as
%   reference.
% Static loads are converted to equivalent impedance.
% Transmission network is represented by algebraic equations connecting
%   generators.
%
%% This function requires the following files:
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
%   QP_dV, calculates the bus voltage deviations 
%       when generator internal voltages fluctuate
%   QP_ddelta, calculates the rotor angle deviation 
%       in classical generator model
%   QP_dw, calculates the rotor speed deviations in classical machine model
%   QP_phantomGenerator, 
%       creates a row of data for a new generator connects to a PV bus
%
%% INPUT ARGUMENTS:
% busData, lineData and generator are matrices decribing the power system
%   The data format is the same as PST by Prof. Chow at RPI
%   Please refer to the data file for details.
% outputResult indicates whether the result will be displayed.
% * outputResult = 1 : display power flow results
% * outputResult = 0 : no output display
% The following input arguments are optional.
% faultyLine: the line that is shorted to ground
% faultDistance: the location of the fault
%   It is the portion (between 0 and 1) of line length from the "FROM BUS" 
%   (which is defined in the data file).
% removeLine: the line that is removed from the system,
%   which can be considered as mis-operated.
% Use 0 if no fault or no line removal.
% delta and omega specify the operating point of the system.
% The system will be linearized around the operating point.
% If delta or omega is not specified, 
%   the stable equilibrium point will be used
%
%% OUTPUT ARGUMENTS:
% * _SysA_, system matrix A.
% * _SysB_Pm_, input matrix B with mechanical power as input.
% * _SysB_V_, input matrix B with controlled voltage magnitude of 
%   voltage source connected to a bus as input.
% * _SysC_V_, output matrix C with bus voltage magnitude as output.
% * _SysC_A_, output matrix C with bus voltage angle as output.
% * _bus_, updated _bus_ matrix with solved power flow data.
% [EXPERIMENT] SysB_X, input matrix B with controlled reactance as input.
%
%% Model Description
% The dynamic model of each generator is given by
%
% $$\dot{\omega_i}=\frac{\omega_R}{2H_i}\left(P_{i}-E_i^2G_i-
% \sum_{j=1,j\neq i}^n E_i E_j Y_{ij}\cos{(\theta_{ij}-\delta_i+\delta_j)}
% -D_i (\omega_i-\omega_R)\right)$$
%
% $$\dot \delta_i = \Omega(\omega_i - \omega_R)$$
% 
% where
% 
% * Subscript $i$ identify different generators.
% * $\omega$ is the rotor speed in pu.
% * $\delta$ is the rotor angle in rad.
% * $\omega_R$ is the nominal rotor speed in pu.
% * $H$ is the inertia of generator.
% * $D$ is the damping coefficient of generator.
% * $P$ is the mechanical power of generator.
% * $E$ is the interal voltage of generator.
% * $Y_{ij}$ is the equivalent transfer admittance between generator $i$
%   and $j$. It is an element in the reduced admittance matrix $Y$, 
%   which only include generator nodes. All other nodes are eliminated. 
% Please see the
%   document of <QP_ReducedYMatrix.html _QP_ReducedYMatrix_> 
%   for more details.
% * $G_i$ is the real part of $Y_{ii}$
% * $\theta_{ij}$ is the angle of $Y_{ii}$
% * $\Omega$ is a constant coefficient converting speed deviation 
%   from pu to rad/s
%
% Noticed that in this function, 
%   all variables are in per unit (pu) except $\delta$,
%   which is in rad.
% In addition, it is assumed that $\omega_R = 1$ pu.
%
% More investigation is needed to estimate the numerical error.

%% Solve for Operating Point
% Power flow and equilibrium point of dynamic equations are solved.
% Please refer to the document of 
% <QP_ClassicalEquilibrium.html _QP_ClassicalEquilibrium_> for details.
% <QP_ClassicalEquilibrium.html _QP_ClassicalEquilibrium_> also calculates
% the reduced admittance matrix.
%
% _bus_ contains the bus data with solved power flow
% (voltages(pu) and angles(deg)).
%
% _delta0_ is the steady state rotor angle (rad).
%
% _P0_ is the steady state mechanical power (pu) output by the generators.
% 
% _E0_ is the steady state internal voltage of generators
% (complex numbers in pu).
%

if nargin == 4
    % No system configuration, 
    % Assume the normal operation at stable equilibrium point
    [ delta0, P0, E0, busData, Y ] = ...
        QP_ClassicalEquilibrium(busData, lineData, generator, outputResult);
    w0 = ones(size(delta0));
elseif nargin >= 7
    % Configure system, create a fault and/or remove a line
    if faultyLine == 0 || faultyLine == removeLine
        % No fault exist
        % Remove one line from the system
        lineData = QP_RemoveLine( lineData, removeLine );
        % Solve post-fault power flow
        [ delta0, P0, E0, busData, Y ] = ...
            QP_ClassicalEquilibrium(busData, lineData, ...
            generator, outputResult);
    else
        % Fault exist, use pre-fault power flow solution
        % Solve pre-fault power flow
        [ delta0, P0, E0, busData, ~ ] = ...
            QP_ClassicalEquilibrium(busData, lineData, ...
            generator, 0);
        % Remove one line from the system
        mLine = QP_RemoveLine( lineData, removeLine );
        % Adjust the index of the faulty line
        if faultyLine == removeLine
            idx = -removeLine;
        elseif faultyLine > removeLine && removeLine ~= 0
            idx = -1;
        else
            idx = 0;
        end
        % Create a short circuit fault
        % Generate a faulty point
        % The faulty point is to be used as input for other functions
        [ fBus, fLine, faultyPoint ] = ...
            QP_ShortToGroundFault(busData, mLine, ...
            faultyLine + idx, faultDistance);
        % Use fault-on Y matrix
        Y = QP_ReducedYMatrix(fBus, fLine, generator, faultyPoint);
    end
    w0 = ones(size(delta0));
    if nargin >= 8
        delta0 = delta;
    end
    if nargin == 9
        w0 = omega;
    end
end
% Voltage Magnitude of each node
EM = abs(E0);
%% Initialization
% Load constant variables. Constant variables are used to represent the
%   column header of the matrix data.
QP_Constants;

%%
% Initialize variables based on inputs.
% Number of Buses
N_bus = size(busData,1);
% Number of Generators
N_gen = size(generator,1);
% Number of generator states
N_state = 2 * N_gen;
% Number of lines
N_line = size(lineData,1);
% Inertia of Generator
H = generator(:,GEN_H) .* generator(:,3)/100;
% Damping Coefficient
D = generator(:,GEN_D) .* generator(:,3)/100;

%% Calculate Reduced Admittance matrix
% The reduced admittance matrix includes only generator nodes
%   (NOT generator bus).
% All other nodes are eliminated.
% Please refer to the document of 
% <QP_ReducedYMatrix.html _QP_ReducedYMatrix_> for details.
% Make sure that the input _bus_ includes solved power flow data.
% Y = QP_ReducedYMatrix(bus, line, generator);

%% Use Perturbation to Estimate Small Signal Model
% A perturbation of 0.01% of steady state value is applied to
%   state variables and input variables, one at a time. 
% The corresponding outputs of the dynamic (state space) equations
%   are used to calculate the linearized state space matrix.
% The perturbation rate cannot be too small as there is some numerical
%   error during the calculation. The perturbation should be significantly
%   larger than the numerical error, but small enough to make the system
%   work near linear.

% Perturbation Rate
perturb_rate = 1e-2;

%%
% Perform a perturbation at each state one at a time to estimate A matrix
SysA = zeros(N_state, N_state);
for i = 1:N_gen
    % Reset values to operating point
    delta = delta0;
    w = w0;
    % Small perturbation on $\delta$
    perturb = 1e-7;
    delta(i) = delta(i) + perturb;
    % Calculate perturbed outputs
    ddelta = QP_ddelta(w);
    dw = QP_dw(P0, EM, Y, delta, D, H, w);
    % Save results to A matrix
    for j = 1:N_gen
        SysA(2*j-1,2*i-1) = ddelta(j) / perturb;
        SysA(2*j,2*i-1) = dw(j) / perturb;
    end
    % Reset values to operating point
    delta = delta0;
    w = w0;
    % Small perturbation on $\omega$
%    perturb = perturb_rate * w(i);
    perturb = 1e-3;
    w(i) = w(i) + perturb;
    % Calculate perturbed outputs
    ddelta = QP_ddelta(w);
    dw = QP_dw(P0, EM, Y, delta, D, H, w);
    % Save results to A matrix
    for j = 1:N_gen
        SysA(2*j-1,2*i) = ddelta(j) / perturb;
        SysA(2*j,2*i) = dw(j) / perturb;
    end
end

%%
% Using mechanical power of generator ($P_i$) as input, 
%   perform a perturbation at each $P_i$ one at a time to estimate B matrix
SysB_Pm = zeros(N_state, N_gen);
for i = 1:N_gen
    % Reset values to operating point
    delta = delta0;
    w = w0;
    Pm = P0;
    % Small perturbation on Pm
    perturb = perturb_rate * Pm(i);
    Pm(i) = Pm(i) + perturb;
    ddelta = QP_ddelta(w);
    dw = QP_dw(Pm, EM, Y, delta, D, H, w);
    % Save results to B matrix
    for j = 1:N_gen
        SysB_Pm(2*j-1,i) = ddelta(j) / perturb;
        SysB_Pm(2*j,i) = dw(j) / perturb;
    end
end

%%
% Using voltage magnitude of controlled voltage source connected to a bus 
%   as input, perform a perturbation at each $P_i$ one at a time to 
%   estimate B matrix.
% Here the bus voltage is controlled indirectly by a voltage source
%   connected to the bus. If there is the generator connected to the
%   bus, the controlled voltage is the generator internal voltage,
%   otherwise, a ideal voltage source is connected to the bus with a
%   small impedance. Such ideal voltage source can be an SVC.
%   

SysB_V = zeros(N_state, N_bus);
for i = 1:N_bus
    % Reset values to operating point
    delta = delta0;
    w = w0;
    Pm = P0;
    Yp = Y;
    Ep = EM;
    % Check if there is a generator connected to the bus
    busNo = busData(i,BUS_NO);
    [r, ~] = find(generator(:, GEN_BUS) == busNo, 1);
    if isempty(r)
        % Create a phantom generator if no generator connected to the bus
        pGen = QP_phantomGenerator(busNo, 0.001, generator);
        % Change the bus to PV bus
        
        % Resolve power flow for more accurate mechanical power
        % Ideally, the PX(1:N_gen) should be the same as P0
        % However, there is some unkonwn errors
        [ ~, PX, ~, ~ ] = ...
        QP_ClassicalEquilibrium(busData, lineData, [generator;pGen], 0);
        Pm = PX(1:N_gen);
        % Re-Calculate the equivalent Y matrix
        Yp = QP_ReducedYMatrix(busData,lineData,[generator;pGen]);
        % Internal voltage for the phantom generator = bus voltage
        % Assume 0 current injection at steady state
        Ep = [EM;busData(i,VOLTAGE)];
        % Angle is also the same
        delta = [delta0;busData(i,ANGLE) / 180 * pi];
        r = N_gen + 1;
    end
    % Small perturbation on voltage
    perturb = perturb_rate * Ep(r);
    % Add perturbation to voltage
    Ep(r) = Ep(r) + perturb;
    % Perturbation Results
    ddelta = QP_ddelta(w);
    dw = QP_dw(Pm, Ep, Yp, delta, D, H, w);
    % Save results to B matrix
    for j = 1:N_gen
        SysB_V(2*j-1,i) = ddelta(j) / perturb;
        SysB_V(2*j,i) = dw(j) / perturb;
    end
end

%% 
% Using controllable line reactance
SysB_X = zeros(N_state, N_line);
for i = 1:N_line
    % Reset values to operating point
    delta = delta0;
    w = w0;
    Pm = P0;
    Ep = EM;
    xline = lineData;
    % Small perturbation on line reactance
    perturb = perturb_rate * lineData(i, LINE_X);
    % Add perturbation to line reactance
    xline(i, LINE_X) = xline(i, LINE_X) + perturb;
    % Recalculate Y matrix
    Yp = QP_ReducedYMatrix(busData,xline,generator);
    % Perturbation Results
    ddelta = QP_ddelta(w);
    dw = QP_dw(Pm, Ep, Yp, delta, D, H, w);
    % Save results to B matrix
    for j = 1:N_gen
        SysB_X(2*j-1,i) = ddelta(j) / perturb;
        SysB_X(2*j,i) = dw(j) / perturb;
    end
end

%% C Matrix, Voltage as output
V0Angle = busData(:,ANGLE) * pi / 180;
V0 = busData(:,VOLTAGE) .* (cos(V0Angle) + 1i .* sin(V0Angle));
SysC_V = zeros(N_bus,N_state);
SysC_A = zeros(N_bus,N_state);
[ ~, ~, ~, Yrn, Yrr ] = ...
    QP_ReducedYMatrix( busData, lineData, generator );
for i = 1:N_gen
    % E = E0;
    E = EM .* (cos(angle(E0)) + 1i * sin(angle(E0)));
    perturb = perturb_rate;
    E(i) = E0(i) * (cos(perturb) + 1i * sin(perturb));
    [ dVM, dVA ] = QP_dV( E, Yrr, Yrn, V0 );
    SysC_V(:,2*i-1) = dVM / perturb;
    SysC_A(:,2*i-1) = dVA / perturb;
end

%% Set Reference Rotor Angle and Transform State Space
% Using angle of generator 1 as reference,
% the following transformation is made for state variables:
% $$\delta_{21} = \delta_2 - \delta_1 $$
T = eye(N_state);
for i = 3:2:N_state
    T(i,1) = -1;
end
SysA = T * SysA * T^-1;
SysB_Pm = T * SysB_Pm;
SysB_V = T * SysB_V;
SysB_X = T * SysB_X;
SysC_V = SysC_V * T^-1;
SysC_A = SysC_A * T^-1;

%%
% Eliminate the reference state variable ($\delta_1$)
SysA = SysA(2:end,2:end);
SysB_Pm = SysB_Pm(2:end,:);
SysB_V = SysB_V(2:end,:);
SysB_X = SysB_X(2:end,:);
SysC_V = SysC_V(:,2:end);
SysC_A = SysC_A(:,2:end);

%% Display System Model
% The following function has not been implemented.
% if outputResult == 1
%     QP_DisplaySmallSignalModel( SysA, SysB_Pm, SysB_V, SysC_V, SysC_A );
% end

end

