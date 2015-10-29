function [ varargout ] = ...
    QP_ElectricNetwork(busData, lineData, generator, scaling, faultyBus)
% This function calculate the small signal model of the electric
%   transmission network.
% This function may not provide independent states, use with caution.
% You can use the minreal function to eliminate the redundant states.
%
% Version $\Epsilon$. Qiu Qin, March 22, 2015. All Rights Reserved.
%
% The order of the states is: 
%   1. bus voltages for PQ buses, 
%   2. load currents, and 
%   3. branch currents (currents of the equivalent inductor ONLY)
%   The order follows the order appear in the data file
%
% The order of the measurements is: 
%   1. bus voltages for PQ buses, 
%   2. load currents, and 
%   3. branch currents going out from the FROM bus
%   4. branch currents going into the TO bus
%   Branch current includes 
%       the currents going through the equivalent capacitor and inductor.
%   There will be two different currents for each branch
%   All current measurements are currents going OUT of the bus.
%   The order follows the order appear in the data file
%
% CAUTION: When there is a fault...
% To create a fault at a transmission line, use the QP_ConfigSystem
%   function to create a new faultyBus before using this function.
% Look into the lineData to see the actual order when a fault exist.
% Note that the voltage of the faultyBus will NOT be a state.
%
% The input for the electric network is the generator internal voltage.
%   Generator transient impedance will be added to the line connected the
%   generator. This line usually represents a transformer.
%
% CAUTION: This program has not been verified for systems other than 9-bus
%   system, though it works on any system.
%
% If scaling is set to 1, the input and output matrices of the system is 
%   scaled, such that the unit for each input/outpu of the system is
%   "percentage of its rms value"

if nargin < 5
    faultyBus = 0;
end
if nargin == 3
    scaling = 0;
end

QP_Constants;

% Number of PQ buses
PQBus = busData(busData(:,BUS_TYPE) == PQ_BUS,BUS_NO);
N_PQ = size(PQBus,1);

% Number of Loads
LDBus = busData(busData(:,Q_LOAD) > 0, BUS_NO);
N_Load = size(LDBus,1);

% Number of Branches
N_Branch = size(lineData,1);

% Number of States
N_State = N_PQ + N_Load + N_Branch;

% Number of Inputs / Voltage Source
N_Input = sum(busData(:,BUS_TYPE) ~= PQ_BUS);

% Number of Outputs
N_Output = N_PQ + N_Load + 2 * N_Branch;

% Nornimal Frequency
f = 60;
w = 2 * pi * f;

% Convert data from PU to SI
% Line parameters
Line_C = lineData(:,LINE_B) / w;
Line_R = lineData(:,LINE_R);
Line_L = lineData(:,LINE_X) / w;
% Load parameters
Load_Z = busData(:, VOLTAGE).^2 ./ ...
    (busData(:, P_LOAD) - 1j * busData(:, Q_LOAD));
Load_R = real(Load_Z);
Load_L = imag(Load_Z) / w;

% Initialize A, B, and C Matrices
SysA = zeros(N_State);
SysB = zeros(N_State, N_Input);
SysC = zeros(N_Output,N_State);

% If a voltage is not a state, i.e. no capacitor connected to the bus, then
%   this voltage can be represented by other states. 
%   Column 1: the voltage that is not a state
%   Column 2: currents going out of the bus
%   Column 3: currents going into the bus
NonStateV = cell(10,3);
idxNSV = 0;

% DEBUG ONLY
CapArray = zeros(N_PQ,1);

% KCL for each bus
for i = 1:N_PQ
    % Lines connected to the bus
    % Lines with currents going out
    LineOut = find(lineData(:,LINE_FROM) == PQBus(i));
    % Lines with currents going in
    LineIn = find(lineData(:,LINE_TO) == PQBus(i));
    % Equivalent capacitor connected to the bus
    Cap = 0.5 * sum(Line_C([LineIn;LineOut]));
    % Row offset in A matrix
    R_Offset = 0;
    % Column offset in A matrix
    C_Offset = N_PQ + N_Load;
    % Set values in A matrix
    CapArray(i) = Cap;
    if Cap > 0
        SysA(R_Offset + i,C_Offset + LineOut) = -1 / Cap;
        SysA(R_Offset + i,C_Offset + LineIn) = 1 / Cap;
    else
        idxNSV = idxNSV + 1;
        NonStateV{idxNSV, 1} = i;
        NonStateV{idxNSV, 2} = LineOut;
        NonStateV{idxNSV, 3} = LineIn;
    end
    % Load connected to the bus
    idxLoad = find(LDBus == PQBus(i),1);
    if ~isempty(idxLoad) && Cap > 0
        SysA(R_Offset + i,N_PQ + idxLoad) = -1 / Cap;
    end
    % Measurement
    SysC(i,i) = 1;
end

% KVL for each Load
for i = 1:N_Load
    % Row offset in A matrix
    R_Offset = N_PQ;
    % Bus index
    idx = find(busData(:,BUS_NO) == LDBus(i),1);
    % State for load bus voltage
    SysA(R_Offset + i, find(PQBus == LDBus(i),1)) ...
        = 1 / Load_L(idx);
    % State for load current
    SysA(R_Offset + i, R_Offset + i) ...
        = - Load_R(idx) / Load_L(idx);
    % Measurement
    SysC(R_Offset + i, R_Offset + i) = 1;
end

% KVL for each branch
for i = 1:N_Branch
    % State Equations
    % Row offset in A matrix
    R_Offset = N_PQ + N_Load;
    % Bus index
    busFrom = lineData(i, LINE_FROM);
    busTo = lineData(i, LINE_TO);
    idxFrom = find(busData(:,BUS_NO) == busFrom,1);
    idxTo = find(busData(:,BUS_NO) == busTo,1);
    % Add generator xd if the branch connects a generator
    if busData(idxFrom, BUS_TYPE) ~= PQ_BUS
        idxGen = find(generator(:,GEN_BUS) == busFrom, 1);
        Line_L(i) = Line_L(i) + generator(idxGen, GEN_XDP) / w;
    end
    if busData(idxTo, BUS_TYPE) ~= PQ_BUS
        idxGen = find(generator(:,GEN_BUS) == busTo, 1);
        Line_L(i) = Line_L(i) + generator(idxGen, GEN_XDP) / w;
    end
    % Voltage at the FROM end
    if busData(idxFrom, BUS_TYPE) == PQ_BUS
        SysA(R_Offset + i, find(PQBus == busFrom,1)) = 1 / Line_L(i);
    else
        idxGen = find(generator(:,GEN_BUS) == busFrom, 1);
        SysB(R_Offset + i, idxGen) = 1 / Line_L(i);
    end
    % Voltage at the TO end
    if busData(idxTo, BUS_TYPE) == PQ_BUS
        SysA(R_Offset + i, find(PQBus == busTo,1)) = -1 / Line_L(i);
    else
        idxGen = find(generator(:,GEN_BUS) == busTo, 1);
        SysB(R_Offset + i, idxGen) = -1 / Line_L(i);
    end
    % Voltage of the resistance
    SysA(R_Offset + i, R_Offset + i) = -Line_R(i) / Line_L(i);
    
    % Measurement Equations
    C_Offset = N_PQ + N_Load;    
    % Current going out from the FROM bus
    R_Offset = N_PQ + N_Load;
    % Lines connected to the FROM bus
    LineOut = find(lineData(:,LINE_FROM) == busFrom);
    LineIn = find(lineData(:,LINE_TO) == busFrom);
    % Indutor Current
    SysC(R_Offset + i, C_Offset + i) = 1;
    % Total capacitance connected to the bus
    Cap = sum(Line_C([LineIn;LineOut]));
    % Load connected to the bus
    %LineOut = [LineOut;find(LDBus == busFrom,1) - N_Load]; %#ok<AGROW>
    if Cap > 0
        Cap = Line_C(i) / Cap;
        SysC(R_Offset + i, C_Offset + LineIn) = ...
            SysC(R_Offset + i, C_Offset + LineIn) + 1 * Cap;
        SysC(R_Offset + i, C_Offset + LineOut) = ...
            SysC(R_Offset + i, C_Offset + LineOut) - 1 * Cap;
    end    
    % Current going out of the TO bus
    R_Offset = N_PQ + N_Load + N_Branch;
    % Lines connected to the TO bus
    LineOut = find(lineData(:,LINE_FROM) == busTo);
    LineIn = find(lineData(:,LINE_TO) == busTo);
    % Inductor Current
    SysC(R_Offset + i, C_Offset + i) = -1;
    % Total capacitance connected to the bus
    Cap = sum(Line_C([LineIn;LineOut]));
    % Load connected to the bus
    %LineOut = [LineOut;find(LDBus == busTo,1) - N_Load]; %#ok<AGROW>
    if Cap > 0
        Cap = Line_C(i) / Cap;
        SysC(R_Offset + i, C_Offset + LineIn) = ...
            SysC(R_Offset + i, C_Offset + LineIn) + 1 * Cap;
        SysC(R_Offset + i, C_Offset + LineOut) = ...
            SysC(R_Offset + i, C_Offset + LineOut) - 1 * Cap;
    end
end

% Process the voltages that are not states
R_Offset = N_PQ + N_Load;
for i = 1:idxNSV
    % The Voltage that is not a State
    VS = NonStateV{i, 1};
    % The EQuation that can be used to represent the voltage
    VarP = sum(SysA(NonStateV{i, 3} + R_Offset,:),1);
    if isempty(VarP)
        VarP = 0;
    end
    VarN = sum(SysA(NonStateV{i, 2} + R_Offset,:),1);
    if isempty(VarN)
        VarN = 0;
    end
    EQ = VarP - VarN;
    % Use other states to represent the voltage
    if EQ(VS) == 0
        disp('WARNING: This function does not work on this system');
    elseif EQ(VS) < 0
        EQ = EQ / abs(EQ(VS));
        EQ(VS) = 0;
    else
        EQ = -EQ / EQ(VS);
        EQ(VS) = 0;
    end
    % Replace the voltage in A matrix
    SysA = SysA + SysA(:,VS) * EQ;
    SysC = SysC + SysC(:,VS) * EQ;
    SysA(:,VS) = 0 * SysA(:,VS);
end

if scaling
    % Calculate Norminal input/output magnitudes
    % Solve power flow
    busData = QP_SolvePowerFlow( busData, lineData, 0 );
    E0 = abs(QP_GenVoltage(generator,busData));
    % Input scaling
    Qu = diag(E0);
    % Bus voltage
    Qy_V = busData(busData(:,BUS_TYPE) == PQ_BUS,VOLTAGE);
    % Load current
    idx = busData(:,Q_LOAD) > 0;
    Qy_L = abs(busData(idx, VOLTAGE) ./ Load_Z(idx));
    % Branch current
    Qy_I = abs( (busData(lineData(:,LINE_FROM),VOLTAGE) ...
        .* exp(1j * busData(lineData(:,LINE_FROM),ANGLE)/180*pi) ...
        - busData(lineData(:,LINE_TO),VOLTAGE) ...
        .* exp(1j * busData(lineData(:,LINE_TO),ANGLE)/180*pi)) ...
        ./ (lineData(:,LINE_R) + 1j * lineData(:,LINE_X)) );
    % Output scaling
    Qy = diag(1./[Qy_V;Qy_L;Qy_I;Qy_I]);
    SysB = SysB * Qu;
    SysC = Qy * SysC;
end

% Process the fault
if faultyBus > 0
    faultyBus = find(PQBus == faultyBus,1);
    SysA = SysA([1:faultyBus-1, faultyBus+1:end],...
        [1:faultyBus-1, faultyBus+1:end]);
    SysB = SysB([1:faultyBus-1, faultyBus+1:end],:);
    SysC = SysC(:,[1:faultyBus-1, faultyBus+1:end]);
end

if nargout == 1
    SYS = ss(SysA,SysB,SysC,0);
    varargout{1} = SYS;
end
if nargout == 3
    varargout{1} = SysA;
    varargout{2} = SysB;
    varargout{3} = SysC;
end
end

