function [ varargout ] = ...
    QP_9B_SSModel( FaultyLine, FaultDistance, RemoveLine )
% Version $\Delta$. Qiu Qin, December 8, 2014. All Rights Reserved.
%
% This function calculates the state space model for the transmission
%   network for the 9-bus system. 
% The transmission network is modeled as an LTI system.
% Generator internal voltages are considered as inputs.
% Voltages and currents at bus 4, 7, and 9 are considered as outputs.
%
% Once there is a fault, the additional state will always be the current
%   going from the bus with a load to the faulted point.
% The original state of the faulted line will be replaced by the current
%   going from the bus with a generator to the faulted point.
% The additional state will be added as the last state of the system.
%
% Once a line is removed, the corresponding state is also removed.
%
% Please refer to QP_9B_StateSpaceEquation.pdf for the order of states.
% 
% This function requires the following files:
%   QP_9B_data3m9b, the data file.
%   QP_9B_LineDefinition, script defining constant variables for
%       transmission lines in the 9-bus system.
%   QP_Constants, script defining constant variables.
%
% INPUT ARGUMENTS:
% FaultyLine: the line number for the faulted line
% FaultDistance: the distance between the fault and the bus
%   with a smaller bus number. e.g. 0.3 for Line 57 indicates
%   a fault is 30% of the line length from bus 5.
% RemoveLine: the line number for the removed line
% Use 0 if no fault or no line removal.
% 
% OUTPUT ARGUMENTS:
% This function can have either 1 or 3 outputs.
%
% ONE OUTPUT
% The output will be a linear system object: SYS
%
% THREE OUTPUTS
% The outputs are SysA, SysB, and SysC,
%   which corresponds to the A, B and C matrices.
%
% EXAMPLE:
% Assuming QP_9B_LineDefinition is loaded, LINEXX can be used to refer a
%  transmission line in the system.
%
% SYS = QP_9B_SSModel( 0, 0, 0 )
%
% Calculates the LTI system for normal operating condition
% 
% [A, B, C ] = QP_9B_SSModel( LINE57, 0.5, LINE78 )
%
% Calculates the A, B and C matrix when LINE57 is shorted at midpoint and
%   LINE78 is removed.
%


% Ignore FaultyLine if it is removed.
if RemoveLine == FaultyLine
    FaultyLine = 0;
end

%% Load Data File
% Using the same data file format as Power System Toolbox by Chow
run('QP_9B_data3m9b.m');
generator = mac_con;
QP_Constants;
QP_9B_LineDefinition;

%% System Parameters
f = 60;
w = 2 * pi * f;
% Line Parameters, R, L, C
SS_C = line(:,LINE_B) / w; %#ok<*NODEF>
SS_R = line(:,LINE_R);
SS_L = line(:,LINE_X) / w;
% Load impedance
LR = [0.6836;1.0255;0.9194];
LL = [0.2735;0.3418;0.3218] / w;
% Generator and Transformer impedance
LG = (generator(:,GEN_XDP) + line(1:3,LINE_X)) / w;

% Adjustment when there is a line removal
if RemoveLine > 0
    SS_C(RemoveLine) = 0;
end
% Shunt admittance at each bus
CC = zeros(1,9);
CC(4) = SS_C(LINE45)/2 + SS_C(LINE46)/2;
CC(5) = SS_C(LINE45)/2 + SS_C(LINE57)/2;
CC(6) = SS_C(LINE69)/2 + SS_C(LINE46)/2;
CC(7) = SS_C(LINE57)/2 + SS_C(LINE78)/2;
CC(8) = SS_C(LINE78)/2 + SS_C(LINE89)/2;
CC(9) = SS_C(LINE69)/2 + SS_C(LINE89)/2;
% Adjustment when there is a line fault
switch FaultyLine
    case LINE45
        CC(4) = FaultDistance * SS_C(LINE45)/2 + SS_C(LINE46)/2;
        CC(5) = (1-FaultDistance) * SS_C(LINE45)/2 + SS_C(LINE57)/2;
    case LINE57
        CC(5) = SS_C(LINE45)/2 + FaultDistance * SS_C(LINE57)/2;
        CC(7) = (1-FaultDistance) * SS_C(LINE57)/2 + SS_C(LINE78)/2;
    case LINE46
        CC(4) = SS_C(LINE45)/2 + FaultDistance * SS_C(LINE46)/2;
        CC(6) = SS_C(LINE69)/2 + (1-FaultDistance) * SS_C(LINE46)/2;
    case LINE69
        CC(6) = FaultDistance * SS_C(LINE69)/2 + SS_C(LINE46)/2;
        CC(9) = (1-FaultDistance) * SS_C(LINE69)/2 + SS_C(LINE89)/2;
    case LINE78
        CC(7) = FaultDistance * SS_C(LINE57)/2 + SS_C(LINE78)/2;
        CC(8) = (1-FaultDistance) * SS_C(LINE78)/2 + SS_C(LINE89)/2;
    case LINE89
        CC(8) = FaultDistance * SS_C(LINE78)/2 + SS_C(LINE89)/2;
        CC(9) = (1-FaultDistance) * SS_C(LINE69)/2 + SS_C(LINE89)/2;
end

%% State Definition
V4 = 1;
V5 = 2;
V6 = 3;
V7 = 4;
V8 = 5;
V9 = 6;
I5A = 7;
I6B = 8;
I8C = 9;
I41 = 10;
I45 = 11;
I46 = 12;
I72 = 13;
I75 = 14;
I78 = 15;
I93 = 16;
I96 = 17;
I98 = 18;
ILF = 19;
%% A Matrix
if FaultyLine > 0
    N_state = 19;
else
    N_state = 18;
end
SysA = zeros(N_state);
SysA(V4,[I41,I45,I46]) = [-1, -1, -1] / CC(4);
SysA(V5,[I5A,I75,I45]) = [-1, 1, 1] / CC(5);
SysA(V6,[I6B,I96,I46]) = [-1, 1, 1] / CC(6);
SysA(V7,[I72,I78,I75]) = [-1, -1, -1] / CC(7);
SysA(V8,[I8C,I78,I98]) = [-1, 1, 1] / CC(8);
SysA(V9,[I93,I98,I96]) = [-1, -1, -1] / CC(9);
SysA(I5A,[I5A,V5]) = [-LR(1), 1] / LL(1);
SysA(I6B,[I6B,V6]) = [-LR(2), 1] / LL(2);
SysA(I8C,[I8C,V8]) = [-LR(3), 1] / LL(3);
SysA(I41,V4) = 1 / LG(1);
SysA(I45,[I45,V4,V5]) = [-SS_R(LINE45), 1, -1] / SS_L(LINE45);
SysA(I46,[I46,V4,V6]) = [-SS_R(LINE46), 1, -1] / SS_L(LINE46);
SysA(I72,V7) = 1 / LG(2);
SysA(I75,[I75,V7,V5]) = [-SS_R(LINE57), 1, -1] / SS_L(LINE57);
SysA(I78,[I78,V7,V8]) = [-SS_R(LINE78), 1, -1] / SS_L(LINE78);
SysA(I93,V9) = 1 / LG(3);
SysA(I96,[I96,V9,V6]) = [-SS_R(LINE69), 1, -1] / SS_L(LINE69);
SysA(I98,[I98,V9,V8]) = [-SS_R(LINE89), 1, -1] / SS_L(LINE89);
% Adjustment when there is a line fault
% IGF: Current from Generator bus to Fault
% ILF: Current from Load bus to Fault
% DGF: Distance from Generator bus to Fault
% DLF: Distance from Load bus to Fault
switch FaultyLine
    case LINE45
        IGF = I45;
        DGF = FaultDistance;
        DLF = 1 - FaultDistance;
        SysA([V4,V5,IGF,ILF],:) = 0*SysA([V4,V5,IGF,ILF],:);
        SysA(V4,[I41,I46,IGF]) = [-1, -1, -1] / CC(4);
        SysA(V5,[I5A,I75,ILF]) = [-1, 1, -1] / CC(5);
        SysA(IGF,[V4,IGF]) = [1, -SS_R(LINE45)*DGF] / (SS_L(LINE45)*DGF);
        SysA(ILF,[V5,ILF]) = [1, -SS_R(LINE45)*DLF] / (SS_L(LINE45)*DLF);
    case LINE57
        IGF = I75;
        DGF = 1 - FaultDistance;
        DLF = FaultDistance;
        SysA([V7,V5,IGF,ILF],:) = 0*SysA([V7,V5,IGF,ILF],:);
        SysA(V7,[I72,I78,IGF]) = [-1, -1, -1] / CC(7);
        SysA(V5,[I5A,I45,ILF]) = [-1, 1, -1] / CC(5);
        SysA(IGF,[V7,IGF]) = [1, -SS_R(LINE57)*DGF] / (SS_L(LINE57)*DGF);
        SysA(ILF,[V5,ILF]) = [1, -SS_R(LINE57)*DLF] / (SS_L(LINE57)*DLF);
    case LINE46
        IGF = I46;
        DGF = FaultDistance;
        DLF = 1 - FaultDistance;
        SysA([V4,V6,IGF,ILF],:) = 0*SysA([V4,V6,IGF,ILF],:);
        SysA(V4,[I41,I45,IGF]) = [-1, -1, -1] / CC(4);
        SysA(V6,[I6B,I96,ILF]) = [-1, 1, -1] / CC(6);
        SysA(IGF,[V4,IGF]) = [1, -SS_R(LINE46)*DGF] / (SS_L(LINE46)*DGF);
        SysA(ILF,[V6,ILF]) = [1, -SS_R(LINE46)*DLF] / (SS_L(LINE46)*DLF);
    case LINE69
        IGF = I96;
        DGF = 1 - FaultDistance;
        DLF = FaultDistance;
        SysA([V9,V6,IGF,ILF],:) = 0*SysA([V9,V6,IGF,ILF],:);
        SysA(V9,[I93,I98,IGF]) = [-1, -1, -1] / CC(9);
        SysA(V6,[I6B,I46,ILF]) = [-1, 1, -1] / CC(6);
        SysA(IGF,[V9,IGF]) = [1, -SS_R(LINE69)*DGF] / (SS_L(LINE69)*DGF);
        SysA(ILF,[V6,ILF]) = [1, -SS_R(LINE69)*DLF] / (SS_L(LINE69)*DLF);
    case LINE78
        IGF = I78;
        DGF = FaultDistance;
        DLF = 1 - FaultDistance;
        SysA([V7,V8,IGF,ILF],:) = 0*SysA([V7,V8,IGF,ILF],:);
        SysA(V7,[I72,I75,IGF]) = [-1, -1, -1] / CC(7);
        SysA(V8,[I8C,I98,ILF]) = [-1, 1, -1] / CC(8);
        SysA(IGF,[V7,IGF]) = [1, -SS_R(LINE78)*DGF] / (SS_L(LINE78)*DGF);
        SysA(ILF,[V8,ILF]) = [1, -SS_R(LINE78)*DLF] / (SS_L(LINE78)*DLF);
    case LINE89
        IGF = I98;
        DGF = 1 - FaultDistance;
        DLF = FaultDistance;
        SysA([V9,V8,IGF,ILF],:) = 0*SysA([V9,V8,IGF,ILF],:);
        SysA(V9,[I93,I96,IGF]) = [-1, -1, -1] / CC(9);
        SysA(V8,[I8C,I78,ILF]) = [-1, 1, -1] / CC(8);
        SysA(IGF,[V9,IGF]) = [1, -SS_R(LINE89)*DGF] / (SS_L(LINE89)*DGF);
        SysA(ILF,[V8,ILF]) = [1, -SS_R(LINE89)*DLF] / (SS_L(LINE89)*DLF);
end
% Adjustment when there is a line removal
switch RemoveLine
    case LINE45
        SysA(V4,I45) = 0;
        SysA(V5,I45) = 0;
    case LINE57
        SysA(V5,I75) = 0;
        SysA(V7,I75) = 0;
    case LINE46
        SysA(V4,I46) = 0;
        SysA(V6,I46) = 0;
    case LINE69
        SysA(V6,I96) = 0;
        SysA(V9,I96) = 0;
    case LINE78
        SysA(V7,I78) = 0;
        SysA(V8,I78) = 0;
    case LINE89
        SysA(V8,I98) = 0;
        SysA(V9,I98) = 0;
end

%% B Matrix, same for all cases
if FaultyLine == 0
    SysB = zeros(18,3);
else
    SysB = zeros(19,3);
end
SysB(10,1) = -1/LG(1);
SysB(13,2) = -1/LG(2);
SysB(16,3) = -1/LG(3);
%% C Matrix
% Output Definition
Y_V4 = 1;
Y_I41 = 2;
Y_I42 = 3;
Y_I43 = 4;
Y_V7 = 5;
Y_I71 = 6;
Y_I72 = 7;
Y_I73 = 8;
Y_V9 = 9;
Y_I91 = 10;
Y_I92 = 11;
Y_I93 = 12;
% Normal
SysC = zeros(12,N_state);
SysC(Y_V4,V4) = 1;
SysC(Y_I41,I41) = 1;
SysC(Y_I42,[I41,I45,I46]) = [-SS_C(LINE45)/(2*CC(4)), ...
    1-SS_C(LINE45)/(2*CC(4)),-SS_C(LINE45)/(2*CC(4))];
SysC(Y_I43,[I41,I45,I46]) = [-SS_C(LINE46)/(2*CC(4)), ...
    -SS_C(LINE46)/(2*CC(4)),1-SS_C(LINE46)/(2*CC(4))];
SysC(Y_V7,V7) = 1;
SysC(Y_I72,I72) = 1;
SysC(Y_I71,[I72,I75,I78]) = [-SS_C(LINE75)/(2*CC(7)), ...
    1-SS_C(LINE75)/(2*CC(7)),-SS_C(LINE75)/(2*CC(7))];
SysC(Y_I73,[I72,I75,I78]) = [-SS_C(LINE78)/(2*CC(7)), ...
    -SS_C(LINE78)/(2*CC(7)),1-SS_C(LINE78)/(2*CC(7))];
SysC(Y_V9,V9) = 1;
SysC(Y_I93,I93) = 1;
SysC(Y_I91,[I93,I96,I98]) = [-SS_C(LINE69)/(2*CC(9)), ...
    1-SS_C(LINE69)/(2*CC(9)),-SS_C(LINE69)/(2*CC(9))];
SysC(Y_I92,[I93,I96,I98]) = [-SS_C(LINE89)/(2*CC(9)), ...
    -SS_C(LINE89)/(2*CC(9)),1-SS_C(LINE89)/(2*CC(9))];
% Adjustment when there is a line fault
switch FaultyLine
    case LINE45
        CX = -SS_C(LINE45) * FaultDistance /(2*CC(4));
        SysC(Y_I42,:) = 0*SysC(Y_I42,:);
        SysC(Y_I42,[I41,I45,I46]) = [CX, 1 + CX, CX];
    case LINE57
        CX = -SS_C(LINE57) * (1 - FaultDistance) / (2*CC(7));
        SysC(Y_I71,:) = 0*SysC(Y_I42,:);
        SysC(Y_I71,[I72,I75,I78]) = [CX, 1 + CX, CX];
    case LINE46
        CX = -SS_C(LINE46) * FaultDistance /(2*CC(4));
        SysC(Y_I43,:) = 0*SysC(Y_I42,:);
        SysC(Y_I43,[I41,I45,I46]) = [CX, CX, 1 + CX];
    case LINE69
        CX = -SS_C(LINE69) * (1 - FaultDistance) / (2*CC(9));
        SysC(Y_I91,:) = 0*SysC(Y_I42,:);
        SysC(Y_I91,[I93,I96,I98]) = [CX, 1 + CX, CX];
    case LINE78
        CX = -SS_C(LINE78) * FaultDistance /(2*CC(7));
        SysC(Y_I73,:) = 0*SysC(Y_I42,:);
        SysC(Y_I73,[I72,I75,I78]) = [CX, CX, 1 + CX];
    case LINE89
        CX = -SS_C(LINE89) * (1 - FaultDistance)/(2*CC(9));
        SysC(Y_I92,:) = 0*SysC(Y_I42,:);
        SysC(Y_I92,[I93,I96,I98]) = [CX, CX, 1 + CX];
end
% Adjustment when there is a line removal
switch RemoveLine
    case LINE45
        SysC(Y_I42,:) = 0 * SysC(Y_I42,:);
        SysC(Y_I43,I45) = 0;
    case LINE57
        SysC(Y_I71,:) = 0 * SysC(Y_I42,:);
        SysC(Y_I73,I75) = 0;
    case LINE46
        SysC(Y_I43,:) = 0 * SysC(Y_I42,:);
        SysC(Y_I42,I46) = 0;
    case LINE69
        SysC(Y_I91,:) = 0 * SysC(Y_I42,:);
        SysC(Y_I92,I96) = 0;
    case LINE78
        SysC(Y_I73,:) = 0 * SysC(Y_I42,:);
        SysC(Y_I71,I78) = 0;
    case LINE89
        SysC(Y_I92,:) = 0 * SysC(Y_I42,:);
        SysC(Y_I91,I98) = 0;
end

% Remove a state if there is a line removal

switch RemoveLine
    case LINE45
        RM_State = I45;
    case LINE57
        RM_State = I75;
    case LINE46
        RM_State = I46;
    case LINE69
        RM_State = I96;
    case LINE78
        RM_State = I78;
    case LINE89
        RM_State = I98;
end
if RemoveLine > 0
    SysA = SysA([1:RM_State-1,RM_State+1:end],...
        [1:RM_State-1,RM_State+1:end]);
    SysB = SysB([1:RM_State-1,RM_State+1:end],:);
    SysC = SysC(:,[1:RM_State-1,RM_State+1:end]);
end

SYS = ss(SysA,SysB,SysC,0);


if nargout == 1
    varargout{1} = SYS;
end
if nargout == 3
    varargout{1} = SysA;
    varargout{2} = SysB;
    varargout{3} = SysC;
end
end

