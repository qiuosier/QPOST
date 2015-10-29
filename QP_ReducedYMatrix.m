function [ Y, varargout ] = ...
    QP_ReducedYMatrix( busData, lineData, generator,...
    faultLocation, faultDistance, removeLine )
% Version $\Epsilon$. Qiu Qin, March 22, 2014. All Rights Reserved.
%
% $\Epsilon$ Revision:
%   Note: Generator leakage impedance is not used.
%   This function will consider generator with different MVA base and
%   convert the data into 100MVA base. 
% 
% This function calculate the reduced admittance matrix, which is the
% equivalent admittance matrix between generators. All nodes without 
% generator are eliminated. 
% PQ loads are converted to equivelant impedance.
%
% This function requires the following files:
%   QP_Constants, script defining constant variables.
%   QP_RemoveLine, Removes a line from the system data
%   QP_ShortToGroundFault, 
%       Creates new data matrices by breaking the faulted line
%
% INPUT ARGUMENTS:
% busData, lineData and generator are matrices decribing the power system
% The data format is the same as PST by Prof. Chow at RPI
% Please refer to the data file for details.
% faultLocation, faultDistance and removeLine are optional
% If only faultLocation is supplied, 
%	it should refer to a bus that is shorted.
% If faultLocation, faultDistance and removeLine are all used,
% faultLocation: the line that is shorted to ground
% faultDistance: the location of the fault
%   It is the percentage of line length from the "FROM BUS" (which is
%   defined in the data file).
% removeLine: the line that is removed from the system,
%   which can be considered as mis-operated.
%
% OUTPUT ARGUMENTS:
% Y = QP_ReducedYMatrix( busData, lineData, generator,...
%   faultLocation, faultDistance, removeLine )
% [ Y, Ynn, Ynr, Yrn, Yrr ] = ...
%    QP_ReducedYMatrix( busData, lineData, generator,...
%   faultLocation, faultDistance, removeLine )
% also returns Ynn, Ynr, Yrn, Yrr.
% Y: the reduced network admittance matrix
% Y can be partitioned as [Ynn, Ynr; Yrn, Yrr];

if nargin == 3
    faultyPoint = 0;
elseif nargin == 4;
    faultyPoint = faultLocation;
else
    faultyLine = faultLocation;
        %% Remove a line
    % Remove one line from the system
    rLine = QP_RemoveLine( lineData, removeLine );
    % Adjust the index of the faulty line
    if faultyLine == removeLine
        idx = -removeLine;
    elseif faultyLine > removeLine && removeLine ~= 0
        idx = -1;
    else
        idx = 0;
    end

    %% Create a short circuit fault
    % Generate a faulty point
    % The faulty point is to be used as input for other functions
    [ fBus, fLine, faultyPoint ] = ...
        QP_ShortToGroundFault(busData, rLine, faultyLine + idx, faultDistance);
    busData = fBus;
    lineData =fLine;
end
%bus = QP_SolvePowerFlow(bus, line, 0);

% Constants indicating columns in bus and line matrix
QP_Constants;
% Number of Buses
N_bus = size(busData,1);
% Number of Lines
N_line = size(lineData,1);
% Number of Generators
N_gen = size(generator,1);
% Number of Loads
[loadBus, ~] = find(abs(busData(:,P_LOAD)+abs(busData(:,Q_LOAD))) > 1e-5);
N_load = size(loadBus, 1);

% Construct admittance matrix
Yrr = zeros(N_bus, N_bus);
for i = 1:N_line
    % Find corresponding bus row
    [rowFBus, ~] = find(busData(:,BUS_NO) == lineData(i,LINE_FROM),1);
    [rowTBus, ~] = find(busData(:,BUS_NO) == lineData(i,LINE_TO),1);
    if isempty(rowFBus) || isempty(rowTBus)
        disp('Error: Bus Not Found!');
        break;
    end
    % Add self series admittance
    yLine = 1 / (lineData(i,LINE_R) + lineData(i,LINE_X) * 1j);
    Yrr(rowFBus,rowFBus) = Yrr(rowFBus,rowFBus) + yLine;
    Yrr(rowTBus,rowTBus) = Yrr(rowTBus,rowTBus) + yLine;
    % Add self shunt admittance
    bLine = 1j *lineData(i, LINE_B) / 2;
    Yrr(rowFBus,rowFBus) = Yrr(rowFBus,rowFBus) + bLine;
    Yrr(rowTBus,rowTBus) = Yrr(rowTBus,rowTBus) + bLine;
    % Add negative mutual admittance
    Yrr(rowFBus,rowTBus) = Yrr(rowFBus,rowTBus) - yLine;
    % Symmetric
    Yrr(rowTBus,rowFBus) = Yrr(rowFBus,rowTBus);
end

% Add generator transient reactance
Ynn = zeros(N_gen,N_gen);
Ynr = zeros(N_gen,N_bus);
Yrn = zeros(N_bus,N_gen);
for i = 1:N_gen
    yGen = 1 / (generator(i, GEN_XDP) / (generator(i,3)/100) * 1j);
    Ynn(i,i) = Ynn(i,i) + yGen;
    GenBus = generator(i, GEN_BUS);
    [rowGenBus, ~] = find(busData(:,BUS_NO) == GenBus,1);
    Yrr(rowGenBus,rowGenBus) = Yrr(rowGenBus, rowGenBus) + yGen;
    Ynr(i, rowGenBus) = Ynr(i, rowGenBus) - yGen;
    Yrn(rowGenBus, i) = Ynr(i, rowGenBus);
end

% Calculate equivalent admittances of load and add them to Y matrix
YLoad = zeros(N_load, 1);
for i = 1:N_load
    LDBus = loadBus(i);
    YLoad(i) = busData(LDBus, P_LOAD) / busData(LDBus, VOLTAGE)^2 - ...
        1j * (busData(LDBus, Q_LOAD) / busData(LDBus, VOLTAGE)^2);
    Yrr(LDBus, LDBus) = Yrr(LDBus, LDBus) + YLoad(i);
end

if faultyPoint > 0
    % Set voltage of the faulty bus to 0
    % removing the corresponding row and column
    Yrr = Yrr([1:faultyPoint-1,faultyPoint+1:end],[1:faultyPoint-1,faultyPoint+1:end]);
    Ynr = Ynr(:,[1:faultyPoint-1,faultyPoint+1:end]);
    Yrn = Yrn([1:faultyPoint-1,faultyPoint+1:end],:);
end

% Calculate Reduced Y matrix based on formulat from Anderson's book pp. 41
Y = Ynn - Ynr * Yrr^-1 * Yrn;

if nargout == 5
    varargout{1} = Ynn;
    varargout{2} = Ynr;
    varargout{3} = Yrn;
    varargout{4} = Yrr;
end

end

