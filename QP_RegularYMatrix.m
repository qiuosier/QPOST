function [ Y ] = QP_RegularYMatrix( busData, lineData )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function calculate the admittance matrix for power flow study.
% Generator and Load impedance is not involved.
%
% Note: This function CANNOT be used when a fault exists.
%
% This function requires the following file:
%   QP_Constants, script defining constant variables
%
% INPUT ARGUMENTS:
% busData and lineData are matrices decribing the power system
% The data format is the same as PST by Prof. Chow at RPT
% Please refer to the data file for details.
% faultLocation
%
% OUTPUT ARGUMENT:
% Y: the admittance matrix



% Constants indicating columns in bus and line matrix
QP_Constants;
% Number of Buses
N_bus = size(busData,1);
% Number of Lines
N_line = size(lineData,1);

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
    % Clean local variable
end
Y = Yrr;

end

