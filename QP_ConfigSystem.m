function [ busDataNew, lineDataNew, faultyBus ] = ...
    QP_ConfigSystem( busData, lineData, ... 
    faultyLine, faultDistance, removeLine )
% Version $\Epsilon$. Qiu Qin, March 23, 2015. All Rights Reserved.
%
% This function will configure the system data for fault or line removal
% WARNING: the data do NOT reflect the fault directly.
%   This function only return the point (as faultyBus) that should be
%   shorted to ground. Since there is no power flow solution for faulty
%   system, the bus data do NOT reflect the fault.
%
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

busDataNew = fBus;
lineDataNew = fLine;
faultyBus = faultyPoint;

end

