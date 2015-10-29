function [ fBus, fLine, faultyPoint ] = ...
    QP_ShortToGroundFault( bus, line, faultyLine, distance )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function returns new data matrices of bus and line, 
%	with the configuration that the faultyLine is broken into two lines,
%	and faultyPoint is a new bus indicating the shorted point.
% This function adds an additional bus to the system,
%	such that the additional bus is the point shorted to ground.
% Note that the data matrices do not carry the short to ground fault.
% This function only creates new data matrices so that the fault can
%	be indicated by refering to bus.
%
% This function requires the following file:
%   QP_Constants, script defining constant variables.
%
% INPUT ARGUMENTS
% bus and line are matrices decribing the power system
% The data format is the same as PST by Prof. Chow at RPI
% Please refer to the data file for details.
% faultLine: the line that is shorted to ground
% distance: the location of the fault
%   It is the percentage of line length from the "FROM BUS" (which is
%   defined in the data file).
%
% OUTPUT ARGUMENTS
% _faultyLine_ is the line(row) number in _line_
% _distance_ is the percentage (between 0 and 1) distance of the fault 
% from the "FROM BUS" of _faultyLine_.

if distance < 0 || distance > 1
    disp('distance must be a real number in (0,1)');
end

QP_Constants;
% Number of Buses
N_bus = size(bus,1);

% Copy
fLine = line;
fBus = bus;

% Return if no fault
if faultyLine == 0
    faultyPoint = 0;
    return;
end

% Check if the fault is at the end of the bus
if distance == 0
    faultyPoint = fLine(faultyLine, LINE_FROM);
    return;
elseif distance == 1
    faultyPoint = fLine(faultyLine, LINE_TO);
    return;
end

% Add a new bus as faulty point
faultyPoint = N_bus + 1;
newBus = [faultyPoint 1.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00 3];
fBus = [fBus;newBus];
% Add the new bus between the faulty line
% Copy the faulty line data to two new sections
line1 = fLine(faultyLine, :);
line2 = fLine(faultyLine, :);
% Add new bus to two new sections
line1(LINE_TO) = faultyPoint;
line2(LINE_FROM) = faultyPoint;
% Change the line parameters (R, X, B) based on distance
line1([LINE_R, LINE_X, LINE_B]) = ...
    line1([LINE_R, LINE_X, LINE_B]) * distance;
line2([LINE_R, LINE_X, LINE_B]) = ...
    line2([LINE_R, LINE_X, LINE_B]) * (1 - distance);
% Replace the faulty line with two sections
fLine = fLine([1:faultyLine-1, faultyLine+1:end],:);
fLine = [fLine; line1; line2];

end

