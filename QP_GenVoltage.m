function [ E ] = QP_GenVoltage( generator, busData )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function calculates generator internal voltage
% 
% Detail descriptions of how to calculate the generator internal voltage
% can be found at:
% Power System Control and Stability, P. M. Anderson and A. A. Fouad
%   pp. 39, eq. 2.61
%
% This function requires QP_Constants.m
%
% INPUT ARGUMENTS:
% busData and generator are matrices decribing the power system
% The data format is the same as PST by Prof. Chow at RPI
% Please refer to the data file for details.
%
% OUTPUT ARGUMENT:
% E: Internal voltages of the generators in the system.
% E is a complex column vector.



% Load Constants
QP_Constants;
% Number of Generators
N_gen = size(generator,1);
% Calculate generator internal voltage
E = zeros(N_gen, 1);
for i = 1:N_gen
    GenBus = generator(i, GEN_BUS);
    [rowGenBus, ~] = find(busData(:,BUS_NO) == GenBus,1);
    E(i) = busData(rowGenBus, VOLTAGE) + ...
        busData(rowGenBus, Q_GEN) * ... 
        generator(i, GEN_XDP) / busData(rowGenBus, VOLTAGE) + ...
        1j * busData(rowGenBus, P_GEN) * ... 
        generator(i, GEN_XDP) / busData(rowGenBus, VOLTAGE);
end

% Calculate the angle of generator internal voltage
%   with is also the rotor angle in classical model
GenBus = generator(:, GEN_BUS);
rowGenBus = zeros(size(GenBus));
for i = 1:length(GenBus)
    [rowGenBus(i), ~] = find(busData(:,BUS_NO) == GenBus(i),1);
end
delta = angle(E) + busData(rowGenBus, ANGLE) / 180 * pi;

E = abs(E) .* (cos(delta) + 1j * sin(delta));

end

