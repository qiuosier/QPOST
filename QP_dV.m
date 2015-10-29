function [ dVM, dVA ] = QP_dV( E, Yrr, Yrn, V0 )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
%
% This function calculates the bus voltage deviations 
%   when generator internal voltages fluctuate.
% 
% This function is used to calculate the C matrix of the small signal
%   model.
%
% INPUT ARGUMENTS:
% E: Generator internal voltage, (Complex number)
% Yrr: Admittance matrix between buses where V0 is measured
% Yrn, Admittance matrix between buese and generators
% V0: Bus voltages at operating point
% Yrr and Yrn can be obtained using:
% [ ~, ~, ~, Yrn, Yrr ] = ...
%     QP_ReducedYMatrix( busData, lineData, generator );
%
% OUTPUT ARGUMENTS:
% dVM: Voltage magnitude deviations
% dVA: Voltage angle deviations

V = - Yrr^-1 * Yrn * E;
dVM = abs(V) - abs(V0);
dVA = angle(V) -angle(V0);
end

