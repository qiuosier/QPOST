function [ pGen ] = QP_phantomGenerator( busNo, xdp, generator )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function creates a row of data for a phantom generator on a pv bus
%
% The generator will only have transient reactance
% No active or reactive power
%
% INPUT ARGUMENTS:
% busNo: the bus to which the generator connects
% xdp: the generator transient reactance
% generator: data matrix describing the existing generators in the system
% 
% OUTPUT ARGUMENT:
% pGen: a row data describing a new generator,
%   which can be appended to the existing generator data


% Base MVA = 100
baseMVA = 100;
GenNo = size(generator,1) + 1;
pGen = generator(1,:);
pGen(1:19) = ...
   [GenNo   busNo   baseMVA 0.000   0.000 ...
        0.      xdp     0       0       0 ...
        0       0       0       0       0 ...
        0       0       0       busNo];

end

