function [ ddelta ] = QP_ddelta( w )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
%
% This function calculates the rotor angle deviation in 
%   classical generator model.
%
% INPUT ARGUMENT
% w: rotor speed in pu, the norminal speed is 1.
%
% OUTPUT ARGUMENT
% ddelta: rotor angle deviation in rad
%
% This function assume a norminal frequency of 60Hz


f = 60;
ddelta = 2 * pi * f * (w - 1);

end

