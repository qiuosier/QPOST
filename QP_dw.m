function [ dw ] = QP_dw(P, E, Y, delta, D, H, w)
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
%
% This function calculate the rotor speed deviations in 
%   classical machine model.
%
% INPUT ARGUMENTS:
% P: Mechanical power input to each generator.
% E:  Generator internal voltage magnitude.
% Y:  Reduced-network admittance matrix.
% D:  Damping coefficient of generator.
% H:  Inertia of generator.
% P, E, D, H should be column vectors.
%
% The size of P can be different from E,
% In such case, this function assumed that the additional elements in E
%   denotes ideal voltage sources.
%
% OUTPUT ARGUMENT:
% w: rotor speed deviations in pu, with nominal as 1 pu.

% Number of Generators
N_gen = size(P,1);
% Number of Nodes
N_node = size(E, 1);

% Initialization
Pm = [P;zeros(N_node - N_gen, 1)];
dw = zeros(N_gen, 1);
E = abs(E);

% Based on Anderson's book, pp. 37, eq. 2.56.
%   Assumed $w_R$ to be 1 pu
for i = 1:N_gen
    dw(i) = Pm(i) - E(i)^2 * real(Y(i,i)) - D(i) * (w(i) - 1);
    for j = 1:N_node
        if j == i 
            continue;
        end
        dw(i) = dw(i) - E(i) * E(j) * abs(Y(i,j)) * ...
            cos(angle(Y(i,j)) - delta(i) + delta(j));
    end
end
dw = dw ./ (2 * H);

end

