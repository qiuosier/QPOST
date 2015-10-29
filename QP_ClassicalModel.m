function [ dX ] = ...
    QP_ClassicalModel( P, E, Y, D, H, X )
% Version $\Delta$. Qiu Qin, December 8, 2014. All Rights Reserved.
%
% This function calculate the rotor angle and rotor speed deviation in 
%   classical generator model.
% 
% INPUT ARGUMENTS:
% P: Mechanical power input to each generator.
% E:  Generator internal voltage magnitude.
% Y:  Reduced-network admittance matrix.
% D:  Damping coefficient of generator.
% H:  Inertia of generator.
% P, E, D, H should be column vectors.
% X = [delta; omega], a column vector
% delta: Rotor angle in Radius
% omega:  Rotor speed in pu, with a nominal speed of 1.
% 
% OUTPUT ARGUMENTS:
% dX = [d_delta; d_omega], a column vector
% d_delta: Rotor angle deviation
% d_omega: Rotor speed deviation
%
% This function is intended to replace QP_dw and QP_ddelta.



% Number of Generators
N_gen = size(P,1);
% Number of Nodes
% Here a node means a generator or an ideal voltage source.
% When the size of E is greater than P0, the additional elements in E are
%   considered as ideal voltage source in the system.
% Ideal voltage source does not have dynamic states but will be treat as 
%   a PV bus, and, a node in the reduced-network impedance matrix.
N_node = size(E, 1);
%
delta = X(1:N_gen);
omega = X(N_gen+1:end);
% Initialization
Pm = [P;zeros(N_node - N_gen, 1)];
dw = zeros(N_gen, 1);
E = abs(E);

% Based on Anderson's book, pp. 37, eq. 2.56.
%   Assumed $w_R$ to be 1 pu
for i = 1:N_gen
    dw(i) = Pm(i) - E(i)^2 * real(Y(i,i)) - D(i) * (omega(i) - 1);
    for j = 1:N_node
        if j == i 
            continue;
        end
        dw(i) = dw(i) - E(i) * E(j) * abs(Y(i,j)) * ...
            cos(angle(Y(i,j)) - delta(i) + delta(j));
    end
end
d_omega = dw ./ (2 * H);

f = 60;
d_delta = 2 * pi * f * (omega - 1);

dX= [d_delta;d_omega];
end

