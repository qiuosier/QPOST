function [ v ] = QP_HankelSV( A, B, C, D )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function calculates the Hankel singular values.
% If the system is unstable, the Hankel singular values are calculated
% using the method described in
%   Balanced realization and model reduction for unstable systems
%   by Kemin Zhou, Gregory Salomon and Eva Wu
%
% INPUT ARGUMENTS:
% A, B, C and D describe the LTI system
%
% OUTPUT ARGUMENT:
% v: a vector of Hankel singular values

% Obtain similarity transformation matrix to diagonalize A
[T, lambda] = eig(A);

% Stable-unstable decomposition
AX = diag(lambda);
n = length(AX);
% k is the last negative eigenvalue after rearrangement
k = 0;
for i = 1:n
    if AX(i) > 0
        % Switch the positive eigenvalue to the back
        for j = i:n
            if AX(j) < 0
                % Switch
                AX([i,j]) = AX([j,i]);
                T(:,[i,j]) = T(:,[j,i]);
                k = i;
                break;
            end
        end
    else
        k = i;
    end
end
% Note: the T matrix obtained above is actully T^-1 in the paper.
T = T^-1;

G = ss(A, B, C, D);
GX = ss2ss(G, T);
GS = modred(GX, k+1:n, 'Truncate');
GNS = modred(GX, 1:k, 'Truncate');

% Solve Lyapunov equations
P1 = lyap(GS.a,GS.b*GS.b');
Q1 = lyap(GS.a',GS.c'*GS.c)';
P2 = lyap(-GNS.a,GNS.b*GNS.b');
Q2 = lyap(-GNS.a',GNS.c'*GNS.c)';

% Calculate Gramians
P = T^-1 * blkdiag(P1,P2) * (T^-1)';
Q = T' * blkdiag(Q1,Q2) * T;

v = real(sqrt(eig(P*Q)));

% The following code is for debug only
% disp(eig(GS.a));
% disp(' ');
% disp(eig(GNS.a));
% End debug only code
end

