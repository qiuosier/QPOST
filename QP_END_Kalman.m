function [xe, ye, varargout] = QP_END_Kalman(A, B, H, P0, Q, R, ...
    x0, u0, z, u )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function implement a simple Kalman filter
% 
% The Kalman filter is implemented based on the algorithm described in
% Welch, Greg, and Gary Bishop. "An introduction to the Kalman filter." 
% (1995). Updated 2006.
% 
% The system model is described by:
% x(k) = Ax(k-1) + Bu(k-1)
% z(k) = Hx(k)
% 
% INPUT ARGUMENTS:
% A, B, and H matrices describe the LTI system
% P0 initial value of the covariance matrix
% Q: process noise covariance
% R: measurement noise covariance
% Q and R remain the same for the filtering process.
% x0 is the initial value for states, which can be an arbitray guess.
% u0 is the initial value for inputs at time t = 0.
% z: measured system outputs.
% z must be a matrix with the same ROW as H
% u: inputs to the system.
% Each COLUMN of z and u is a sample at a time.
% The number of samples in u and z should be the same.
% 
% OUTPUT ARGUMENTS:
% When using the this function, the number of outputs can be 2, 3, or 4.
% One outputs: xe, ye;
% Three outputs: xe, ye, P;
% Four outputs: xe, ye, P, L;
% where
% xe: the estimated/predicted state
% ye: the estimated/predicted measurement
% P: the covariance matrix
% L: the likelihood (that the data match the system)

N_output = size(z,1);         % Number of outputs
N_sample = size(z,2);         % Number of samples
N_state = size(A,1);          % Number of states

xe = zeros(N_state, N_sample);
ye = zeros(N_output, N_sample);
P = zeros(N_state,N_state,N_sample);
L = zeros(1,N_sample);
%% Kalman Filter
% First time update
% Update state
xe(:,1) = A * x0 + B * u0;
% Update covariance
P(:,:,1) = A * P0 * A' + Q;

for k = 1:N_sample
    %% Measurement update
    % Residual covariance
    S = H * P(:,:,k) * H' + R;
    % START DEBUG ONLY
    if rcond(S) < 1e-5
        XXX = 1; %#ok<NASGU>
    end
    % END DEBUG ONLY
    % Filter gain
    M = P(:,:,k) * H' * S^-1;
    % Predict measurement
    ye(:,k) = H * xe(:,k);
    % Measurement residual
    v = z(:,k) - ye(:,k);
    % Update state
    xe(:,k) = xe(:,k) + M * v;
    % Update covariance
    P(:,:,k) = (eye(N_state) - M * H) * P(:,:,k);
    % Likelihood
    L(k) = 1/sqrt(norm(2*pi*S)) * exp(-0.5 * v' * S^-1 * v);
    %% Time update
    if k < N_sample
        % Predict state
        xe(:,k+1) = A * xe(:,k) + B * u(:,k);
        % Predict covariance
        P(:,:,k+1) = A * P(:,:,k) * A' + Q;
    end
end
% Set output variables
if nargout == 3
    varargout{1} = P;
end
if nargout == 4
    varargout{1} = P;
    varargout{2} = L;
end

end

