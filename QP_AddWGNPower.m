function [ noisySignal ] = QP_AddWGNPower( signal, PWR )
%
% Version $\Epsilon$. Qiu Qin, Feburary 5, 2015. All Rights Reserved.
%
% Add white Gaussian noise based on Central Limit Thorem Method
% "Simulation of Communication Systems"
%   by M. C. Jeruchim, P. Balaban and K. S. Shanmugan

% Each row is considered as a channel
%
% PWR is the noise power

NX = 20;

noise = zeros(size(signal));
for i = 1:NX
    u = rand(size(signal));
    noise = noise + u;
end
% For uniform randoms in [0,1], mu = 0.5 and var = 1/12
% Adjust noise so mu = 0 and var = 1
noise = noise - NX/2;
noise = noise * sqrt(12/NX);

% Adjust noise based on SNR
N = size(signal,1);
for i = 1:N
    L = size(signal,2);
    noisePWR = sqrt(sum(noise(i,:).^2)) / L;
    scaleFactor = PWR/noisePWR;
    noise(i,:) = scaleFactor * noise(i,:);
end

noisySignal = signal + noise;

end

