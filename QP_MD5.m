function [ md5 ] = QP_MD5( busData, lineData, generator )
% This function calculates the md5 hash for the power system data
% The input data are re-shaped as a single row string

% This function is implemented using the algorithm in
%   http://en.wikipedia.org/wiki/MD5

% There is something wrong in this function.
% The result is different from standard md5
% However, this problem does not effect the results in QPOST
% since this hash function is already far more complicated than needed.

% START PREPROCESSING
% Reshape the inputs
busData = reshape(busData, 1, numel(busData));
lineData = reshape(lineData, 1, numel(lineData));
generator = reshape(generator, 1, numel(generator));

% Concatenate the matrice
Msg = [busData, lineData, generator];

% Convert to string
Msg = num2str(Msg);

% Convert to ASCII
Msg = double(Msg);

% Convert to binary
Msg = dec2bin(Msg,8);

% % Make each binary as 8 bit
% n = size(Msg,1);
% if size(Msg,2) < 8
%     o = 8 - size(Msg,2);
%     Msg = [repmat('0',n,o), Msg];
% end

% Fill binary so that its length is divisible by 512
% Add a single '1' bit
n = size(Msg,1);
Msg = [Msg; '10000000'];
% Add '0' bits
if mod(n + 1,64) < 56
    k = 56 - mod(n + 1,64);
elseif mod(n + 1,64) > 56
    k = 64 - mod(n + 1,64) + 56;
else
    k = 0;
end
Msg = [Msg; repmat('0',k,8)];
% Add "original length in bits mod (2^64)"
% MODIFIED, since 2^64 is a large number, mod is skipped for performance
l = dec2bin(n * 8);
Msg = [Msg; reshape([repmat('0',1,64 - length(l)),l]',8,8)'];

% Reshape binary rows with 512 column
Msg = reshape(Msg', 512, (n+1+k+8)/64)';
% END OF PREPROCESSING

s = zeros(1,64);
K = zeros(1,64);
modulo = 2^32;

% s specifies the per-round shift amounts
s(1:16) = [7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22];
s(17:32)= [5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20];
s(33:48)= [4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23];
s(49:64)= [6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21];

% Use binary integer part of the sines of integers (Radians) as constants
for i = 1:64
    K(i) = mod(floor(abs(sin(i)) * (2^32)), modulo);
end

% Initialize variables
a0 = hex2dec('67452301');
b0 = hex2dec('efcdab89');
c0 = hex2dec('98badcfe');
d0 = hex2dec('10325476');
N = size(Msg,1);
assumedtype = 'uint32';
for j = 1:N
    M = reshape(Msg(j,:)',32,16)';
    A = a0;
    B = b0;
    C = c0;
    D = d0;
    for i = 0:63
        if i >= 0 && i <= 15
            F = bitor(bitand(B, C, assumedtype), ...
                (bitand(bitcmp(B, assumedtype), D, assumedtype)));
            g = i;
        elseif i >= 16 && i <= 31
            F = bitor(bitand(D, B, assumedtype), ...
                (bitand(bitcmp(D, assumedtype), C, assumedtype)) ...
                ,assumedtype);
            g = mod(5 * i + 1, 16);
        elseif i >= 32 && i<= 47
            F = bitxor(bitxor(B, C, assumedtype), D, assumedtype);
            g = mod(3 * i + 5, 16);
        elseif i >= 48 && i<= 63
            F = bitxor(C, bitor(B, bitcmp(D, assumedtype), ...
                assumedtype), assumedtype);
            g = mod(7 * i, 16);
        end
        dTemp = D;
        D = C;
        C = B;
        
        Mg = bin2dec(M(g+1,:));
        R = A + F + K(i+1) + Mg;
        R = mod(R, modulo);
        R = dec2bin(R,32);
        R = R([s(i+1)+1:32, 1:s(i+1)]);
        B = B + bin2dec(R);
        B = mod(B, modulo);
        
        A = dTemp;
    end
    a0 = mod(a0 + A, modulo);
    b0 = mod(b0 + B, modulo);
    c0 = mod(c0 + C, modulo);
    d0 = mod(d0 + D, modulo);
end
md5 = [dec2hex(a0,8),dec2hex(b0,8),dec2hex(c0,8),dec2hex(d0,8)];
end

