function [ SysFailed ] = QP_CheckSystemFailure( delta )
% Version $\Delta$. Qiu Qin, December 8, 2014. All Rights Reserved.
%
% This function checks whether a system has failed based on the rotor angle
%   difference.
%
% This function checks the angle differences of generators
% The system is assumed to be failed if any angle difference is greater
%   than 2*pi
% This is just a rough assessment, NOT A RELIABLE METHOD.
% DO NOT use this function for precise calculation
%
% INPUT ARGUMENT:
% delta: a tall matrix, each row corresponds to a set of angle at
%   a time.
%
% OUTPUT ARGUMENT:
% SysFailed: indicate whether the system failed.
%   0 = System Not Failed, 1 = System Failed
%



Diff = delta - diag(delta(:,1)) * ones(size(delta));
if max(max(Diff)) > 2*pi
    SysFailed = 1;
else
    SysFailed = 0;
end

end

