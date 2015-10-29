function [ rLineData ] = QP_RemoveLine( lineData, removeLine )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function removes a line from the lineData
%
% INPUT ARGUMENTS:
% lineData: Data matrix describing the line in the system
% removeLine: The line to be removed (line/row number)
%
% OUTPUT ARGUMENT:
% rLineData: Data matrix with the line removed.



rLineData = lineData([1:removeLine-1,removeLine+1:end],:);

end

