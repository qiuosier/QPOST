function [ ResampledData ] =...
    QP_ResampleData( TimeSeriesData, SamplePerCycle )
% Version $\Delta$. Qiu Qin, December 9, 2014. All Rights Reserved.
% 
% This function resamples a timeseries data and 
%	provides sample data at a fixed interval.
%
% INPUT ARGUMENTS:
% TimeSeriesData: The data to be resampled
% SamplePerCycle: Number of samples per cycle in the resampled data
% 	Frequency of 60 Hz is used in this function, i.e. 60 cycles/second
%
% OUTPUT ARGUMENT:
% ResampledData: The resampled timeseries

%% Resample Data
% Data Length (seconds)
EndTime = TimeSeriesData.Time(end);
% 
ResampleTime = 0:1/60/SamplePerCycle:EndTime;
%
ResampledData = resample(TimeSeriesData, ResampleTime);

end

