%getActBands.m
%
% [binMask, Nactframes] = getBinMask(timeFreq)
%
% This function outputs a binary mask of active time-frequency units in a
% time-frequency (time x frequency) signal. It consideres all time-frequency 
% units to be active that have non-nan components.
%
% Inputs:
% timeFreq - Time-frequency signal (works also for container-frequency signals)
% 
% Outputs:
% binMask - binary mask of active time-frequency units
% Nactframes - #active frames per frequency-channel
%
% by johannes Käsbach, DTU, CAHR, johk@elektro.dtu.dk, 10. December 2015
%
%---------History-----------
% - based on getActBands.m function that estimates the active bands in a
%   time-frequency matrix
%
%---------------------------

function [binMask, Nactframes] = getBinMask(timeFreq)

% check inputs
if nargin<1||isempty(timeFreq)
    error('Please provide a time-frequency signal!')
end

% Estimate binary mask
binMask = isnan(timeFreq);

% vector with #active frames per frequency-channel
Nactframes = sum(~isnan(timeFreq),1);

end

%EOF