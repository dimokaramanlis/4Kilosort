function dataRAW = gpufilter(buff, ops, chanMap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% filter this batch of data after common average referencing with the
% median
% buff is timepoints by channels
% chanMap are indices of the channels to be kep
% ops.fs and ops.fshigh are sampling and high-pass frequencies respectively
% if ops.fslow is present, it is used as low-pass frequency (discouraged)

% set up the parameters of the filter
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high'); % the default is to only do high-pass filtering at 150Hz
end

dataRAW = gpuArray(buff); % move int16 data to GPU
dataRAW = dataRAW';
dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
dataRAW = dataRAW(:, chanMap); % subsample only good channels

% subtract the mean from each channel
dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

% CAR, common average referencing by median
if getOr(ops, 'CAR', 1)
    dataRAW = dataRAW - median(dataRAW, 2); % subtract median across channels
end

% next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
dataRAW = filter(b1, a1, dataRAW); % causal forward filter
dataRAW = flipud(dataRAW); % reverse time
dataRAW = filter(b1, a1, dataRAW); % causal forward filter again
dataRAW = flipud(dataRAW); % reverse time back
end

