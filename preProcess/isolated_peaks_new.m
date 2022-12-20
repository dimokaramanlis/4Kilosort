function [row, col, mu] = isolated_peaks_new(S1, ops, c2c)

loc_range = getOr(ops, 'loc_range', [5 4]);
long_range = getOr(ops, 'long_range', [30 6]);
Th = ops.spkTh;
nt0 = ops.nt0;

smin = my_min_dk(S1, loc_range, c2c);
%smin = my_min(S1, loc_range, [1 2]);
peaks = single(S1<smin+1e-3 & S1<Th);

sum_peaks = my_sum_dk(peaks, long_range, c2c);
%sum_peaks = my_sum(peaks, long_range, [1 2]);
peaks = peaks .* (sum_peaks<1.2) .* S1;

% exclude temporal buffers
peaks([1:nt0 end-nt0:end], :) = 0;

% exclude edge channels 
% noff = 8;
% peaks(:, [1:noff end-noff+ [1:noff]]) = 0;

[row, col, mu] = find(peaks);

mu = - mu;
