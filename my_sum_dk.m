function S1 = my_sum_dk(S1, sig, c2c)
% takes an extra argument which specifies which dimension to filter on
% extra argument can be a vector with all dimensions that need to be
% smoothed, in which case sig can also be a vector of different smoothing
% constants

S1 = movsum(S1,2*sig(1)+1);

% [xc,yc] = meshgrid(1:16,1:16);
% c2c = sqrt((xc(:) - xc(:)').^2 + (yc(:) - yc(:)').^2);

mask = gpuArray(single(c2c < sig(2)+1));


S1 = S1*mask;

end