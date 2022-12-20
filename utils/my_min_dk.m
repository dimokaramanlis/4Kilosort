function S1 = my_min_dk(S1, sig, c2c)
% takes an extra argument which specifies which dimension to filter on
% extra argument can be a vector with all dimensions that need to be
% smoothed, in which case sig can also be a vector of different smoothing
% constants

S1 = movmin(S1,2*sig(1)+1);

% [xc,yc] = meshgrid(1:16,1:16);
% c2c = sqrt((xc(:) - xc(:)').^2 + (yc(:) - yc(:)').^2);

mask = c2c < sig(2)+1;

Nmax = max(sum(mask));
maskind = repmat((1:size(mask,1))', [1 Nmax]);
for ichan = 1:size(mask,1)
    currneighbors =  find(mask(ichan,:));
    maskind(ichan,1:numel(currneighbors)) = currneighbors;
end

% S1 = S1(:, maskind);
% S1 = reshape(S1, [size(S1,1), size(mask,1), size(maskind,2)]);
% S1 = min(S1, [], 3);

Sout = S1;
for ii = 1:size(maskind,2)
    Sout = min(Sout, S1(:, maskind(:,ii)));
end
S1 = Sout;


end