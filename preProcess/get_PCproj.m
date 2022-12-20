function Us = get_PCproj(S1, row, col, wPCA, maskMaxChans,nt0min, c2c)

[~, nChan] = size(S1);

dt =-nt0min + (1:size(wPCA,1));

inds = repmat(row', numel(dt), 1) + repmat(dt', 1, numel(row));

clips = reshape(S1(inds, :), numel(dt), numel(row), nChan);


% mask = repmat([1:nChan], [numel(row) 1]) - repmat(col, 1, nChan);
% Mask(1,:,:) = abs(mask)<maskMaxChans;

cmat = c2c < maskMaxChans;
Mask(1,:,:) = cmat(col,:);

% pcaSpikes=zeros(numel(row),numel(dt));
% for ii=1:numel(row)
%    pcaSpikes(ii,:)=gather(S1(inds(:,ii),col(ii))); 
% end

clips = bsxfun(@times, clips , Mask);

Us = wPCA' * reshape(clips, numel(dt), []);
Us = reshape(Us, size(wPCA,2), numel(row), nChan);

Us = permute(Us, [2 3 1]);