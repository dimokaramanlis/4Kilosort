function isduplicate = find_duplicates2(st3, Xsim, fs, estas, coords)
% this function merges clusters based on template correlation
% however, a merge is veto-ed if refractory period violations are introduced

dt = 1/1000;
similarityThres = 0.75;


Nk = size(Xsim,1);
Xsim = Xsim - diag(diag(Xsim)); % remove the diagonal of ones

% st = st3;
% [~, iclust] = sort(st(:,2),'ascend');
% st = st(iclust, :);
st = sortrows(st3, 2, 'ascend');

% sort by firing rate first
nspk = accumarray(st3(:,2), 1, [Nk, 1], @sum);
aspk = accumarray(st3(:,2), st3(:,3), [Nk, 1], @mean);

[~, isort] = sort(nspk,'descend'); % we traverse the set of neurons in ascending order of firing rates

%SORT SPIKES BY UNIT! THIS WAY ONE CAN INDEX FASTER
spktoind = [0;cumsum(nspk)];
isduplicate = false(Nk,1);
nduplicates = 0;
% require 0.8 similarity, sig. lower amplitude, bigger group contains
% almost all spikes of smaller group

msg = []; tic;
for j = 1:Nk
    if isduplicate(j)
        continue
    end
    s1 = st(spktoind(isort(j))+1:spktoind(isort(j)+1),1)/fs;

    if numel(s1)~=nspk(isort(j))
        fprintf('lost track of spike counts') %this is a check for myself to make sure new cluster are combined correctly into bigger clusters
    end
    % sort all the pairs of this neuron, discarding any that have fewer spikes than 5 %
    [ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk' < numel(s1)), 'descend'); %.* (nspk'>numel(s1)*0.05) 
    %[ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk' < numel(s1)), 'descend');
    ienu = find(ccsort<similarityThres, 1) - 1; % find the first pair which has too low of a correlation

	%lagcurr = lagsim(isort(j),:);
    % for all pairs above 0.5 correlation
    for k = 1:ienu
        s2 = st(spktoind(ix(k))+1:spktoind(ix(k)+1),1)/fs;
        [K, ~, ~, ~, ~] = ccg(s1, s2, 2*size(estas,1), dt);
        if max(K)/numel(s2) > 0.5 && aspk(ix(k))/aspk(isort(j)) < 0.8
            isduplicate(ix(k)) = true;
            nduplicates = nduplicates + 1;
%             subplot(1,3,1);cla;
%             plotSpikeTemplate(estas(:,:,isort(j)), coords,'k')
%             plotSpikeTemplate(estas(:,:,ix(k)), coords,'r')
%             subplot(1,3,2)
%             plot(K)
%             subplot(1,3,3)
%             plot(s1,st(spktoind(isort(j))+1:spktoind(isort(j)+1),3),'k.',...
%                 s2,st(spktoind(ix(k))+1:spktoind(ix(k)+1),3),'r.')
%             pause;
        end
    end


    %--------------------------------------------------------------------------
    if mod(j, 100) == 0 || j == Nk
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('units checked: %d/%d, Nduplicates: %d, time elapsed: %2.2f s...\n', ...
            j, Nk, nduplicates, toc); 
        fprintf(msg);
    end
    %--------------------------------------------------------------------------
end

% [irow, icol] = find(tril(Xsim) > similarityThres);
% 
% msg = []; tic;
% 
% for ipair = 1:numel(icol)
%     
%     s1 = st((spktoind(irow(ipair))+1):spktoind(irow(ipair)+1),1)/fs;
%     s2 = st(spktoind(icol(ipair))+1:spktoind(icol(ipair)+1),1)/fs;
%     [K, ~, ~, ~, ~] = ccg(s1, s2, 50, dt);
%     
%     minspk = min(nspk(icol(ipair)), nspk(irow(ipair)));
%     %--------------------------------------------------------------------------
%     % label adjacent pair
%     if  sum(K(50:52))/minspk > 0.15 %0.2
%         adjmat(irow(ipair), icol(ipair)) = true;
%     end
%     %--------------------------------------------------------------------------
%     if mod(ipair, 500) == 0 || ipair == numel(icol)
%         fprintf(repmat('\b', 1, numel(msg)));
%         msg = sprintf('pairs checked: %d/%d, Nhits: %d, time elapsed: %2.2f s...\n', ...
%             ipair, numel(icol), nnz(adjmat), toc); 
%         fprintf(msg);
%     end
%     %--------------------------------------------------------------------------
% end
% 
% adjmat = adjmat | adjmat'; % expand adjacency matrix


end
