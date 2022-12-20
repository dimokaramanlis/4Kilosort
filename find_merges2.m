function st3 = find_merges2(st3, Xsim, fs,lagsim)
% this function merges clusters based on template correlation
% however, a merge is veto-ed if refractory period violations are introduced

dt = 1/1000;

Nk = size(Xsim,1);
Xsim = Xsim - diag(diag(Xsim)); % remove the diagonal of ones

st = st3;

[~, iclust] = sort(st(:,2),'ascend');
st = st(iclust, :);

% sort by firing rate first

nspk = accumarray(st(:,2), 1, [Nk, 1], @sum);
[~, isort] = sort(nspk,'descend'); % we traverse the set of neurons in descending order of firing rates

spktoind = [0;cumsum(nspk)];

fprintf('merging clusters...\n')
nmerges = 0;
ismerged = false(Nk,1);

% merge pieces with at least 5% Nspikes, and also misaligned with lag 1

msg = []; tic;
for j = 1:Nk
    
    if ismerged(j)
        continue
    end
    s1 = st(spktoind(isort(j))+1:spktoind(isort(j)+1),1)/fs;

    if numel(s1)~=nspk(isort(j))
        fprintf('lost track of spike counts') %this is a check for myself to make sure new cluster are combined correctly into bigger clusters
    end
    % sort all the pairs of this neuron, discarding any that have fewer spikes than 10 %
    [ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk'>numel(s1)*0.05) .* (nspk' < numel(s1)), 'descend');
    ienu = find(ccsort<.5, 1) - 1; % find the first pair which has too low of a correlation

	lagcurr = lagsim(isort(j),:);
	
    % for all pairs above 0.5 correlation
    for k = 1:ienu
        %s2 = st(st==ix(k), 1)/ops.fs; % find the spikes of the pair
        simlag = lagcurr(ix(k));
        %intlag = 501+(-abs(simlag):abs(simlag));
		if abs(simlag) > 3
			continue;
		end
		
        s2 = st(spktoind(ix(k))+1:spktoind(ix(k)+1),1)/fs;

        
        % compute cross-correlograms, refractoriness scores (Qi and rir), and normalization for these scores
        [K, Qi, Q00, Q01, rir] = ccg(s1, s2, 500, dt);
        Q = min(Qi/(max(Q00, Q01))); % normalize the central cross-correlogram bin by its shoulders OR by its mean firing rate
        R = min(rir); % R is the estimated probability that any of the center bins are refractory, and kicks in when there are very few spikes
        

        if K(501) == max(K)
            continue;
        end
		
        if Q<.2 && R<.05 % if both refractory criteria are met
            i = isort(j);
            % now merge j into i and move on
            smallspkind  = st3(:,2) == ix(k);
            st3(smallspkind,2) = i;
            if simlag~=0
                st3(smallspkind,1) = st3(smallspkind,1) - simlag;
            end

            st3(st3(:,2)==ix(k),2) = i; % simply overwrite all the spikes of neuron j with i (i>j by construction)
            nspk(i) = nspk(i) + nspk(ix(k)); % update number of spikes for cluster i
            %fprintf('merged %d into %d \n', isort(j), i)
            ismerged(ix(k)) = true;
            % YOU REALLY SHOULD MAKE SURE THE PC CHANNELS MATCH HERE
            nmerges = nmerges + 1;
            break; % if a pair is found, we don't need to keep going (we'll revisit this cluster when we get to the merged cluster)
        end
       
    end
    
    %--------------------------------------------------------------------------
    if mod(j, 100) == 0 || j == Nk
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('units checked: %d/%d, Nmerges: %d, time elapsed: %2.2f s...\n', ...
            j, Nk, nmerges, toc); 
        fprintf(msg);
    end
    %--------------------------------------------------------------------------
end

% %sort back
% [~, isortback] = sort(st(:,1),'ascend');
% st = st(isortback, :);

end
