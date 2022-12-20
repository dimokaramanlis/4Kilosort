function rez = mergeChannelChunks(allrez,chunckInds)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% define all pairs to be checked, determine overlap region, find units
% that are in that region, compare waveforms and number of spikes

Nchunks = numel(allrez);
NfiltChunk   = cellfun(@(x) size(x.simScore,1),allrez);
cumtemps = [0;cumsum(NfiltChunk)];

NspikesChunk = cellfun(@(x) size(x.st3,1),allrez);
NspikesTOT   = sum(NspikesChunk);
cumspikes = [0;cumsum(NspikesChunk)];
%============================================================
% combine spikes

st3 = zeros(NspikesTOT, 3);
for ichunk = 1:Nchunks
    istart = cumspikes(ichunk)+1;
    iend   = cumspikes(ichunk+1);
    currspikes = allrez{ichunk}.st3(:,1:3);
    currspikes(:,2) = currspikes(:,2) + cumtemps(ichunk);
    st3(istart:iend,:) = currspikes;
end
% sort spikes
st3 = sortrows(st3, 1);
% remove duplicate spikes
st3 = unique(st3, 'rows');
% re-index templates
[~, ~, ic] = unique(st3(:,2));
st3(:,2)   = ic;
Nfilt = max(ic);
% remove templates with very few spikes (e.g. < 20)
Nspikes = accumarray(st3(:,2),1,[Nfilt 1], @sum);
spkthres1 = 100;
irem    = find(Nspikes < spkthres1);
st3(ismembc(st3(:,2), irem), :) = [];
% re-index templates
[~, ~, ic] = unique(st3(:,2));
st3(:,2)   = ic;
Nfilt = max(ic);
fprintf('Removed %d units with less than %d spikes, %d units remaining\n',  numel(irem), spkthres1,Nfilt);
%============================================================
ops = allrez{1}.ops;
Nchan = ops.Nchan;
estas   = zeros(ops.nt0, Nchan, Nfilt,'single');
normfac = zeros(Nfilt,1,'single');

dt =-ops.nt0min + (1:ops.nt0);
Nbatch      = ops.Nbatch;
NT  	= ops.NT;
batchstart = 0:NT:NT*Nbatch;
%============================================================
fid = fopen(ops.fproc, 'r');

fprintf('Merging sorted results... '); msg = [];

for ibatch = 1:4:Nbatch


    offset = 2 * Nchan*batchstart(ibatch);
    fseek(fid, offset, 'bof');
    dat = fread(fid, [NT Nchan], '*int16');
    % move data to GPU and scale it
    dat = double(gpuArray(dat));
    dat = dat / ops.scaleproc;
    
    if ibatch == 1
        ioffset = 0;
    else
        ioffset= ops.ntbuff;
    end
    stmin = (NT-ops.ntbuff)*(ibatch-1) - ioffset;
    stmax = stmin + NT;

    currid = st3(:,1)>stmin & st3(:,1)<stmax;
    stcurr = st3(currid,:);
    
    if isempty(stcurr), continue, end

    stimescurr = stcurr(:,1) - stmin;
    irem = stimescurr+max(dt) > NT | stimescurr + min(dt) <1;
    stcurr(irem, :) = [];
    stimescurr(irem) = [];

%     tic;
%     spkshapes = reshape(dat(stimescurr + dt,:), [size(stimescurr,1), numel(dt)*Nchan]);
%     
%     nspkcurr = accumarray(stcurr(:,2),1,[Nfilt,1],@sum);
%     [iun,~,ic_curr] = unique(stcurr(:,2));
%     currsparse = sparse(1:numel(stimescurr), ic_curr, 1,numel(stimescurr), max(ic_curr));
%     currsum = spkshapes'*currsparse;
%     currsum = reshape(currsum, [numel(dt),Nchan,max(ic_curr)]);
%     estas(:,:,iun) = estas(:,:,iun) + currsum;
%     toc;
    
    [iun,~,ic_curr] = unique(stcurr(:,2));
    spgpu = sparse(1:numel(stimescurr), ic_curr, gpuArray(1),numel(stimescurr), max(ic_curr));
    gpucurr = gpuArray.zeros([numel(dt), Nchan, numel(iun)],'single');
    for ii = 1:numel(dt)
        gpucurr(ii, :, :) = dat(stimescurr + dt(ii),:)'*spgpu;
    end
    estas(:,:,iun) = estas(:,:,iun) + gather(gpucurr);

    nspkcurr = accumarray(stcurr(:,2),1,[Nfilt,1],@sum);
    normfac = normfac + nspkcurr;

    % print message
    if ops.verbose
        fprintf(repmat('\b', 1, numel(msg)));
        msg             = sprintf('Time %2.0f min, batch %d/%d \n', ...
            toc/60, ibatch,Nbatch);        
        fprintf(msg);
    end

end
fclose(fid);

%normalize eSTAs, make sure no NaNs are produced!!!
estas = estas./reshape(normfac, [1,1,Nfilt]);

% remove ugly templates or templates with no spikes
%qual = sq(median(range(estas./max(abs(estas),[],[1 2]),1),2));
qual = sq(median(std(estas./max(abs(estas),[],[1 2]),1),2));
%irem = qual > 0.1;
irem    = find(qual > 0.1 | normfac == 0);
st3(ismembc(st3(:,2), irem), :) = [];
% re-index templates
[ia, ~, ic] = unique(st3(:,2));
estas = estas(:,:,ia);
normfac = normfac(ia);
st3(:,2)   = ic;
Nfilt = max(ic);
fprintf('Removed %d units with noise templates, %d units remaining\n', numel(irem), Nfilt);

% calculate similarities
[Xsim,lagsim] = compareSpikeTemplatesCMOS(estas);
Xsim(isnan(Xsim)) = 0;
Xsim = Xsim - diag(diag(Xsim));

% remove more units
Nspikes = accumarray(st3(:,2),1,[Nfilt 1], @sum);
irem    = find(max(Xsim,[],2) < 0.5 & Nspikes < 1000);
st3(ismembc(st3(:,2), irem), :) = [];
% re-index templates
[ia, ~, ic] = unique(st3(:,2));
estas = estas(:,:,ia);
normfac = normfac(ia);
st3(:,2)   = ic;
Xsim   = Xsim(ia,ia);
lagsim = lagsim(ia,ia);

Nfilt = max(ic);


fprintf('Removed %d units with less than 1000 spikes and low similarity, %d units remaining\n', numel(irem), Nfilt);

% do easy merges
st4 = find_merges2(st3, Xsim, ops.fs,lagsim);

[~, ~, stidfin]   = unique(st4(:,2));
clustids = accumarray(stidfin, st3(:,2), [],@(x) {unique(x)});
mergesperunit = cellfun(@numel, clustids) - 1;
Nunitsf     = numel(clustids);
st4(:,2) = stidfin;

% instead of adding up templates, use the one from the first merger

final_temps = estas(:,:, cellfun(@(x) x(1), clustids));


% final_temps = zeros(numel(dt), Nchan, Nunitsf, 'single');
% simpleunitids = [clustids{mergesperunit == 0}];
% final_temps(:,:, mergesperunit == 0) = estas(:,:, simpleunitids);
% mergeIds = find(mergesperunit>0);
% 
% for ii = 1:numel(mergeIds)
%     currtemp = clustids{mergeIds(ii)};
%     currspk  = normfac(currtemp);
%     avgtemp  = reshape(estas(:,:,currtemp), [numel(dt)*Nchan, numel(currtemp)]) *currspk/sum(currspk);
%     final_temps(:,:,mergeIds(ii)) = reshape(avgtemp, [numel(dt),Nchan]);
% end





% recalculate similarities
[Xsim2,lagsim2] = compareSpikeTemplatesCMOS(final_temps);
Xsim2(isnan(Xsim2)) = 0;
Nfilt = size(Xsim2,1);

% do easy duplicate removal

isdupli = find_duplicates2(st4, Xsim2, ops.fs, final_temps, [allrez{1}.xc(:) allrez{1}.yc(:)]);

st4(ismembc(st4(:,2), find(isdupli)), :) = [];
% re-index templates
[ia, ~, ic] = unique(st4(:,2));
final_temps = final_temps(:,:,ia);
st4(:,2)   = ic;
Xsim2   = Xsim2(ia,ia);
lagsim2 = lagsim2(ia,ia);

Nfilt = max(ic);
fprintf('%d units remaining, good luck with the curation\n', Nfilt);

% calculate PC components
% for spikes that were missed, give them the average component 

% extract low-rank templates from elec images
[W, U, mu, UtU, nu] = decompose_dWU(ops, final_temps, ops.Nrank,allrez{1}.ops.kcoords);

% recalculate similarities based on low-rank decomposition



% edit W for phy
W               = cat(1, zeros(ops.nt0 - 2*ops.nt0min - 1, Nfilt, ops.Nrank), W);

% sort spike times
[~, isort] = sort(st4(:,1),'ascend');
st4 = st4(isort,:);

% sort pairwise template
rez = struct();
rez.simScore = Xsim2;
rez.W   = W;
rez.U   = U;
rez.mu  = mu;
rez.dWU = final_temps;
rez.ops = ops;
rez.xcoords = allrez{1}.xcoords;
rez.ycoords = allrez{1}.ycoords;
rez.Wrot    = allrez{1}.Wrot;
rez.st3     = st4;
rez.chunckInds = chunckInds;


[~, iNgsort] = sort(rez.simScore, 1, 'descend');
maskTT = zeros(size(iNgsort,1), 'single');
rez.iNeigh = iNgsort(1:rez.ops.nNeigh, :);
for i = 1:size(iNgsort,1)
    maskTT(rez.iNeigh(:,i),i) = 1;
end
rez.cProj = zeros(size(rez.st3,1),rez.ops.nNeigh,'single');

% sort best channels
[~, iNch]       = sort(abs(rez.U(:,:,1)), 1, 'descend');
maskPC          = zeros(rez.ops.Nchan, size(iNgsort,1), 'single');
rez.iNeighPC    = iNch(1:rez.ops.nNeighPC, :);
for i = 1:size(iNgsort,1)
    maskPC(rez.iNeighPC(:,i),i) = 1;
end



rez.cProjPC = zeros(size(rez.st3,1),3,rez.ops.nNeigh,'single');


end