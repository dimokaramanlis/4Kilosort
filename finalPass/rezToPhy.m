
function  rezToPhy(rez, savePath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%
% spikeTimes will be in samples, not seconds


outputs = {'amplitudes.npy', 'channel_map.npy', 'channel_positions.npy', 'pc_features.npy', ...
           'pc_feature_ind.npy', 'similar_templates.npy', 'spike_clusters.npy', 'spike_templates.npy', ...
           'spike_times.npy', 'templates.npy', 'templates_ind.npy', 'template_features.npy', ...
           'template_feature_ind.npy', 'whitening_mat.npy', 'whitening_mat_inv.npy'};

fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
    fname = fs(i).name;
    % don't delete .npy files which have nothing to do with us
    if find(strcmp(fname, outputs))
        delete(fullfile(savePath, fname));
    end
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end

% clean up rez duplicates
% [~,uniqueIdxs,~] = unique(rez.st3,'rows');
% rez.st3 = rez.st3(uniqueIdxs,:);
% rez.cProj = rez.cProj(uniqueIdxs,:);
% rez.cProjPC = rez.cProjPC(uniqueIdxs,:,:);

spikeTimes = uint64(rez.st3(:,1));
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end

amplitudes = single(rez.st3(:,3));

Nchan = rez.ops.Nchan;

xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;
Nfilt = size(rez.dWU,3);
templates = zeros(Nchan, nt0, Nfilt, 'single');
for iNN = 1:Nfilt
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatureInds = uint32(rez.iNeigh);
pcFeatureInds       = uint32(rez.iNeighPC);

whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
whiteningMatrixInv = whiteningMatrix^-1;


if ~isempty(savePath)
    
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy')); clear spikeTimes;
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy')); clear amplitudes;
    writeNPY(templates, fullfile(savePath, 'templates.npy'));   clear templates;
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
    writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));
    
    writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    writeNPY(rez.cProj, fullfile(savePath, 'template_features.npy'));
    if rez.ops.lowmem
         savePrincipalComponentsKS1(rez, fullfile(savePath, 'pc_features.npy'))
    else
        writeNPY(rez.cProjPC, fullfile(savePath, 'pc_features.npy'));
    end

    writeNPY(whiteningMatrix,    fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
    
    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end
    
     %make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');
        
        [~, fname, ext] = fileparts(rez.ops.fbinary);
        
        fprintf(fid,['dat_path = ''',fname ext '''\n']);
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez.ops.fs,1)
            fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
        end
        fprintf(fid,'hp_filtered = False\n');
        fprintf(fid,'template_scaling = 20.0\n');
        fclose(fid);
    end
    
end
