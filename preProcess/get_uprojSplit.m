function [uproj] = get_uprojSplit(rez,chunckInds)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ops = rez.ops;
Nchan = numel(chunckInds);
Nbatch      = rez.temp.Nbatch;

NT  	= ops.NT;
batchstart = 0:NT:NT*Nbatch;
wPCA = ops.wPCA;

uproj = zeros(2e6,  size(wPCA,2) * Nchan, 'single');
ops.GPU = 1;
Nmax = 12000;

fid = fopen(ops.fproc, 'r');

mind = diff(rez.yc(chunckInds));
mind = min(mind(mind>0));
c2c = sqrt((rez.xc(chunckInds)-rez.xc(chunckInds)').^2 + (rez.yc(chunckInds)-rez.yc(chunckInds)').^2)/mind;

i0 = 0;

for ibatch = 1:ops.nskip:Nbatch

    offset = 2 * ops.Nchan*batchstart(ibatch);
    fseek(fid, offset, 'bof');
    dat = fread(fid, [NT ops.Nchan], '*int16');
    dat = dat(:, chunckInds);

    % move data to GPU and scale it
    if ops.GPU
        dataRAW = gpuArray(dat);
    else
        dataRAW = dat;
    end
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
    % find isolated spikes
    [row, col, mu] = isolated_peaks_new(dataRAW, ops, c2c);
    
    if numel(row)>Nmax
        iuse = randperm(numel(row),Nmax);
        row = row(iuse);
        col = col(iuse);
        mu  = mu(iuse);
    end
    
    % find their PC projections
    uS = get_PCproj(dataRAW, row, col, wPCA, ops.maskMaxChannels,ops.nt0min, c2c);
    %uS = permute(uS, [2 1 3]); % done permutation before
    uS = reshape(uS,numel(row), Nchan * size(wPCA,2));
    
    if i0+numel(row)>size(uproj,1)
        uproj(2e6 + size(uproj,1), 1) = 0;
    end
    
    uproj(i0 + (1:numel(row)), :) = gather_try(uS);
    i0 = i0 + numel(row);
    
end
fclose(fid);

uproj(i0+1:end, :) = [];

        
end

