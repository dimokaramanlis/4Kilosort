function uprojsize = get_uproj_hdd(rez)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ops = rez.ops;
Nchan = ops.Nchan;
Nbatch      = rez.temp.Nbatch;

NT  	= ops.NT;
batchstart = 0:NT:NT*Nbatch;
wPCA = ops.wPCA;

%uproj = zeros(2e6,  size(wPCA,2) * Nchan, 'single');
ops.GPU = 0;

fid   = fopen(ops.fproc, 'r');
fidWr = fopen(ops.uprojpath, 'W');

i0 = 0;

for ibatch = 1:ops.nskip:Nbatch

    offset = 2 * ops.Nchan*batchstart(ibatch);
    fseek(fid, offset, 'bof');
    dat = fread(fid, [NT ops.Nchan], '*int16');
 
    % move data to GPU and scale it
    dataRAW = gpuArray(dat);
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
    % find isolated spikes
    [row, col, mu] = isolated_peaks_new(dataRAW, ops);
    
    % find their PC projections (IN CPU)
    uS = get_PCproj(gather(dataRAW), gather(row), gather(col), wPCA, ops.maskMaxChannels,ops.nt0min);
    uS = permute(uS, [2 1 3]);
    uS = reshape(uS,numel(row), Nchan * size(wPCA,2));
    
%     if i0+numel(row)>size(uproj,1)
%         uproj(2e6 + size(uproj,1), 1) = 0;
%     end
    
    fwrite(fidWr, uS', 'single');
    
    i0 = i0 + numel(row);
    
end
fclose(fid);
fclose(fidWr);

uprojsize = [i0,  Nchan * size(wPCA,2)];
        
end

