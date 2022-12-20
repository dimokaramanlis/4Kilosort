function wPCA = extractPCfromSnippets(rez, nPCs)

ops = rez.ops;

% Nchan 	= ops.Nchan;

Nbatch      = rez.temp.Nbatch;

NT  	= ops.NT;
batchstart = 0:NT:NT*Nbatch;

mind = diff(rez.yc);
mind = min(mind(mind>0));
c2c = sqrt((rez.xc-rez.xc').^2 + (rez.yc-rez.yc').^2)/mind;

% extract the PCA projections
CC = zeros(ops.nt0);
fid = fopen(ops.fproc, 'r');

for ibatch = 1:100:Nbatch

    offset = 2 * ops.Nchan*batchstart(ibatch);
    fseek(fid, offset, 'bof');
    dat = fread(fid, [NT ops.Nchan], '*int16');

    % move data to GPU and scale it
    dataRAW = gpuArray(dat);
 
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
    
    % find isolated spikes
    [row, col, mu] = isolated_peaks_new(dataRAW, ops,c2c);
    
    clips = get_SpikeSample(dataRAW, row, col, ops, 0);
    
    c = sq(clips(:, :));
    CC = CC + gather(c * c')/1e3;
    
end
fclose(fid);

[U Sv V] = svdecon(CC);

wPCA = U(:, 1:nPCs);

wPCA(:,1) = - wPCA(:,1) * sign(wPCA(ops.nt0min,1));
