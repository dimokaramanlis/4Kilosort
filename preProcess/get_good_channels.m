function igood = get_good_channels(ops, chanMap)

Nbatch = ops.Nbatch;
twind = 0;
NchanTOT = ops.NchanTOT;
NT = ops.NT;
Nchan = numel(chanMap);
lrange = getOr(ops, {'long_range'}, [30 6]);

% load data into patches, filter, compute covariance
if isfield(ops,'fslow') && ops.fslow<ops.fs/2 && ops.filter
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end


fid = fopen(ops.fbinary, 'r');
% irange = [NT/8:(NT-NT/8)];

ibatch = 1;
ich    = gpuArray.zeros(5e4,1, 'int16');
maxall = zeros(Nbatch, NchanTOT);
k = 0;
ttime = 0;

while ibatch<=Nbatch
    offset = twind + 2*NchanTOT*NT* (ibatch-1);
    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NT], '*int16');
        
    if isempty(buff)
        break;
    end
    
    datr    = gpufilter(buff, ops, chanMap); % apply filters and median subtraction
    
    maxall(ibatch, :) = max(datr,[], 1);
    % determine any threshold crossings
    %datr = datr./std(datr,1,1);
    currmad = mad(datr, 1, 1);
    datr = datr./currmad;
    
    mdat = my_min(datr, lrange(1), 1);

    ind = find(datr<mdat+1e-3 & datr<ops.madTh);
    [xi, xj] = ind2sub(size(datr), ind);
    xj(xi<ops.nt0 | xi>NT-ops.nt0) = [];

    if k+numel(xj)>numel(ich)
        ich(2*numel(ich)) = 0;        
    end        
    ich(k + [1:numel(xj)]) = xj;
    
    k = k + numel(xj);
    
    ibatch = ibatch + ceil(Nbatch/100); % skip every 100 batches
    ttime = ttime + size(datr,1)/ops.fs; % keep track of total time where we took spikes from    
end
fclose(fid);

maxel   = mean(maxall);
inormal = maxel < median(maxel) + 40*mad(maxel,1);

ich = ich(1:k);

nc = histcounts(ich, .5 + [0:Nchan]); % count how many spikes each channel got
nc = nc/ttime; % divide by total time to get firing rate

%igood = nc>=getOr(ops, 'minfr_goodchannels', .1);

igood = nc>=getOr(ops, 'minfr_goodchannels', .1) & ops.bininfo.saturationPcnt(3:end)' < ops.min_saturpcnt & inormal;

% chanMap = chanMap(nc>.1);
fprintf('found %d threshold crossings in %2.2f minutes of data \n', k, ttime/60)
fprintf('found %d bad channels \n', sum(~igood))
fprintf('removed %d channels with high amp events\n', sum(~inormal))


