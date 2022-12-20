function Wrot = get_whitening_matrix(rez)

ops = rez.ops;
Nbatch = ops.Nbatch;
twind = 0;
NchanTOT = ops.NchanTOT;
NT = ops.NT;
NTbuff = ops.NTbuff;
chanMap = ops.chanMap;
Nchan = rez.ops.Nchan;
xc = rez.xc;
yc = rez.yc;

% load data into patches, filter, compute covariance
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2 && ops.filter
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end

fid = fopen(ops.fbinary, 'r');
CC = gpuArray.zeros( Nchan,  Nchan, 'single');


% irange = [NT/8:(NT-NT/8)];

ibatch = 1;
while ibatch<=Nbatch    
    offset = max(0, twind + 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NTbuff], '*int16');
        
    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end
    
    datr    = gpufilter(buff, ops, chanMap); % apply filters and median subtraction
    
    CC        = CC + (datr' * datr)/NT;    
    
    ibatch = ibatch + ops.nSkipCov;
end
CC = CC / ceil((Nbatch-1)/ops.nSkipCov);

fclose(fid);
fprintf('Channel-whitening filters computed. \n');

if ops.whiteningRange<Inf
    ops.whiteningRange = min(ops.whiteningRange, Nchan);
    Wrot = whiteningLocal(gather_try(CC), yc, xc, ops.whiteningRange);
else
    [E, D] 	= svd(CC);
    D       = diag(D);
%     eps 	= mean(D); %1e-6;
    eps 	= 1e-6;
    f  = mean((D+eps) ./ (D+1e-6));
%     fprintf('%2.2f ', f)
    Wrot 	= E * diag(f./(D + eps).^.5) * E';
end
Wrot    = ops.scaleproc * Wrot;
