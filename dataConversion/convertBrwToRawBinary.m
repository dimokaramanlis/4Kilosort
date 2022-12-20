function [bininfo] = convertBrwToRawBinary(ops, targetpath)
% first attempt for cmos writing

%--------------------------------------------------------------------------
tic;
% get brw filenames
h5filenames = dir([ops.root,filesep,'*.brw']);
[~, reindex]=sort(str2double(regexp(({h5filenames(:).name}),'\d+','match','once')));
h5filenames={h5filenames(reindex).name}'; Nfiles=numel(h5filenames);
%--------------------------------------------------------------------------
% get information about the recording time
% config files for dataloading

stimdata = cell(numel(h5filenames),1);

%--------------------------------------------------------------------------
imcd=1;
h5pathname = [ops.root,filesep,h5filenames{imcd}]; %get mcd path

fs       = h5read(h5pathname, '/3BRecInfo/3BRecVars/SamplingRate');
bitDepth = h5read(h5pathname, '/3BRecInfo/3BRecVars/BitDepth');
minVolt  = h5read(h5pathname, '/3BRecInfo/3BRecVars/MinVolt');
maxVolt  = h5read(h5pathname, '/3BRecInfo/3BRecVars/MaxVolt');
isinv    = h5read(h5pathname, '/3BRecInfo/3BRecVars/SignalInversion');
chs      = h5read(h5pathname, '/3BRecInfo/3BMeaStreams/Raw/Chs');
NchanTOT = length(chs.Row);

newRange=2^15*[-1 1]; multFact=range(newRange)/(maxVolt-minVolt);
satVal  = 2^single(bitDepth) - 1;
replVal = uint16(-minVolt*(satVal+1)/(maxVolt-minVolt));
%--------------------------------------------------------------------------
stimsamples = zeros(numel(h5filenames),1);
for imcd = 1:numel(h5filenames)
    h5pathname = fullfile(ops.root, h5filenames{imcd}); %get mcd path
    info = h5info(h5pathname, '/3BData/Raw');
    stimsamples(imcd) = info.Dataspace.Size(1)/NchanTOT;
end
%--------------------------------------------------------------------------
bininfo.stimsamples = stimsamples;
bininfo.NchanTOT = NchanTOT - 2;
bininfo.fs = fs;
fprintf('Total length of recording is %2.2f min...\n',sum(stimsamples)/fs/60);
%--------------------------------------------------------------------------
analogpath = fullfile(ops.root,'ks_sorted','analog');
if ~exist(analogpath, 'dir')
    mkdir(analogpath);
end
%--------------------------------------------------------------------------
fprintf('Saving .brw data as binary .dat...\n');

maxSamples = 20e5;% chunk size

fidOut = fopen(targetpath, 'W'); %using W (capital), makes writing ~4x faster

NsampleSaturations = zeros(NchanTOT, 1, 'single');

for ifile = 1:Nfiles
    
    h5pathname = fullfile(ops.root, h5filenames{ifile}); %get mcd path
    Nsamples   = stimsamples(ifile);
    Nchunks    = ceil(Nsamples/maxSamples);

    [~,analogpathname] = fileparts(h5filenames{ifile});
    analogpathname     = fullfile(analogpath, [analogpathname '.dat']);
    fidanalog          = fopen(analogpathname, 'W'); %using W (capital), makes writing ~4x faster

    for ichunk = 1:Nchunks
        offset        = max(0, (maxSamples * (ichunk-1)));
        sampstoload   = min(Nsamples-offset, maxSamples);
        
        dat = h5read(h5pathname, '/3BData/Raw', offset * NchanTOT + 1, sampstoload * NchanTOT);
        
        % KEEP TRACK OF SATURATED VALUES
        indssat = (dat == satVal); 
        NsampleSaturations = NsampleSaturations + sum(reshape(indssat, NchanTOT, sampstoload), 2);
        
        dat(indssat) = replVal; % replace saturated values
        
        dat = convertADC(dat, maxVolt, minVolt, bitDepth);
        dat = reshape(dat, [NchanTOT, sampstoload]);
        %dat = int16(dat(3:end, :) * multFact); % skip first two useless channels
        dat = int16(dat * multFact); % skip first two useless channels
        
        fwrite(fidanalog,  dat(2,:), 'int16');
        fwrite(fidOut, dat(3:end,:), 'int16');
    end

    fclose(fidanalog); %close analog path

    %report status
    fprintf('Time %3.0f min. Brw files processed %d/%d \n', toc/60, ifile,Nfiles);   
end
fclose(fidOut);


bininfo.saturationPcnt = NsampleSaturations/sum(stimsamples);
%--------------------------------------------------------------------------
end

function volts = convertADC(value, MaxVolt, MinVolt, BitDepth)
    % volts = convertADC(value, MaxVolt, MinVolt, BitDepth)
    span = MaxVolt - MinVolt;
    n = 2^double(BitDepth); % number of possible values
    
   volts = MinVolt + single(value) * span/n;
end
