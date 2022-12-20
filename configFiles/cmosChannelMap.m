

function channelmap = cmosChannelMap(ops, electrodedist, savingpath, varargin)
%CREATE252CHANNELMAPFILE
%Doesn't reorder channels, the actual reordering is happening through convertMcdToRawBinaryCAR

if nargin>4
    excludeanalog = varargin{1};
else
    excludeanalog = true;
end

% h5filenames = dir([ops.root,filesep,'*.brw']);
% [~, reindex]=sort(str2double(regexp(({h5filenames(:).name}),'\d+','match','once')));
% h5filenames={h5filenames(reindex).name}'; 
% 
% h5pathname = [ops.root,filesep,h5filenames{1}]; %get mcd path
% chs        = h5read(h5pathname, '/3BRecInfo/3BMeaStreams/Raw/Chs');
% 
% xcoords = electrodedist * double(chs.Col(:)); 
% ycoords = electrodedist * double(chs.Row(:));

[ycoords,xcoords] = meshgrid(0:63,63:-1:0);

ycoords = electrodedist * ycoords;
xcoords = electrodedist * xcoords;

if excludeanalog
    xcoords(1:2)=[];
    ycoords(1:2)=[];
end
chanMap = 1:numel(xcoords); 
chanMap0ind = chanMap - 1;
kcoords = ones(size(xcoords));
connected = true(size(xcoords));

%-------------------------------------------------------------------------- 
if nargin > 2
    if not(exist(fullfile(savingpath, 'chanMap.mat'),'file'))
        save(fullfile(savingpath, 'chanMap.mat'), 'chanMap',...
            'chanMap0ind','connected', 'xcoords', 'ycoords', 'kcoords');
    else
        disp('we found channel map from previous run!');
    end
end

channelmap.chanMap = chanMap;
channelmap.chanMap0ind = chanMap0ind;
channelmap.connected = connected;
channelmap.xcoords = xcoords;
channelmap.ycoords = ycoords;
channelmap.kcoords = kcoords;
channelmap.numelectrodes = numel(xcoords);
channelmap.electrodedist = electrodedist;

end