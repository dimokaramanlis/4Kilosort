function  deleteBrwFiles(ops)
%UNTITLED Summary of this function goes here

h5filenames = dir([ops.root,filesep,'*.brw']);
[~, reindex]=sort(str2double(regexp(({h5filenames(:).name}),'\d+','match','once')));
h5filenames={h5filenames(reindex).name}'; Nfiles=numel(h5filenames);

for ifile = 1:Nfiles
    h5pathname = fullfile(ops.root, h5filenames{ifile}); %get mcd path
    delete(h5pathname);
end

end