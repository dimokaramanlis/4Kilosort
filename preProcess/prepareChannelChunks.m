function [rez, Egroups] = prepareChannelChunks(rez)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf('Time %3.0f min. Extracting 3 PCs from data...\n', toc/60)
wPCA    = extractPCfromSnippets(rez, 3);
rez.ops.wPCA = wPCA;

% 
% 
% Nchuncks = ceil(rez.ops.Nchan/(Nebacth-Noverlap));
% 
% coords = [rez.xc rez.yc];



coords = [rez.xc(:) rez.yc(:)];
% irem = randperm(4096,2500);
% %coords(3000:end, :) = [];
% coords(irem, :) = [];

Nebacth = 800;
Noverlap = 200;
if size(coords,1) < Nebacth
    Egroups{1} = 1:size(coords,1);
    fprintf('Too few channels to split the data. Continuing.,,\n')
    return;
end

[~, pscores] = pca(coords);
[psorted, isort] = sort(pscores(:,1),'ascend');
Nchuncks = ceil((size(coords,1) + Noverlap)/Nebacth);
if Nchuncks < 3
    Nsize = ceil(size(coords,1)/2);
    Egroups{1} = isort(1:(Nsize+Noverlap/2));
    Egroups{2} = isort(Nsize-Noverlap/2:end);       
else
    Nsize = ceil((size(coords,1)-2*Noverlap)/Nchuncks);
    Egroups = cell(Nchuncks,1);
    
    Egroups{1} = isort(1:(Nsize+Noverlap));
    iend = Nsize+Noverlap;
    for ii = 2:Nchuncks
        istart      = iend-Noverlap/2+1;
        iend        = min(istart+Nsize+Noverlap,size(coords,1));
        Egroups{ii} = isort(istart:iend);
    end
    
end
for ii = 1:Nchuncks
    Egroups{ii} = sort(Egroups{ii},'ascend');
end
fprintf('Split the data into %d chunks. Continuing.,,\n',Nchuncks)

% colsall = [0.5 0.2 0.7; 0.3 0.6 0.9; 0.2 0.5 0.2; 0.8 0.3 0.3];
% clf; hold on;
% for ii = 1:Nchuncks
%     plot(coords(Egroups{ii,:},1)+randn(1)*10.1, coords(Egroups{ii,:},2)+randn(1)*10.1, 'o','Color',colsall(ii,:))
% end

% 
% Egroups{1} = isort(1:(Nsize+Noverlap));
% Egroups{Nchuncks} = isort(end-(Nsize+Noverlap-1):end);
% 
% elecrem   = psorted(Nsize+Noverlap/2+1:end-Nsize-Noverlap/2);
% isortrem = isort(Nsize+Noverlap/2+1:end-Nsize-Noverlap/2);
% 
% if Nchuncks-2>1
%     [centids,centroids] = kmeans(elecrem,Nchuncks-2,'Replicates',100);
%     [~,sortedcent] = sort(centroids, 'ascend');
%     for ii = 1:numel(centroids)
%         currinds      = isortrem(centids==sortedcent(ii));
%         [~, iminhigh] = min(abs(psorted-max(elecrem(centids==sortedcent(ii)))));
%         [~, iminlow] = min(abs(psorted-min(elecrem(centids==sortedcent(ii)))));
%         Egroups{1+ii} = [isort(iminlow-Noverlap/2:iminlow-1); currinds; isort(iminhigh+1:iminhigh+Noverlap/2)];
%     end
% else
%     Egroups{2} = isort(Nsize-Noverlap/2:end-(Nsize-Noverlap/2));
% 
% end

% colsall = [0.5 0.2 0.7; 0.3 0.6 0.9; 0.2 0.5 0.2; 0.8 0.3 0.3];
% clf;
% for ii = 1:Nchuncks
%     plot(coords(Egroups{ii,:},1)+randn(1)*0.1, coords(Egroups{ii,:},2)+randn(1)*0.1, 'o','Color',colsall(ii,:))
% end
% %%
% [~, idir] = max(mad(coords, 1));
% 
% [centids,centroids] = kmeans(coords(:,idir),Nchuncks,'Replicates',100);
% %Ncents = accumarray(centids, 1, [Nchuncks 1],@sum);
% %[~, isort] = sort(Ncents,'descend');
% 
% chunckInds = zeros(Nebacth, Nchuncks);
% for ii = 1:Nchuncks
%     [~, isort] = sort(abs(coords(:,idir) - centroids(ii)),'ascend');
%     chunckInds(:, ii) = sort(isort(1:Nebacth),'ascend');
% end


end