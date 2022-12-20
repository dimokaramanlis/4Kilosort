function [] = plotSpikeTemplate(spk_template, coords, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    col = 'k';
else 
    col = varargin{1};
end

plotfac = 0.45;
spk_template = squeeze(spk_template);

xvals = coords(:,1); yvals = coords(:,2);

alldist = sqrt( (xvals - xvals').^2 + (yvals - yvals').^2 );
mindist = min(alldist(alldist(:)>0));

tvec = linspace(-mindist*plotfac, mindist*plotfac, size(spk_template,1));
ymax = max(abs(spk_template(:)));

xdata = xvals' + tvec';
xdataall = [xdata; nan(1,size(xdata,2))];
xdataall = repmat(xdataall, [1 size(spk_template,3)]);

ydata = yvals' + spk_template*plotfac*mindist/ymax;
ydataall = [ydata(:,:); nan(1,size(ydata,2)*size(ydata,3))];

% axis equal; 
% xlim([min(xvals)-mindist max(xvals)+mindist]); 
% ylim([min(yvals)-mindist max(yvals)+mindist]); 
l = line(xdataall(:),ydataall(:),'LineStyle','-','Marker','none', 'Color',col);


end

