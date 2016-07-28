function getDataFromPlot(plotTitle)
% open the figure
figName = sprintf('%s.fig',plotTitle);
open(figName);

% get the figure information
h = gcf;                                %current figure handle
axesObjs = get(h, 'Children');          %axes handles
% dataObjs = axesObjs.Children;
dataObjs = get(axesObjs, 'Children');   %handles to low-level graphics objects in axes

% get the data points
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
zdata = get(dataObjs, 'ZData');

% save the data pulled from the plot
saveName = sprintf('%s.mat', plotTitle);
save(saveName, 'xdata', 'ydata', 'zdata');
end