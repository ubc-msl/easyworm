% Function that converts top graph plot into analysable data
function savetop
h = findobj(gca,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
datapoints(:,1) = transpose(x{2,1});
datapoints(:,2) = transpose(y{2,1});

[filename, pathname] = uiputfile('*.mat','Save Workspace Variables As');
newfilename = fullfile(pathname, filename);
save(newfilename, 'datapoints');
