addpath('sampleData')
data = loadSampleData('sampleData/allData.txt');
subjects = unique(data.text(:,1));
ergo = 'KE';
oxygenLevel = '1';

hold all

deltaLow = zeros(length(subjects),1);
deltaHigh = zeros(length(subjects),1);
allW = [];
allO2 = [];
allSubjects = [];

for i = 1:length(subjects)
    wData = extractData(data, {'Subject', 'Ergo', 'Treat'}, {subjects{i}, ergo, oxygenLevel}, 'LOAD_W');
    O2Data =  extractData(data, {'Subject', 'Ergo', 'Treat'}, {subjects{i}, ergo, oxygenLevel}, 'VO2');
    
    deltaLow(i) = (O2Data(2)-O2Data(1))/(wData(2)-wData(1));
    deltaHigh(i) = (O2Data(end-1)-O2Data(end))/(wData(end-1)-wData(end));
    
    allW = [allW; wData/max(wData)];
    allO2 = [allO2; O2Data/max(O2Data)];
    allSubjects = [allSubjects; i*ones(length(wData),1)];
end

deltas = [deltaLow; deltaHigh];
groups = [ones(length(deltaLow),1); 2 * ones(length(deltaLow),1)];


groupNames = {'low W', 'high W'};
h = boxplot(deltas, groups, 'Orientation',  'vertical',  'Colors', 'k', 'Widths', 1, 'OutlierSize',4, 'Symbol', 'ko');
set(gca,'xticklabel', groupNames)  
ylim([0 35])
ylabel('dO2/dW')

color = [67 116 160; 151 185 224]/255;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    currentColor = color((mod(j,2)+1),:);
    a = patch(get(h(j),'XData'),get(h(j),'YData'), currentColor);
    uistack(a,'bottom');
end
uistack(a,'bottom');
p = ranksum(deltas(groups==1),deltas(groups==2));
text(1, 32, sprintf('p=%2.4f', p))

color2 = [215 86 40]/256;
hold all
xvals = 0.2*rand(length(deltaLow),1);
xvals = xvals-mean(xvals);
scatter(1 + xvals, deltaLow, 20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
scatter(2 + xvals, deltaHigh, 20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)

% for i = 1:length(deltaLow)
%     plot([1 2]+xvals(i), [deltaLow(i) deltaHigh(i)])
% end

