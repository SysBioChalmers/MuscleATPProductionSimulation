close all

subjectData = [
    9.3830   11.4079
    8.5441   11.9159
    8.0908   12.6074
    9.2035   11.7251
   10.6617   14.2856
];

subjectError = [
    0.3575    0.1829
    0.1913    0.1055
    0.1825    0.1198
    0.2197    0.2070
    0.1161    0.3981];

resultsCplx1Max = [
    0.4210
    0.4058
    0.3938
    0.4045
    0.6247
    ];

slope1 = subjectData(:,1);
slope2 = subjectData(:,2);

subjectData = [slope1; slope2];
groups = [ones(5,1); 2*ones(5,1)];

groupNames = {'low W', 'high W'};
h = boxplot(subjectData, groups, 'Orientation',  'vertical',  'Colors', 'k', 'Widths', 1);
set(h(7,:),'Visible','off')
set(gca,'xticklabel', groupNames)  
ylim([0 17])
ylabel('dO2/dW')

color = [67 116 160; 151 185 224]/255;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    currentColor = color((mod(j,2)+1),:);
    a = patch(get(h(j),'XData'),get(h(j),'YData'), currentColor);
    uistack(a,'bottom');
end
%uistack(a,'bottom');
p = ranksum(subjectData(groups==1),subjectData(groups==2));
%[h p] = ttest(deltas(groups==1),deltas(groups==2));

text(1, 4, sprintf('p=%2.2e', p))
ylim([0 16])
color2 = [215 86 40]/256;
hold all
rng('default') % set random generator to same value each time
xvals = 0.5*rand(5,1);
xvals = xvals-mean(xvals);

for i = 1:length(slope1)
    tmp = plot([1 2] + xvals(i), [slope1(i) slope2(i)], '-', 'color', color2);
    tmp.Color(4) = 0.5;
end

scatter(1 + xvals, slope1, 20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
scatter(2 + xvals, slope2, 20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)


text(0.3, median(slope1), sprintf('%2.1f', median(slope1)))
text(1.3, median(slope2), sprintf('%2.1f', median(slope2)))

%%
figure()

h = boxplot(resultsCplx1Max, 'Orientation',  'vertical', 'Whisker', 5, 'Colors', 'k', 'Widths', 1);
set(h(7,:),'Visible','off')
set(gca,'xticklabel', '')  
ylim([0 17])
ylabel('complex I max/vO2 max')

color = [67 116 160; 151 185 224]/255;
h = findobj(gca,'Tag','Box');
a = patch(get(h,'XData'),get(h,'YData'), color(2,:));
uistack(a,'bottom');

text(1, 4, sprintf('p=%2.2e', p))
ylim([0 1])
xlim([0 2])
color2 = [215 86 40]/256;
hold all
rng('default') % set random generator to same value each time
xvals = 0.5*rand(5,1);
xvals = xvals-mean(xvals);

scatter(1 + xvals, resultsCplx1Max,20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)


text(0.3, median(resultsCplx1Max), sprintf('%2.1f', median(resultsCplx1Max)))


%%
addpath('src2')
Wvalues = [median(slope1) median(slope2)];
Ovalue = 2 * mlToMol(Wvalues)/3600;
ATPvalue = (1/27)/1000;
ATPvalue./Ovalue



