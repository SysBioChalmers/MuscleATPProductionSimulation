close all
hold all
addpath('src1')
addpath('src2')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

color2 = [215 86 40]/256;

O2corr = false;
subjects = {'subject1';'subject2';'subject3';'subject4';'subject5'};
results = zeros(length(subjects),2);
resultsError = zeros(length(subjects),2);
resultsCplx1Max = zeros(length(subjects),1); 
resultsRQ = zeros(length(subjects),2); 
for i = 1:length(subjects)
    subplot(3,2,i)
    [delta, se, p, vO2max, intercept, RQ] = fitSubject(subjects{i}, O2corr);
    results(i,:) = delta';
    resultsError(i,:) = se';
    resultsCplx1Max(i) = intercept(2,1)/vO2max;
    resultsRQ(i,:) = RQ;
end

%%
subplot(3,2,6)
hold all

for i = 1:size(resultsRQ,1)
    plot([1 2], resultsRQ(i,:), 'ko-', 'markerfacecolor', 'k')
end

nonNormal = or(lillietest(resultsRQ(:,1)), lillietest(resultsRQ(:,2)));

if nonNormal
    [p h] = ranksum(resultsRQ(:,1), resultsRQ(:,2), 'tail','left');
else
    [h p] = ttest(resultsRQ(:,1), resultsRQ(:,2));
end
text(1.5, 0.75, sprintf('p=%2.2e',p))

text(0.3, median(resultsRQ(:,1)), sprintf('%2.2f', median(resultsRQ(:,1))))
text(1.3, median(resultsRQ(:,2)), sprintf('%2.2f', median(resultsRQ(:,2))))
xlim([0 3])
ylim([0.7 1])
figure()
hold all
xdisplace = rand(length(subjects),1)*0.2;
xdisplace = xdisplace-mean(xdisplace);

for i = 1:length(subjects)
   errorbar([1 2]+xdisplace(i), results(i,:), resultsError(i,:), 'ko-', 'markerfacecolor', 'k');
end

mresults = mean(results);
effectSize = mean(results(:,2)./results(:,1));
effectError = std(results(:,2)./results(:,1));
%plot([1 2], mresults, 'linewidth', 2, 'color', color2);
plot([0.5 2.2], mresults(1) * [1 1], '--', 'linewidth', 2, 'color', color2);
plot([0.5 2.2], mresults(2) * [1 1], '--', 'linewidth', 2, 'color', color2);
plot([2.2 2.2], mresults, '--', 'linewidth', 2, 'color', color2);


text(2.2, sum(mresults)/2, sprintf('%2.2fx\\pm%2.2f',effectSize,effectError))

xlim([0.5 3])
xticks([1 2])
xticklabels({'low', 'high'})
ylim([7 16])

p = ranksum(results(:,1),results(:,2));
text(1.5, 7.5, sprintf('p=%2.2f',p))
