close all
hold all

factor = 5 * 60*60 * 10^-12 * 10^3 * 10^3; %per hour per gdw 

data = factor*[
    33.45	2.11
    76.41	4.58
    105.99	5.63
    ];

exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;

hold all     
bar(data(:,1), 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
errorbar(1:3,data(:,1), data(:,2),'k.')
set(gca, 'XTick', [1 2 3])
set(gca, 'XTickLabel', {'Complex I', 'Complex II', 'Max Capacity'})
ylabel('vO2 [mmol/gdw/h]')
ylim([0, 3])
xtickangle(45)
set(findall(gcf,'-property','FontSize'),'FontSize',15)