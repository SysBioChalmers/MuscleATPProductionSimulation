clf
addpath('src2')
addpath('src3')

muscleMass = 40; %Kg
muscleWater = 0.8; 
blood = 5; %L
lactateReservoir = muscleMass*muscleWater + blood; %l
%lactateEfficency = mean([(29/32) (19/22)]); <-theoretical
lactateEfficency = 1.35/1.57; %empirical

timeSpan = [0/3600, 70/60];
interestingTime = [-3/60 12/60];
vMax = 6.6;
basalFlux = mmLactate(0.9, vMax);

%estimate lactate flux:
lactateFluxIn = 1.8;
basalAmount = 0.9 * lactateReservoir;

tIn = linspace(timeSpan(1),timeSpan(2),500);

odeParam = @(T,X) lactateOde(T, X, lactateReservoir, basalFlux, vMax, lactateFluxIn);
[t,yv] = ode23s(odeParam, tIn, basalAmount);

empiricalVO2 = 1350 + 1800 * (1-exp(-(60*t-(9.3/60))/(24.9/60)));
empiricalVO2(t<9.3/3600) = 1345;
lactateResults = yv/lactateReservoir;

vLactateOut = mmLactate(lactateResults, vMax)-basalFlux;
oxygenRequirement = molToMl(vLactateOut * 3); %ml/mol lactate
lactateFraction = oxygenRequirement./empiricalVO2;

%modifiedVO2 = empiricalVO2 .* (1-lactateFraction) + empiricalVO2 .* lactateFraction/lactateEfficency;
modifiedVO2 = empiricalVO2./((1-lactateFraction) + lactateEfficency .* lactateFraction);


%%
subplot(2,2,1)
hold all

plot(60*tIn, lactateResults, 'linewidth', 3)


lactateData = [
    -2.5    0.9
    4.5     3.9
    8.8     4.3];

hold all
scatter(lactateData(1:3,1), lactateData(1:3,2), 'filled');

plot(60*interestingTime, (lactateResults(end)*[1 1]), '--', 'color', [0    0.4470    0.7410])
plot(60*interestingTime, (lactateData(1,2)*[1 1]), '--', 'color', [0.8500    0.3250    0.0980])

ylim([0 10])


xlabel('time [min]')
ylabel('lactate concentration [mmol/l]')
%title('Lactate dynamics')
xlim(interestingTime*60)

% plot(60*6 * [1 1], lactateResults(i,end)*[0 1], 'k--')

subplot(2,2,2)
hold all
plot(60*tIn, modifiedVO2, 'linewidth', 3)

% legend(wattLabels, 'location', 'se')
% legend boxoff

%%

data = [
    -2.5 1.35   0.05
    3   3.23    0.1 %con
    10  3.48	0.1];

data(:,2:3) = data(:,2:3) * 1000;

hold all
scatter(data(:,1), data(:,2), 'filled');

errorbar(data(:,1), data(:,2), data(:,3), 'k.');

plot(60*interestingTime, modifiedVO2(end)*[1 1],  '--', 'color', [0    0.4470    0.7410])
plot(60*interestingTime, (data(1,2)*[1 1]), '--', 'color', [0.8500    0.3250    0.0980])

empiricalVO2 = 1345 + 1767 * (1-exp(-(60*t-(9.3/60))/(24.9/60)));
empiricalVO2(t<9.3/3600) = 1345;
plot(t*60, empiricalVO2,'k--');

xlim(interestingTime*60)
ylim([1000 4000])
xlabel('time [min]')
ylabel('vO2 [ml/min]')
%title('Oxygen dynamics')


% figure
% plot(t*60,lactateFraction)
% xlim([0 15])

basalAmount = 14.0 * lactateReservoir;

tIn = linspace(timeSpan(1),timeSpan(2),500);

odeParam = @(T,X) lactateOde(T, X, lactateReservoir, basalFlux, vMax, lactateFluxIn);
[t,yv] = ode23s(odeParam, tIn, basalAmount);

empiricalVO2 = 1350 + 1800 * (1-exp(-(60*t-(9.3/60))/(24.9/60)));
empiricalVO2(t<9.3/3600) = 1345;
lactateResults = yv/lactateReservoir;

vLactateOut = mmLactate(lactateResults, vMax);
oxygenRequirement = molToMl(vLactateOut * 3); %ml/mol lactate
lactateFraction = oxygenRequirement./empiricalVO2;
lactateFraction(lactateFraction>1) = 1;

modifiedVO2 = empiricalVO2./((1-lactateFraction) + lactateEfficency .* lactateFraction);

%%
subplot(2,2,3)
hold all

plot(60*tIn, lactateResults, 'linewidth', 3)


lactateData = [
    -2.5	14.0
    4.5     11.4
    8.8     9.6
    ];

hold all
scatter(lactateData(1:3,1), lactateData(1:3,2), 'filled');

plot(60*interestingTime, (lactateResults(end)*[1 1]), '--', 'color', [0    0.4470    0.7410])
plot(60*interestingTime, (lactateData(1,2)*[1 1]), '--', 'color', [0.8500    0.3250    0.0980])

ylim([0 20])


xlabel('time [min]')
ylabel('lactate concentration [mmol/l]')
%title('Lactate dynamics')
xlim(interestingTime*60)

% plot(60*6 * [1 1], lactateResults(i,end)*[0 1], 'k--')

subplot(2,2,4)
hold all
plot(60*tIn, modifiedVO2, 'linewidth', 3)

% legend(wattLabels, 'location', 'se')
% legend boxoff



data = [
    -2.5	1.57	0.04
    3	3.66	0.1
    10	3.65	0.11
    ];

data(:,2:3) = data(:,2:3) * 1000;

hold all
scatter(data(:,1), data(:,2), 'filled');

errorbar(data(:,1), data(:,2), data(:,3), 'k.');

plot(60*interestingTime, modifiedVO2(end)*[1 1],  '--', 'color', [0    0.4470    0.7410])
plot(60*interestingTime, (data(1,2)*[1 1]), '--', 'color', [0.8500    0.3250    0.0980])

empiricalVO2 = 1574 + 1987 * (1-exp(-(60*t-(9.6/60))/(22.3/60)));
empiricalVO2(t<9.6/3600) = 1574;
plot(t*60, empiricalVO2,'k--');

xlim(interestingTime*60)
ylim([1000 4000])
xlabel('time [min]')
ylabel('vO2 [ml/min]')
%title('Oxygen dynamics')

