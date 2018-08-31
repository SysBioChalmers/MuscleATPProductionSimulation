clf
addpath('src2')
addpath('src3')

muscleMass = 40; %Kg
muscleWater = 0.8; 
blood = 5; %L
lactateReservoir = muscleMass*muscleWater + blood; %l
lactateTreshold = 240; %W
maxLactateFlux = 4; %mol/h

wattMax = 350; %W
wattSamples = [250 300 350];
wattSampleVO2 = [3230 3750 4400]; %ml/min

timeSpan = [0, 60/60];
interestingTime = [0 15/60];
vMax = 8;
basalFlux = 1;



watt = linspace(0,wattMax);
lactateFlux = 1*(watt>lactateTreshold);
lactateFlux(lactateFlux==1) = linspace(0,maxLactateFlux,sum(lactateFlux));
lactateConcentration = convertVlacToConc(lactateFlux + basalFlux, 0, vMax);
basalAmount = lactateReservoir * lactateConcentration(1);

lactateSamples = interp1(watt,lactateFlux,wattSamples);



tIn = linspace(timeSpan(1),timeSpan(2));
lactateResults = zeros(length(wattSamples), length(tIn));


for i = 1:length(wattSamples)
    odParam = @(T,X) lactateOde(T, X, lactateReservoir, basalFlux, vMax, lactateSamples(i));
    [t,yv] = ode23s(odParam, tIn, basalAmount);
    lactateResults(i,:) = yv/lactateReservoir;
end


vO2results = zeros(length(wattSamples), length(tIn));

for i = length(wattSamples):-1:1
    vLactateOut = mmLactate(lactateResults(i,:), vMax)-basalFlux;
    vLactateIn = lactateSamples(i);
    netLactate = vLactateIn-vLactateOut;
    netOxygen = netLactate * 6 * 1.5/30;
    oxygenDeficit = molToMl(netOxygen); %ml/mol lactate
    vO2results(i,:) = wattSampleVO2(i) - oxygenDeficit;
end


%%
subplot(2,2,1)
plot(watt, lactateFlux, 'linewidth', 3)
title('steady state lactate flux')
xlabel('W')
ylabel('vLactate')

subplot(2,2,2)
hold all
plot(watt, lactateConcentration, 'linewidth', 3)
title('steady state lactate concentration')
xlabel('W')
ylabel('lactate concentration [mM]')
%plot(wattSamples, lactateResults(:,end), 'o')

subplot(2,2,3)
hold all

wattLabels = [];
for i = size(lactateResults,1):-1:1
    plot(60*tIn, lactateResults(i,:), 'linewidth', 3)
    wattLabels{i} = sprintf('%2.0f Watt', wattSamples(i));    
end
wattLabels = wattLabels(end:-1:1);
legend(wattLabels)
legend boxoff

for i = 1:size(lactateResults,1)
    plot(60*interestingTime, (lactateResults(i,end)*[1 1]), 'k-')
end

xlabel('time [min]')
ylabel('lactate concentration [mmol/l]')
title('Lactate dynamics')
xlim(60*interestingTime)

% plot(60*6 * [1 1], lactateResults(i,end)*[0 1], 'k--')

subplot(2,2,4)
hold all
for i = size(vO2results,1):-1:1
    plot(60*tIn, vO2results(i,:), 'linewidth', 3)
end

% legend(wattLabels, 'location', 'se')
% legend boxoff

for i = 1:size(vO2results,1)
    plot(60*interestingTime, vO2results(i,end)*[1 1], 'k-')
end
% plot(60*6 * [1 1], vO2results(i,end)*[0.7 1], 'k--')

xlim(60*interestingTime)
xlabel('time [min]')
ylabel('vO2 [ml/min]')
title('Oxygen dynamics')

