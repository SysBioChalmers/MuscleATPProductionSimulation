addpath('src2')
addpath('src3')

muscleMass = 40;
muscleWater = 0.8;
blood = 5;
lactateReservoir= 40*0.8 + 5;
lactateTreshold = 250;
wattMax = 350;
timeSteps = 3;

timeSpan = [0, 20*60];
interestingTime = [0 6*60];
vMax = 9;
basalFlux = 1;



watt = linspace(0,wattMax);
lactateFlux = 1*(watt>lactateTreshold);
lactateFlux(lactateFlux==1) = linspace(0,4,sum(lactateFlux));
lactateConcentration = convertVlacToConc(lactateFlux, basalFlux, vMax);
basalConcentration = lactateConcentration(1);

wattSamples = linspace(lactateTreshold, wattMax, timeSteps);
lactateSamples = interp1(watt,lactateFlux,wattSamples);

tIn = linspace(timeSpan(1),timeSpan(2));
lactateResults = zeros(length(wattSamples), length(tIn));


for i = 1:length(wattSamples)
    odParam = @(T,X) lactateOde(T, X, lactateReservoir, basalFlux, lactateSamples(i));
    [t,yv] = ode23s(odParam, tIn, basalConcentration);
    lactateResults(i,:) = yv;
end


baseVO2 = 500; %ml/min;
vO2results = zeros(length(wattSamples), length(tIn));
for i = length(wattSamples):-1:1
    vLactate = mmLactate(lactateResults(i,:));
    nonLactateVO2 = baseVO2 + wattSamples(i)*10; %ml/min/W
    additionalVO2 = 1000*molToMl(vLactate * 3/3600); %ml/mol lactate
    vO2results(i,:) = nonLactateVO2 + additionalVO2;
end




subplot(2,2,1)
plot(watt, lactateFlux, 'linewidth', 3)
title('steady state lactate flux')
xlabel('W')
ylabel('vLactate')

subplot(2,2,2)
plot(watt, lactateConcentration, 'linewidth', 3)
title('steady state lactate concentration')
xlabel('W')
ylabel('[lactate]')

subplot(2,2,3)
hold all

wattLabels = [];
for i = size(lactateResults,1):-1:1
    plot(tIn, lactateResults(i,:), 'linewidth', 3)
    wattLabels{i} = sprintf('%2.0f Watt', wattSamples(i));    
end
wattLabels = wattLabels(end:-1:1);
legend(wattLabels)
legend boxoff

for i = 1:size(lactateResults,1)
    plot(interestingTime, (lactateResults(i,end)*[1 1]), 'k--')
end

xlabel('time [s]')
ylabel('[lactate]')
title('Lactate dynamics')
xlim(interestingTime)

subplot(2,2,4)
hold all
for i = size(vO2results,1):-1:1
    plot(tIn, vO2results(i,:), 'linewidth', 3)
end

legend(wattLabels, 'location', 'se')
legend boxoff

for i = 1:size(vO2results,1)
    plot(interestingTime, vO2results(i,end)*[1 1], 'k--')
end


xlim(interestingTime)
xlabel('time [s]')
ylabel('vO2')
title('Oxygen dynamics')