load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;


%Random sampling settings
nrOfSamples = 100;
nrOfTimesteps = 10;
pertubationLevel = 0.2; %+- 20%


%Parameterset
dwMuscle = 20*(1-0.792); %kg muscle
muscleRatio = 0.55; %type1
vO2perDryweight = 1.38 * 2.6; %Specific activity Complex IV
complex1Ratio = 33.45/76.41;
m2Efficency = 0.5;
vO2max = 1000;
internalO2 = 10;
maintainance = 5.5;
internalWork = 2;
peripheralFA = 0.004;
peripheralLactateCapacity = 2;
FAFactor = 0.8*16.55/76.41;
type1tresh = 9; %type 1 is activated first, 10 mol ATP corresponds to 20% of Wmax
type2target = 0.55; %type 2 dominates at the highest W rates

%Make reference condition
model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);
settings = setSimulationSettings(nrOfTimesteps, false);

%Reactions to plot
readouts = getTransport(model, {'O2', 'CO2', 'glycogen'}, 'sb', 's');

%apply oxygen constraint
model = limitOxygenDelivery(model, internalO2);
[referenceATP, referenceSolution] = runFullModel(model, settings);

%Store parameters
parameterList = zeros(14,1);
parameterList(1) = maintainance;
parameterList(2) = internalWork;
parameterList(3) = dwMuscle;
parameterList(4) = muscleRatio;
parameterList(5) = vO2perDryweight;
parameterList(6) = m2Efficency;
parameterList(7) = complex1Ratio;
parameterList(8) = vO2max;
parameterList(9) = peripheralFA;
parameterList(10) = peripheralLactateCapacity;
parameterList(11) = FAFactor;
parameterList(12) = type1tresh;
parameterList(13) = type2target;
parameterList(14) = internalO2;

%Result vectors
ATPResuls = zeros(nrOfSamples,settings.timeSteps); 
vO2Resuls = zeros(nrOfSamples,settings.timeSteps); 
vCO2Resuls = zeros(nrOfSamples,settings.timeSteps);
vglycogenResuls = zeros(nrOfSamples,settings.timeSteps);


for i = 1:nrOfSamples
    i
    randValues = (1-pertubationLevel) + 2 * pertubationLevel * rand(length(parameterList),1);
    perturb = parameterList .* randValues;
    model = setupSimulation(superModel, perturb(1), perturb(2), perturb(3),  perturb(4), perturb(5), perturb(6), perturb(7), perturb(8), perturb(9), perturb(10), perturb(11), perturb(12), perturb(13));
    model = limitOxygenDelivery(model, perturb(14));
    [ATPrate, fullSolution] = runFullModel(model, settings);
    ATPResuls(i,:) = ATPrate;
    vO2Resuls(i,:) = fullSolution(:, readouts(1));
    vCO2Resuls(i,:) = fullSolution(:, readouts(2));
    vglycogenResuls(i,:) = fullSolution(:, readouts(3));    
end

%Remove the oxygen tradeoff psuedo metabolite after simulation since it 
%interferes with data retrival
model.S(findIndex(model.mets, 'oxygenTradeof'),:) = 0;

%save('randomSampling')

%%
%Plot sampling result
subplot(3,2,1)
hold all
plot(ATPResuls', -vO2Resuls', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), -mean(vO2Resuls), 'r', 'linewidth', 2)
plot(referenceATP, -referenceSolution(:, readouts(1)), 'k', 'linewidth', 3)

title('vO2')
xlabel('mmol ATP')


subplot(3,2,2)
hold all
plot(ATPResuls', vCO2Resuls', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), mean(vCO2Resuls), 'r', 'linewidth', 2)
plot(referenceATP, referenceSolution(:, readouts(2)), 'k', 'linewidth', 3)
title('vCO2')
xlabel('mmol ATP')

subplot(3,2,3)
hold all
plot(ATPResuls', -vglycogenResuls', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), -mean(vglycogenResuls), 'r', 'linewidth', 2)
plot(referenceATP, -referenceSolution(:, readouts(3)), 'k', 'linewidth', 3)
title('vGlycogen')
xlabel('mmol ATP')

subplot(3,2,4)
hold all
plot(ATPResuls', -(vCO2Resuls./vO2Resuls)', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), -mean(vCO2Resuls./vO2Resuls), 'r', 'linewidth', 2)
plot(referenceATP, -referenceSolution(:, readouts(2))./referenceSolution(:, readouts(1)), 'k', 'linewidth', 3)
title('RQ')
xlabel('mmol ATP')

subplot(3,2,5)
hold all
dO2 = molToMl(-(vO2Resuls(:,end-4)-vO2Resuls(:,1)));
dW = molToW(1000*(ATPResuls(:,end-4)-ATPResuls(:,1)));
dO2dw1 = dO2./dW
histogram(dO2dw1, 20)
plot(median(dO2dw1) * [1 1], [0 30], 'r')
text(median(dO2dw1), 20, sprintf('Median %2.1f', median(dO2dw1)));
xlim([9 21])
title('0-55% wmax')
ylabel('nr of subjects')
xlabel('dvO2/dW')

subplot(3,2,6)
hold all
dO2 = molToMl(-(vO2Resuls(:,end)-vO2Resuls(:,end-4)));
dW = molToW(1000*(ATPResuls(:,end)-ATPResuls(:,end-4)));
dO2dw2 = dO2./dW
histogram(dO2dw2, 20)
plot(median(dO2dw2) * [1 1], [0 30], 'r')
text(median(dO2dw2), 20, sprintf('Median %2.1f', median(dO2dw2)));
xlim([9 21])
title('55%-100% wmax')
ylabel('nr of subjects')
xlabel('dvO2/dW')
