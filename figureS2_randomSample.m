load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;


%Random sampling settings
nrOfSamples = 50;
pertubationLevel = 0.2; %+- 20%


%Parameterset
settings = [];
dwMuscle = 20*(1-0.792); %kg muscle
m1Ratio = 0.55;
trainingEffect = 2.25;
vO2perDryweight = 1.38 * 1.2 * 2.25;
m2Efficency = 0.5;
vO2max = 11.9;
maintainance = 4.3;
internalWork = 2.4;
peripheralFA = 0.02;
peripheralLactateCapacity = 1.8;
FAFactor = 16.55/76.41 * 0.6;



%Make reference condition
model = setupSimulation(model, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);
settings = setSimulationSettings(20, false);
[referenceATP, referenceSolution] = runFullModel(model, settings);

%Store parameters
parameterList = zeros(11,1);
parameterList(1) = maintainance;
parameterList(2) = internalWork;
parameterList(3) = dwMuscle;
parameterList(4) = m1Ratio;
parameterList(5) = vO2perDryweight;
parameterList(6) = m2Efficency;
parameterList(7) = complex1Ratio;
parameterList(8) = vO2max;
parameterList(9) = peripheralFA;
parameterList(10) = peripheralLactateCapacity;
parameterList(11) = FAFactor;

%Result vectors
ATPResuls = zeros(nrOfSamples,settings.timeSteps); 
vO2Resuls = zeros(nrOfSamples,settings.timeSteps); 
vCO2Resuls = zeros(nrOfSamples,settings.timeSteps);
vglycogenResuls = zeros(nrOfSamples,settings.timeSteps);

%Reactions to plot
readouts = getTransport(model, {'O2', 'CO2', 'glycogen'}, 'sb', 's');

for i = 1:nrOfSamples
    i
    randValues = (1-pertubationLevel) + 2 * pertubationLevel * rand(length(parameterList),1);
    perturb = parameterList .* randValues;
    modelPerturbed = setupSimulation(model, perturb(1), perturb(2), perturb(3),  perturb(4), perturb(5), perturb(6), perturb(7), perturb(8), perturb(9), perturb(10), perturb(11));
    [ATPrate, fullSolution] = runFullModel(model, settings);
    ATPResuls(i,:) = ATPrate;
    vO2Resuls(i,:) = fullSolution(:, readouts(1));
    vCO2Resuls(i,:) = fullSolution(:, readouts(2));
    vglycogenResuls(i,:) = fullSolution(:, readouts(3));    
end
%save('randomSampling')

%%
%Plot sampling result
subplot(1,4,1)
hold all
plot(ATPResuls', -vO2Resuls', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), -mean(vO2Resuls), 'r', 'linewidth', 2)
plot(referenceATP, -referenceSolution(:, readouts(1)), 'k', 'linewidth', 3)

title('vO2')
xlabel('mmol ATP')


subplot(1,4,2)
hold all
plot(ATPResuls', vCO2Resuls', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), mean(vCO2Resuls), 'r', 'linewidth', 2)
plot(referenceATP, referenceSolution(:, readouts(2)), 'k', 'linewidth', 3)
title('vCO2')
xlabel('mmol ATP')

subplot(1,4,3)
hold all
plot(ATPResuls', -vglycogenResuls', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), -mean(vglycogenResuls), 'r', 'linewidth', 2)
plot(referenceATP, -referenceSolution(:, readouts(3)), 'k', 'linewidth', 3)
title('vGlycogen')
xlabel('mmol ATP')

subplot(1,4,4)
hold all
plot(ATPResuls', -(vCO2Resuls./vO2Resuls)', 'color', 0.5*[1 1 1])
plot(mean(ATPResuls), -mean(vCO2Resuls./vO2Resuls), 'r', 'linewidth', 2)
plot(referenceATP, -referenceSolution(:, readouts(2))./referenceSolution(:, readouts(1)), 'k', 'linewidth', 3)
title('RQ')
xlabel('mmol ATP')


