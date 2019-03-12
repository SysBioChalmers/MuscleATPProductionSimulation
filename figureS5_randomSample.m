load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;

%Random sampling settings
nrOfSamples = 100;
nrOfTimesteps = 10;
pertubationLevel = 0.2; %+- 20%

subjectId = 'subject1';

param = loadParm(subjectId, true);

%Make reference condition
model = setupSimulation(model, param);
settings = setSimulationSettings(nrOfTimesteps, false);
[referenceATP, referenceSolution] = runFullModel(model, settings);

%Reactions to plot
readouts = getTransport(model, {'O2', 'CO2', 'glycogen'}, 'sb', 's');

%Result vectors
ATPResuls = zeros(nrOfSamples,settings.timeSteps); 
vO2Resuls = zeros(nrOfSamples,settings.timeSteps); 
vCO2Resuls = zeros(nrOfSamples,settings.timeSteps);
vglycogenResuls = zeros(nrOfSamples,settings.timeSteps);

%Store reference parameters
paramRef = param;
fns = fieldnames(param);

for i = 1:nrOfSamples
    i
    randValues = (1-pertubationLevel) + 2 * pertubationLevel * rand(length(fns),1);
    param = paramRef;
    for j = 1:length(randValues)
        param.(fns{j}) = param.(fns{j}) .* randValues(j);
    end
    
    model = setupSimulation(superModel, param);
    [ATPrate, fullSolution] = runFullModel(model, settings);
    ATPResuls(i,:) = ATPrate;
    vO2Resuls(i,:) = fullSolution(:, readouts(1));
    vCO2Resuls(i,:) = fullSolution(:, readouts(2));
    vglycogenResuls(i,:) = fullSolution(:, readouts(3));    
end

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
dO2dw1 = dO2./dW;
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
