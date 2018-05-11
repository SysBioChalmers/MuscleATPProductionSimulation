load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;


%Random sampling settings
nrOfSamples = 50;
pertubationFactor = 1.001; %.1 difference


%Parameterset
settings = [];
dwMuscle = 20*(1-0.792); %kg muscle
m1Ratio = 0.55;
vO2perDryweight = 1.38 * 1.2 * 2.5;
complex1Ratio = 33.45/76.41;
m2Efficency = 0.5;
vO2max = 11.9;
maintainance = 4;
internalWork = 2;
peripheralFA = 0.01;
peripheralLactateCapacity = 3;
FAFactor = 16.55/76.41 * 0.7;


settings.timeSteps = 10;
settings.pfba = false;

%Make reference condition
[referenceATP, referenceSolution] = setupAndRunSimulation(model, settings, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);

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
parameterNames = {'maintainance'
                  'internalWork'
                  'dwMuscle'
                  'm1Ratio'
                  'vO2perDryweight'
                  'm2Efficency'
                  'complex1Ratio'
                  'vO2max'
                  'peripheralFA'
                  'peripheralFAsynth'
                  'FAFactor'};




%Result vectors
results = zeros(length(parameterList),settings.timeSteps);

%Reactions to plot
readouts = getTransport(model, {'O2'}, 'sb', 's');

%Reference result
resultReference = referenceATP./referenceSolution(:, readouts)';

for i = 1:length(parameterList)
    perturb = parameterList;
    perturb(i) = perturb(i)*pertubationFactor;
    [ATPrate, fullSolution] = setupAndRunSimulation(model, settings, perturb(1), perturb(2), perturb(3),  perturb(4), perturb(5), perturb(6), perturb(7), perturb(8), perturb(9), perturb(10), perturb(11));
    results(i,:) = ATPrate./fullSolution(:,readouts)';
end

%%
%Plot sampling result
fprintf('parameterName')
for i = 1:length(referenceATP)
    fprintf('\t%2.2f', referenceATP(i)/max(referenceATP)); 
end
fprintf('\n'); 

for i = 1:length(parameterList)
    percentCange = resultReference./results(i,:);
    E = -(1-percentCange)./(1-pertubationFactor);
    fprintf('%s', parameterNames{i}); 
    for j = 1:length(E)
        fprintf('\t%2.5f', E(j)); 
    end
    fprintf('\n'); 
end

