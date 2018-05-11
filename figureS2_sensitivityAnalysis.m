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
trainingEffect = 2.25;
vO2perDryweight = 1.38 * 1.2 * 2.25;
m2Efficency = 0.5;
vO2max = 11.9;
maintainance = 4.3;
internalWork = 2.4;
peripheralFA = 0.02;
peripheralLactateCapacity = 1.8;
FAFactor = 16.55/76.41 * 0.6;

settings.timeSteps = 3;
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
results = zeros(length(parameterList),4);

%Reactions to plot
readouts = getTransport(model, {'O2', 'CO2', 'glycogen'}, 'sb', 's');

%Reference result
resultReference = zeros(1,4);
resultReference(1) = referenceATP(end);
resultReference(2:4) = referenceSolution(end, readouts);
resultReference(5) = resultReference(3)/-resultReference(2);

for i = 1:length(parameterList)
    perturb = parameterList;
    perturb(i) = perturb(i)*pertubationFactor;
    [ATPrate, fullSolution] = setupAndRunSimulation(model, settings, perturb(1), perturb(2), perturb(3),  perturb(4), perturb(5), perturb(6), perturb(7), perturb(8), perturb(9), perturb(10), perturb(11));
    results(i,1) = ATPrate(end);
    results(i,2:4) = fullSolution(end, readouts);  
    results(i,5) = results(i,3)/-results(i,2);
end

%%
%Plot sampling result
fprintf('parameterName\tWatt\tvO2\tvCo2\tvglycogen\tRQ\n')

for i = 1:length(parameterList)
    percentCange = resultReference./results(i,:);
    E = -(1-percentCange)./(1-pertubationFactor);
    fprintf('%s\t%2.5f\t%2.5f\t%2.5f\t%2.5f\t%2.5f\n', parameterNames{i}, E(1), E(2), E(3), E(4), E(5)); 
end

