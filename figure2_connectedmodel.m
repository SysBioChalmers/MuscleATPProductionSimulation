addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;
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

settings.timeSteps = 20;
settings.pfba = false;

[ATPrate, fullSolution] = setupAndRunSimulation(model, settings, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);

%%
close all
clf
compareWithSampleData(model, 'mikael', ATPrate, fullSolution, 1)
figure()
plotMetaboliteList = {'O2', 'glycogen', 'L-lactate', 'CO2', 'palmitate', 'stearate'};
plotFullSolutionInternalBlood(model, ATPrate, fullSolution, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
figure()
simData = plotFullSolutionAandB(model, ATPrate, fullSolution, 'sb', {'sm1', 'sm2'});


