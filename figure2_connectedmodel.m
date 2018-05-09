addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;
settings = [];

%Make Uncoupling favorable over ATP degradation
model = configureSMatrix(model, 10, 'HMR_7638_m3', 'H+[cm3]');
model = configureSMatrix(model, -10, 'HMR_7638_m3', 'H+[mm3]');

dwMuscle = 20*(1-0.792); %kg muscle
m1Ratio = 0.55;
vO2perDryweight = 1.38 * 1.2 * 2.4;
complex1Ratio = 33.45/76.41;
m2Efficency = 0.5;
vO2max = 11.9;
maintainance = 4.3;
internalWork = 2.4;
peripheralFA = 0.01;
peripheralFAsynth = 0.2;
FAFactor = 16.55/76.41 * 0.7;

settings.timeSteps = 20;
settings.pfba = false;

[ATPrate, fullSolution] = setupAndRunSimulation(model, settings, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralFAsynth, FAFactor);

%%
close all
clf
compareWithSampleData(model, 'mikael', ATPrate, fullSolution, 1)
figure()
plotMetaboliteList = {'O2', 'glycogen', 'L-lactate', 'CO2', 'palmitate'};
plotFullSolutionInternalBlood(model, ATPrate, fullSolution, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
figure()
simData = plotFullSolutionAandB(model, ATPrate, fullSolution, 'sb', {'sm1', 'sm2'});
figure()
ucp3 = findIndex(model.rxns, 'HMR_7638_m3');
plot(ATPrate, fullSolution(:,ucp3));


