load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;
settings = [];

dwMuscle = 3.6149*(1-0.792); %kg muscle
m1Ratio = 0.55;
vO2perDryweight = 1.38 * 1.2 * 2.4;
complex1Ratio = 33.45/76.41;
m2Efficency = 0.5;
vO2max = 11.9;
maintainance = 4.3;
internalWork = 2.4/2;
peripheralFA = 0.01;
peripheralLactateCapacity = 0.55;
FAFactor = 16.55/76.41 * 0.7;

settings.timeSteps = 5;
settings.pfba = true;

[ATPrate, fullSolution] = setupAndRunSimulation(model, settings, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);


%%
close all
clf
compareWithKneeData(model, ATPrate, fullSolution, 'KE', 0)
figure()
plotFullSolutionInternalBlood(model, ATPrate, fullSolution, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
% figure()
% ucp3 = findIndex(model.rxns, 'HMR_7638_m1');
% plot(ATPrate, fullSolution(:,ucp3))
% 

