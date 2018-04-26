load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;
settings = [];

mScale = 37*(1-0.8); %kg muscle
m1Ratio = 0.55;
c1 = 2; %mmol O2/gdw
c2 = c1/2; %mmol O2/gdw

model = addInternalConstraints(model, mScale, m1Ratio, c1, c2);

%Exchange fluxes:
model = addTransportConstraints(model,   's',   {'O2'}, [-12], [0]);
model = addTransportConstraints(model, 'sm1', {'palmitate', 'O2', 'L-lactate'}, [-0.03 -1000 -1000 0], [0 0 1000]);
model = addTransportConstraints(model, 'sm2', {'palmitate', 'O2' , 'L-lactate'}, [-0.03 -1000 0], [0 0 1000]);
model = addTransportConstraints(model, 'sm3', {'palmitate', 'L-lactate'}, [-0.03 -1000],  [0.18 0]);
%Special constraints
model = addSpecializedConstraints(model, 1000, 7.5, 1000);

settings = addExchangeMedum(settings);
settings.timeSteps = 20;

%quad, lincomb, uncoupling, lactate
[fullSolution, ATPrate] = optimizeFluxes(model, settings, 'AandB');
%[fullSolution2, ATPrate] = optimizeFluxes(model, settings, 'AandB');
%AandB = mixSolutions(fullSolution1, fullSolution2, 0.4, 0.55);
%fullSolutionMod = crossAndScale(fullSolution, ATPrate, 23, 36, 1);

fullSolutionMod=fullSolution;
close all
clf
plotFullSolutionInternalBlood(model, ATPrate, fullSolutionMod, 'sb', {'sm1', 'sm2', 'sm3'});
figure()
compareWithSampleData(model, 'mikael', ATPrate, fullSolutionMod, 1)
figure()
simData = plotFullSolutionAandB(model, ATPrate, fullSolutionMod, 'sb', {'sm1', 'sm2'});
%figure()
%compareWithKneeData(simData, 'BI', 10)

