load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;
settings = [];

dwMuscle = 20*(1-0.792); %kg muscle
m1Ratio = 0.55;
trainingEffect = 2.2;

c1 = 1.38 * 1.2 * trainingEffect; %mmol O2/gdw (compensating for fiber difference and training
c2 = c1*0.51; %mmol O2/gdw

model = addInternalConstraints(model, dwMuscle, m1Ratio, c1, c2);


%Exchange fluxes:
model = addTransportConstraints(model,   's',   {'O2'}, [-11.9], [0]);
model = addTransportConstraints(model, 'sm1', {'palmitate', 'L-lactate'}, [-1000 -1000], [0 1000]);
model = addTransportConstraints(model, 'sm2', {'palmitate', 'L-lactate'}, [-1000 0], [0 1000]);
model = addTransportConstraints(model, 'sm3', {'palmitate', 'L-lactate'}, [-0.025 -1000],  [1000 0]);

%Special constraints
model = addSpecializedConstraints(model, 1000, 4.3+2.4, 1000);

settings = addExchangeMedum(settings);
settings.timeSteps = 20;

%Make Uncoupling favorable over ATP degradation
model = configureSMatrix(model, 10, 'HMR_7638_m3', 'H+[cm3]');
model = configureSMatrix(model, -10, 'HMR_7638_m3', 'H+[mm3]');

%quad, lincomb, uncoupling, lactate
[fullSolution, ATPrate] = optimizeFluxes(model, settings, 'AandB');
%[fullSolution2, ATPrate] = optimizeFluxes(model, settings, 'AandB');
%AandB = mixSolutions(fullSolution1, fullSolution2, 0.4, 0.55);
%fullSolutionMod = crossAndScale(fullSolution, ATPrate, 23, 36, 1);

fullSolutionMod=fullSolution;
close all
clf
plotMetaboliteList = {'O2', 'glycogen', 'L-lactate', 'CO2', 'palmitate'};
plotFullSolutionInternalBlood(model, ATPrate, fullSolutionMod, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
figure()
compareWithSampleData(model, 'mikael', ATPrate, fullSolutionMod, 1)
figure()
simData = plotFullSolutionAandB(model, ATPrate, fullSolutionMod, 'sb', {'sm1', 'sm2'});
%figure()
%compareWithKneeData(simData, 'BI', 10)
figure()
ucp3 = findIndex(model.rxns, 'HMR_7638_m3');
plot(ATPrate, fullSolution(:,ucp3));


