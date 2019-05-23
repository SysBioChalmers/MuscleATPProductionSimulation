addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;

subjectId = 'subject1';

complex1Bypass = true;

param = loadParm(subjectId, true);

settings = setSimulationSettings(30, false);

model = setupSimulation(model, param);
[ATPrate, fullSolution] = runFullModel(model, settings);
%saveFluxes('fluxes/incremental.txt', model, fullSolution', ATPrate, 30);
%should be run with the parsimonius flag set as true in settings

if complex1Bypass == false
    param.HMR_6921 = 1000;
    model2 = setupSimulation(model, param);
    [ATPrateWOCIBP, fullSolutionWOCIBP] = runFullModel(model2, settings);
end


%%
close all
clf
if complex1Bypass == false
   appendCplx1bypass(model, ATPrateWOCIBP, fullSolutionWOCIBP); 
end
compareWithSampleData(model, subjectId, ATPrate, fullSolution)


figure()
plotMetaboliteList = {'O2', 'glycogen', 'L-lactate', 'CO2', 'palmitate'};
plotFullSolutionInternalBlood(model, ATPrate, fullSolution, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
complexIbypassM1 = getComplexIbypass(model, fullSolution, 'HMR_6921_m1');
complexIbypassM2 = getComplexIbypass(model, fullSolution, 'HMR_6921_m2');

subplot(1,3,1)
hold all
plot(1000*molToW(ATPrate(complexIbypassM1-1)) * [1 1], 9*[0 1],'k')
subplot(1,3,2)
hold all
plot(1000*molToW(ATPrate(complexIbypassM2-1)) * [1 1], 9*[0 1],'k')
%figure()
simData = plotFullSolutionAandB(model, ATPrate, fullSolution, 'sb', {'sm1', 'sm2'});
%figure();
%plotFiberDistribution(model, ATPrate, fullSolution, muscleRatio)
%figure();
%plotPathwayMap(model, ATPrate, fullSolution, settings)
figure();
plotSubstrateDistribution(model, ATPrate, fullSolution)

