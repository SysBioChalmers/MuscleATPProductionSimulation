addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;

dwMuscle = 20*(1-0.792); %kg muscle
muscleRatio = 0.55; %type1
vO2perDryweight = 1.38 * 2.6; %Specific activity Complex IV
complex1Ratio = 33.45/76.41;
m2Efficency = 0.5;
vO2max = 1000;
internalO2 = 10;
maintainance = 5.5;
internalWork = 2;
peripheralFA = 0.004;
peripheralLactateCapacity = 2;
FAFactor = 0.8*16.55/76.41;
type1tresh = 9; %type 1 is activated first, 10 mol ATP corresponds to 20% of Wmax
type2target = 0.55; %type 2 dominates at the highest W rates

model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);
settings = setSimulationSettings(20, false);

% We can add a lactate accumulation flux to simulate non steady state
% conditions, e.g. at a maximum of 2 mol/hour 
% lactateFlux=getBounds(model, {'L-lactate[s]'});
% model.ub(lactateFlux) = 1000;
% internalLactate = getTransport(model, {'L-lactate'}, 'sb', 's');
% model.ub(internalLactate) = 2;

%oxygen delivery to muscle may be limiting
model = limitOxygenDelivery(model, internalO2);

[ATPrate, fullSolution] = runFullModel(model, settings);

%Remove the oxygen tradeoff psuedo metabolite after simulation since it 
%interferes with data retrival
model.S(findIndex(model.mets, 'oxygenTradeof'),:) = 0;

%%
close all
clf
compareWithSampleData(model, 'mikael', ATPrate, fullSolution)
figure()
plotMetaboliteList = {'O2', 'glycogen', 'L-lactate', 'CO2', 'palmitate', 'stearate'};
plotFullSolutionInternalBlood(model, ATPrate, fullSolution, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
complexIbypassM1 = getComplexIbypass(model, fullSolution, 'HMR_6921_m1');
complexIbypassM2 = getComplexIbypass(model, fullSolution, 'HMR_6921_m2');

subplot(1,3,1)
hold all
plot(1000*molToW(ATPrate(complexIbypassM1-1)) * [1 1], 9*[0 1],'k')
subplot(1,3,2)
hold all
plot(1000*molToW(ATPrate(complexIbypassM2-1)) * [1 1], 9*[0 1],'k')
figure()
simData = plotFullSolutionAandB(model, ATPrate, fullSolution, 'sb', {'sm1', 'sm2'});

figure();
plotFiberDistribution(model, ATPrate, fullSolution, muscleRatio)