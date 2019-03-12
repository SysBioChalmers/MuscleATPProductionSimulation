load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;

tryWOCmplx1 = true;
subject = 'subject1';


%Parameterset (initial guess)
param = [];
param.dwMuscle = 20*(1-0.792); %kg muscle
param.muscleRatio = 0.55; %type1
param.m2Efficency = 0.5;
param.type1tresh = 9; %type 1 is activated first, 10 mol ATP corresponds to 20% of Wmax
param.type2target = 0.55; %type 2 dominates at the highest W rates
param.vO2max = 12;
param.maintainance = 5.5;
param.internalWork = 2;
param.peripheralFA = 0;
param.peripheralLactateCapacity = 0;
param.lactateBuffering = 1.3; %Lactate accumulation flux for non steady state
param.scalingFactor = 2.9;
param.HMR_6914 = 1.25; %Specific activity Complex IV
param.HMR_6921 = 0.88; %Specific activity Complex I
param.HMR_6911 = 0.07; %Specific activity ETF

fitList = {'scalingFactor', 'maintainance', 'HMR_6911'};



%Before
model = setupSimulation(model, param);
settings = setSimulationSettings(30, false);
[ATPrateBefore, fullSolutionBefore] = runFullModel(model, settings);

%After
fitedParam = fitModel(model, param, fitList, subject);
saveParam(subject, fitedParam);

model = setupSimulation(model, fitedParam);
settings = setSimulationSettings(30, false);
[ATPrateAfter, fullSolutionAfter] = runFullModel(model, settings);


if tryWOCmplx1
    %No complex 1 bypass
    fitedParamNocmplx1 = param;
    fitedParamNocmplx1.HMR_6921 = 1000;
    fitedParamNocmplx1 = fitModel(model, fitedParamNocmplx1, fitList, subject);
    model = setupSimulation(model, fitedParamNocmplx1);
    settings = setSimulationSettings(30, false);
    [ATPrateNoCmplx1, fullSolutionNoCmplx1] = runFullModel(model, settings);
end

%%
close all
clf
compareWithSampleData(model, subject, ATPrateBefore, fullSolutionBefore)
figure()
compareWithSampleData(model, subject, ATPrateAfter, fullSolutionAfter)

if tryWOCmplx1
    figure()
    compareWithSampleData(model, subject, ATPrateNoCmplx1, fullSolutionNoCmplx1)
end

if tryWOCmplx1 == false
    fitedParamNocmplx1 = fitedParam;
end

for i = 1:length(fitList)
    fprintf('%s\t%2.2f\t%2.2f\t%2.2f\n', fitList{i}, param.(fitList{i}), fitedParam.(fitList{i}), fitedParamNocmplx1.(fitList{i}))
end


