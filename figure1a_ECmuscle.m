load('model/reducedModel')
addpath('src')

model = mapDataToRxns(model, 'data/RxnAndSA.txt');

%Set specific activity to high
model.specificActivity(findIndex(model.rxns, 'HMR_6921')) = 3;

model = addSpecificActivityConstraint(model, 0.35, 0.1*0.5, 60);
model = addReversedReactions(model);

allowUncoupling = true;
complex1Bypass = false;
removeLactate = false;

close all

minimalMedia = {
    'O2[s]'
    'glycogen[s]'
    'H2O[s]'
    };
minimalFlux= [-1000
              -1000
              -1000];

growthRates = linspace(0, 9, 100);

reactionNumbers = getBounds(model, minimalMedia);
model = setParam(model, 'obj', reactionNumbers(2), 1);

objectiveFunction = {'human_ATPMaintainance'};


model = configureModel(model, minimalMedia, minimalFlux);

%Remove lactate
if removeLactate
     model = setParam(model, 'ub', 'HMR_9135', 0); 
end

%bypass of complex 1
if complex1Bypass == true
    model = setParam(model, 'ub', 'HMR_6921', 0); 
end

%block native uncoupling
if allowUncoupling
    model = setParam(model, 'ub', 'HMR_7638', 1000); 
else
    model = setParam(model, 'ub', 'HMR_7638', 0); %Direct uncoupling
end

%model = setParam(model, 'ub', 'HMR_6916', 0); 

fullSolution = runChemostatExperiment(model, growthRates, objectiveFunction);

plotExchange = {
    'O2[s]'
    'glycogen[s]'
    'L-lactate[s]'
    };

plotFullSolution(model, growthRates, fullSolution, plotExchange);

figure()
%results = runChemostatExperimentYields(model, growthRates, minimalMedia);
fullSolution2 = searchForNOMA(model, growthRates, fullSolution);

plotFullSolution(model, growthRates, fullSolution2, objectiveFunction);

%printFluxesValues(model, fullSolution)

