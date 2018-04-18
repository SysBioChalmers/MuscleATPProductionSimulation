load('model/reducedModel')
addpath('src')

%rename subsystems
%Pyruvate to Ac-Coa towards TCA
model.subSystems{findIndex(model.rxns, 'HMR_4137')} = 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism';
model.subSystems(ismember(model.subSystems, 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism')) = {'Tricarboxylic acid cycle'};
model.subSystems(ismember(model.subSystems, 'Alanine, aspartate and glutamate metabolism')) = {'Tricarboxylic acid cycle'};



model = mapDataToRxns(model, 'data/RxnAndSA.txt');

%Set specific activity to high
model.specificActivity(findIndex(model.rxns, 'HMR_6921')) = 3;

rxnsWithValues = model.specificActivity>0;
histogram(log10(model.specificActivity(rxnsWithValues)), 6)

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
% 
% figure()
% %results = runChemostatExperimentYields(model, growthRates, minimalMedia);
% fullSolution2 = searchForNOMA(model, growthRates, fullSolution);
% 
% plotFullSolution(model, growthRates, fullSolution2, objectiveFunction);
% 
% %printFluxesValues(model, fullSolution)

%%
massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightValue = model.b(massConstraintRow,2);
weightRow = 1000 * full(model.S(massConstraintRow,:));
fluxDistribution = fullSolution(99,:);

%Normalize per half glucose
fluxDistribution = -0.5*fluxDistribution/fluxDistribution(reactionNumbers(2));
fluxAndMass = fluxDistribution .*weightRow;
equations = constructEquations(model);
[subSystem, indx] = sort(model.subSystems);

for i = 1:length(subSystem)
   curRxn =  indx(i);
   rxnName = model.rxns{curRxn};
   eqnName = equations{curRxn};
   if fluxAndMass(curRxn)>10^-5
       fprintf('%s\t%s\t%s\t%2.2f\t%2.2f\n', rxnName, eqnName, subSystem{i}, fluxDistribution(curRxn), fluxAndMass(curRxn))
   end
end



%%
fprintf('\n')
allSubs = unique(subSystem);
for i = 1:length(allSubs)
    includedRxns = ismember(model.subSystems, allSubs{i});
    sumOfMass = sum(fluxAndMass(includedRxns));
    if sumOfMass>10^-5
       fprintf('%s\t%2.2f\n', allSubs{i},  sumOfMass)
    end
end


