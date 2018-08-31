load('model/reducedModel')
addpath('src1')

%rename subsystems
%Pyruvate to Ac-Coa towards TCA

model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = mapProteomToRxns(model, 'data/RxnAndProtein.txt');


%Set specific activity to high
rxnsWithValues = model.specificActivity>0;
histogram(log10(model.specificActivity(rxnsWithValues)), 6)

model = addSpecificActivityConstraint(model, 0.5, 0.0444*0.72, 60);
model = addReversedReactions(model);

allowUncoupling = true;
complex1Bypass = false;
removeLactate = false;

massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightValue = model.b(massConstraintRow,2);
weightRow = full(model.S(massConstraintRow,:));

close all

minimalMedia = {
    'O2[s]'
    'glycogen[s]'
    'H2O[s]'
    };
minimalFlux= [-1000
              -1000
              -1000];

growthRates = linspace(0, 22.1, 100);

reactionNumbers = getBounds(model, minimalMedia);
model = setParam(model, 'obj', reactionNumbers(2), 1);

objectiveFunction = {'human_ATPMaintainance'};
model = setParam(model, 'lb', 'human_ATPMaintainance', 0); 

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
fullMass = fullSolution .* repmat(weightRow, size(fullSolution,1),1);
massOfFlux = sum(fullMass,2);

plotExchange = {
    'O2[s]'
    'glycogen[s]'
    'L-lactate[s]'
    };
plotFullSolution(model, growthRates, fullSolution, plotExchange);
ylim([0 15])


%%

fluxDistribution = fullSolution(2,:);

%Normalize per half glucose
fluxDistribution = -0.5*fluxDistribution/fluxDistribution(reactionNumbers(2));
fluxAndMass = 1000 * fluxDistribution .* weightRow;
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

respirationChain = {
'HMR_6921'
'HMR_4652'
'HMR_6918'
'HMR_6914'
'HMR_6916'};
name = {'I', 'II', 'III', 'IV', 'V'};

for i = 1:length(respirationChain)
    curRxn = findIndex(model.rxns, respirationChain{i});
    proteinMass = model.proteinMass(curRxn);
    predictedMass = fluxAndMass(curRxn)/1000;
    
    if predictedMass == 0
        curRxn = findIndex(model.rxns, [respirationChain{i} '_back']);
        predictedMass = fluxAndMass(curRxn)/1000;        
    end
        
    if sumOfMass>10^-5
       fprintf('%s\t%2.6f\t%2.6f\n', name{i}, predictedMass, proteinMass)
    end
end

