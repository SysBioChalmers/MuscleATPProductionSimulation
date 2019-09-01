close all
hold all
load('model/reducedModel')
addpath('src1')
saturation = 0.5;

model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = mapProteomToRxns(model, 'data/RxnAndProtein.txt');

%We do not trust proteomics for 6-PHOSPHOFRUCTOKINASE
model.proteinMass(findIndex(model.rxns,  'HMR_4379')) = -1;

[model, constraindRxns] = addProteinConstrant(model, saturation);

minimalMedia = {
    'O2[s]'
    'octanoic acid[s]'    
    'glycogen[s]'
    'H2O[s]'
    };

minimalFlux= [-1000
              -1000
              -1000
              -1000];


reactionNumbers = getBounds(model, minimalMedia);

model = setParam(model, 'lb', reactionNumbers, minimalFlux);

objectiveFunction = {'human_ATPMaintainance'};
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);
solution = solveLinMin(model);
maxATP = -solution.f;
atpRates = linspace(0, maxATP, 200);

model = setParam(model, 'obj', reactionNumbers(3), 1);
fullSolution = runChemostatExperiment(model, atpRates, objectiveFunction);

uBounds = model.ub(constraindRxns);
uFlux = max(abs(fullSolution(:,constraindRxns)));
eqns = constructEquations(model, constraindRxns);
enzymeUsage = uFlux'./uBounds;
rxnNames = model.rxns(constraindRxns);


%%
hold all
[enzymeUsage, indx] = sort(enzymeUsage, 'ascend');
eqns = eqns(indx);

plot(enzymeUsage, 1:length(eqns), 'o')
yticks(1:length(eqns))
yticklabels(eqns)

vO2max = solution.x(reactionNumbers(1));

grid on

ylim([0 length(eqns)]+0.5)

xlim([10^-4 1])
set(gca, 'xscale', 'log')
xlabel('enzyme usage')
%%
respirationChain = {
'HMR_6921'
'HMR_4652'
'HMR_6918'
'HMR_6914'
'HMR_6916'};
name = {'I', 'II', 'III', 'IV', 'V'};
solutions = find(sum(abs(fullSolution),2)>0);
solutions = solutions(end);
fluxDist = fullSolution(solutions,:);

for i = 1:length(respirationChain)
    curRxn = findIndex(model.rxns, respirationChain{i});
    proteinMass = model.proteinMass(curRxn);
    specificActivity = model.specificActivity(curRxn);
    curFlux = fluxDist(curRxn);
    
    if curFlux < 0
        curFlux = -curFlux;
    end
    
    predictedMass = curFlux/(60*specificActivity*0.5);  
        
    if proteinMass>0
       fprintf('%s\t%2.6f\t%2.6f\n', name{i}, predictedMass, proteinMass)
    end
end

