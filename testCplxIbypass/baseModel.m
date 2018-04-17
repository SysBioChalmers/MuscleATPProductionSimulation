%load('reducedModel')
load('../model/muscleModel')
addpath('../src')

%model = changeBiomassEquation(model, {'human_growthMaintainance[c]'}, []);
objectiveFunction = {'human_ATPMaintainance'};
ATPProduction = findIndex(model.rxns, objectiveFunction);


%model.lb(findIndex(model.rxns, 'HMR_5998')) = -1000; %open lactate export


glycogenPhos = createRXNStuct(model, 'GlycogenPhosphorylase', 'glycogen[c] + Pi[c] => glucose-1-phosphate[c]', 0, 1000, 'Starch and sucrose metabolism');
model=addRxns(model,glycogenPhos,3,'c',false);

%Minimal media simulation
inMedia = {
    'glycogen[s]'
    'O2[s]'
    'H+[s]'    
};

outMedia = {
    'HCO3-[s]'
    'CO2[s]' 
    'H+[s]'
    'H2O[s]'
};


[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');


%Remove ProtPool

reactionIn = getBounds(model, inMedia);
reactionOut = getBounds(model, outMedia);

model = setParam(model, 'lb', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', reactionOut, 1000);
model = setParam(model, 'lb', reactionIn, -1000);

model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);

model = setParam(model, 'lb', reactionIn(1), -1);

model = setParam(model, 'ub', reactionOut(1), 0);

model = setParam(model, 'obj', objectiveFunction, 1);

solutionA = solveLinMin(model, true)

%Complex 1 knock out
model = setParam(model, 'lb', 'HMR_6921', 0);
model = setParam(model, 'ub', 'HMR_6921', 0);

solutionB = solveLinMin(model, true)
%printAffected(model, solutionA, solutionB);


%glycerol-3-phosphate shunt knock
model = setParam(model, 'lb', 'HMR_0483', 0);
model = setParam(model, 'ub', 'HMR_0483', 0);
solutionC = solveLinMin(model, true)
%printAffected(model, solutionB, solutionC)


%ETF-ubiquinone oxidoreductase knock
model = setParam(model, 'lb', 'HMR_6911', 0);
model = setParam(model, 'ub', 'HMR_6911', 0);

solutionD = solveLinMin(model, true)
printAffected(model, solutionC, solutionD)



