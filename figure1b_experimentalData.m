close all
hold all
load('model/reducedModel')
addpath('src1')
saturation = 0.5;


model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = mapProteomToRxns(model, 'data/RxnAndProtein.txt');

%We do not trust proteomics for 6-PHOSPHOFRUCTOKINASE
model.proteinMass(findIndex(model.rxns,  'HMR_4379')) = -1;

%We do not trust proteomics for TPI
model.proteinMass(findIndex(model.rxns,  'HMR_4391')) = -1;

%We do not trust proteomics for GPD2
model.proteinMass(findIndex(model.rxns,  'HMR_0483')) = -1;


[model, constraindRxns] = addProteinConstrant(model, saturation);



minimalMedia = {
    'O2[s]'
    'glycogen[s]'
    'H2O[s]'
    };

plotExchange = {
    'O2[s]'
    'glycogen[s]'
    'L-lactate[s]'
    };

minimalFlux= [-1000
              -1000
              -1000];

growthRates = linspace(0, 8.5, 100);

reactionNumbers = getBounds(model, minimalMedia);

objectiveFunction = {'human_ATPMaintainance'};

model = setParam(model, 'obj', reactionNumbers(2), 1);

fullSolution = runChemostatExperiment(model, growthRates, objectiveFunction);

plotFullSolution(model, growthRates, fullSolution, plotExchange);
ylim([0, 2])
figure()

model = configureModel(model, minimalMedia, minimalFlux);

model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);

solution = solveLin(model);




%plotFullSolution(model, growthRates, fullSolution, plotExchange);
hold all


uBounds = model.ub(constraindRxns);
uFlux = abs(solution.x(constraindRxns));
eqns = constructEquations(model, constraindRxns)
enzymeUsage = saturation * uFlux./uBounds;
[enzymeUsage, indx] = sort(enzymeUsage);
eqns = eqns(indx);

plot(enzymeUsage, 1:length(eqns), 'o-')
yticks(1:length(eqns))
yticklabels(eqns)

vO2max = solution.x(reactionNumbers(1));

%Simulate max flux through complex 1 and 2
model = setParam(model, 'lb', reactionNumbers(2), 0);

%Simulate Complex 1
glycogenPhos = createRXNStuct(model, 'Complex1Flux', 'NAD+[m] => NADH[m] + H+[m]', 0, 1000, 'Simulate complex I');
model=addRxns(model,glycogenPhos,3,'m',false);

model.ub(findIndex(model.rxns, 'Complex1Flux')) = model.ub(findIndex(model.rxns, 'HMR_6921'));
solution = solveLin(model);
Complex1O2 = -solution.x(reactionNumbers(1));

%Simulate Complex 1 + 2
glycogenPhos = createRXNStuct(model, 'Complex2Flux', 'fumarate[m] => succinate[m]', 0, 1000, 'Simulate complex II');
model=addRxns(model,glycogenPhos,3,'m',false);
model.ub(findIndex(model.rxns, 'Complex1Flux')) = model.ub(findIndex(model.rxns, 'HMR_4652'));
solution = solveLin(model);
Complex2O2 = -solution.x(reactionNumbers(1));


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



%%
figure()




factor = 1/(1-0.792) * 60*60 * 10^-12 * 10^3 * 10^3; %per hour per gdw 

data = factor*[
    33.45	2.11
    76.41	4.58
%    105.99	5.63
    ];

exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;

hold all     
bar([1 3], data(:,1), 0.5, 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
errorbar([1 3], data(:,1), data(:,2),'k.')
bar([2 4], [Complex1O2 Complex2O2], 0.5, 'FaceColor', exMap(2,:), 'EdgeColor', 'none')

set(gca, 'XTick', [1 2 3 4])
set(gca, 'XTickLabel', {'Complex I', '(model)', 'Complex I+II', '(model)'})
ylabel('vO2 [mmol/gdw/h]')
ylim([0, 2])
xtickangle(45)
set(findall(gcf,'-property','FontSize'),'FontSize',15)