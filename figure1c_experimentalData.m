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
%model.proteinMass(findIndex(model.rxns,  'HMR_0483')) = -1;


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


reactionNumbers = getBounds(model, minimalMedia);

model = setParam(model, 'lb', reactionNumbers, minimalFlux);

objectiveFunction = {'human_ATPMaintainance'};
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);
solution = solveLinMin(model);
maxATP = -solution.f;
atpRates = linspace(0, maxATP, 200);

model = setParam(model, 'obj', reactionNumbers(2), 1);
fullSolution = runChemostatExperiment(model, atpRates, objectiveFunction);

plotFullSolution(model, atpRates, fullSolution, plotExchange);


oxygenFlux = fullSolution(:, findIndex(model.rxns, 'HMR_9048'));
DeltaO2vsDeltaATP = (diff(oxygenFlux)./diff(atpRates'));
DeltaO2vsDeltaATP(20)/DeltaO2vsDeltaATP(2)
breakPoints = find(abs(diff(DeltaO2vsDeltaATP))>0.001);

hold all
plot(atpRates(breakPoints(1)) * [0 1], -oxygenFlux(breakPoints(1))*[1 1], 'k--')
plot(atpRates(breakPoints(1)) * [1 1], -oxygenFlux(breakPoints(1))*[0 1], 'k--')

plot(atpRates(breakPoints(end)) * [0 1], -oxygenFlux(breakPoints(end))*[1 1], 'k--')
plot(atpRates(breakPoints(end)) * [1 1], -oxygenFlux(breakPoints(end))*[0 1], 'k--')

xlim([0 8])
ylim([0 2])

%%
ylim([0, 2])
figure()

%model = configureModel(model, minimalMedia, minimalFlux);


%plotFullSolution(model, growthRates, fullSolution, plotExchange);

hold all

uBounds = model.ub(constraindRxns);
uFlux = max(abs(fullSolution(:,constraindRxns)));
eqns = constructEquations(model, constraindRxns);
enzymeUsage = saturation * uFlux'./uBounds;
[enzymeUsage, indx] = sort(enzymeUsage);
eqns = eqns(indx);

plot(enzymeUsage, 1:length(eqns), 'o-')
yticks(1:length(eqns))
yticklabels(eqns)

vO2max = solution.x(reactionNumbers(1));


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
    16.55   1.41
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
bar([1 2 3], data(:,1), 0.8, 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
errorbar([1 2 3], data(:,1), data(:,2),'k.')
%bar([2 4], [Complex1O2 Complex2O2], 0.5, 'FaceColor', exMap(2,:), 'EdgeColor', 'none')

set(gca, 'XTick', [1 2 3 4])
set(gca, 'XTickLabel', {'Fat', 'Complex I', 'Complex I+II'})
ylabel('vO2 [mmol/gdw/h]')
ylim([0, 2])
xlim([0.5 3.5])
xtickangle(45)
set(findall(gcf,'-property','FontSize'),'FontSize',15)