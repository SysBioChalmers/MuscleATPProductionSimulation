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
complexIrxn = findIndex(model.rxns, 'HMR_6921') == find(constraindRxns);
fprintf('CI flux: %1.3f, bound: %1.3f\n', uFlux(complexIrxn), uBounds(complexIrxn))
GPDreaction = findIndex(model.rxns, 'HMR_0479') == find(constraindRxns);
fprintf('GPD flux: %1.3f, bound: %1.3f\n', uFlux(GPDreaction), uBounds(GPDreaction))

%%
hold all
xdata = log10(uBounds);
ydata = log10(uFlux');
xvals = linspace(min(xdata), max(xdata));
%plot([-2, 3], [-2, 3], 'color', [0.5 0.5 0.5], 'linewidth', 2, 'HandleVisibility','off');
fill([-2, 3 -2], [-2, 3 3], [0.85 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility','off');


fillColors = [106 189 69
              237 34 36
              58 83 164
              100 100 100
            ]/256;

subSystemCategories = {'Glycolysis / Gluconeogenesis', 'Tricarboxylic acid cycle', 'Oxidative phosphorylation', 'Other'};        
subSystemShortnames = {'Glycolysis', 'TCA', 'OXPHOS', 'Other'};    
allSubs = model.subSystems(constraindRxns);
allSubs(not(ismember(allSubs,subSystemCategories))) = {'Other'}; 

fprintf('------\n')
for i = 1:length(xdata)
   curSub = find(contains(subSystemCategories, allSubs{i}));
   scatter(xdata(i), ydata(i), 'filled', 'MarkerFaceColor', fillColors(curSub,:), 'MarkerEdgeColor', 'none', 'HandleVisibility','off')

   if abs(ydata(i)-xdata(i))<0.01
       text(xdata(i), ydata(i), eqns{i});
       fprintf('%s\t%2.2f\t(%s)\n', rxnNames{i}, uBounds(i), eqns{i})
   end    
end
fprintf('------\n')

notOther = not(ismember(allSubs,'Other'));
[c1, p1] = corr(xdata(notOther), ydata(notOther));
[c2, p2] = corr(xdata, ydata);
[p1 p2]

text(1.5, 0, sprintf('r=%2.2f\np=%2.1e', c1, p1))
linearModel = fitlm(xdata(notOther), ydata(notOther));
km = linearModel.Coefficients.Estimate;
plot(xvals, km(1) + km(2)*xvals, 'k', 'linewidth', 2);


axis equal
xlim([-2 3])
ylim([-2 3])
xlabel('log10(Vmax)')
ylabel('log10(Flux)')

for i = 1:length(subSystemShortnames)
    h(i) = scatter(-5,-5, 'MarkerFaceColor', fillColors(i,:), 'MarkerEdgeColor', 'none');
end
legend(h, subSystemShortnames, 'location', 'NW')
legend boxoff



