load('model/reducedModel')
addpath('src1')

model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = mapProteomToRxns(model, 'data/RxnAndProtein.txt');
model = addSpecificActivityConstraint(model, 0.5, 1000, 60);
massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightRow = full(model.S(massConstraintRow,:));

close all

minimalMedia = {
    'O2[s]'
    'glycogen[s]'
    'H2O[s]'
    };

minimalFlux= [-1000
              -1
              -1000];

model = configureModel(model, minimalMedia, minimalFlux);
objectiveFunction = {'human_ATPMaintainance'};

model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);

%solutionNames = {'Aerobic', 'Complex I bypass', 'Uncoupling', 'Fermentative'};
solutionNames = {'Aerobic', 'Complex I bypass', 'Fermentative', 'Uncoupling'};
solutions = zeros(4,length(model.rxns));

%Remove lactate
res = solveLinMin(model);
solutions(1,:) = res.x;

%bypass of complex 1
model = setParam(model, 'ub', 'HMR_6921', 0); 
res = solveLinMin(model);
solutions(2,:) = res.x;

%uncoupling
model = setParam(model, 'ub', 'HMR_6916', 0); 
res = solveLinMin(model);
solutions(4,:) = res.x;

%Remove oxygen
model = setParam(model, 'ub', 'HMR_6914', 0); 
res = solveLinMin(model);
solutions(3,:) = res.x;

fullMass = abs(solutions) .* repmat(weightRow, size(solutions,1),1);
massOfFlux = sum(fullMass,2);

atpFlux = solutions(:,findIndex(model.rxns, objectiveFunction));
glyFlux = -solutions(:,findIndex(model.rxns, 'HMR_9728'));
O2flux = -solutions(:,findIndex(model.rxns, 'HMR_9048'));
TCAflux = solutions(:,findIndex(model.rxns, 'HMR_4152'));

glyYields = atpFlux./glyFlux;
O2ATP = atpFlux-TCAflux-3*glyFlux;
O2Yields = O2ATP./O2flux/2;
O2Yields(end) = 0; %assume that 0/0 is = 0
massYield = atpFlux./massOfFlux;

exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;

subplot(3,1,1)
hold all
for i = 3:-1:1
    barh(i, glyYields(i), 'FaceColor', exMap(i,:), 'LineStyle', 'none');
end
legend(fliplr(solutionNames(1:3)))
legend boxoff
ylabel('ATP/glycogen')
set(gca,'ytick',[])
ylim([0.5 3.5])

subplot(3,1,2)
hold all
for i = 3:-1:1
    barh(i, massYield(i), 'FaceColor', exMap(i,:), 'LineStyle', 'none');
end
set(gca,'ytick',[])
ylim([0.5 3.5])

ylabel('ATP/protein')

subplot(3,1,3)
hold all
for i = 3:-1:1
    barh(i, O2Yields(i), 'FaceColor', exMap(i,:), 'LineStyle', 'none');
end
ylabel('P/O ratio')
set(gca,'ytick',[])
ylim([0.5 3.5])

%%
load('model/reducedModel')
model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = addSpecificActivityConstraint(model, 0.5, 0.0444*0.72, 60);
model = addReversedReactions(model);

massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightValue = model.b(massConstraintRow,2);
weightRow = full(model.S(massConstraintRow,:));

minimalFlux= [-1000
              -1000
              -1000];

growthRates = linspace(0, 22.1, 100);

reactionNumbers = getBounds(model, minimalMedia);
model = setParam(model, 'obj', reactionNumbers(2), 1);

objectiveFunction = {'human_ATPMaintainance'};
model = setParam(model, 'lb', 'human_ATPMaintainance', 0); 

model = configureModel(model, minimalMedia, minimalFlux);

fullSolution = runChemostatExperiment(model, growthRates, objectiveFunction);
fullMass = fullSolution .* repmat(weightRow, size(fullSolution,1),1);
massOfFlux2 = sum(fullMass,2);

%%
figure()
hold all
ATPperGly = growthRates./-fullSolution(:,reactionNumbers(2))';
ATPperProtein = growthRates./massOfFlux2';

ATPperProtein(1) = 0;
ATPperGly(1) = ATPperGly(2);

area(ATPperProtein, ATPperGly, 'FaceColor', [17 115 187]/256, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(ATPperProtein, ATPperGly, 'linewidth', 2)
ylim([0 36])
xlim([0 701])
ylabel('ATP/glycogen')
xlabel('mmol ATP/g protein/h')
scatter(massYield, glyYields, 50, 'filled')
for i = 1:length(glyYields)
    text(massYield(i), glyYields(i), solutionNames{i})
end
