close all
hold all

load('model/reducedModel')
addpath('src1')

model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = addSpecificActivityConstraint(model, 0.5, 0.03, 60);
model = addReversedReactions(model);
massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightRow = full(model.S(massConstraintRow,:));

[glyYields, massYield, solutionNames, solutions] = calculateSubSolutions(model);
[ATPperGly,ATPperProtein] = getParetoCurve(model);

model.ub(findIndex(model.rxns, 'HMR_0483_back')) = 0;
[ATPperGlyNoCIbypass,ATPperProteinNoCIbypass] = getParetoCurve(model);

ATPperGly(end+1) = 0;
ATPperProtein(end+1) = ATPperProtein(end);


ATPperGlyNoCIbypass(isnan(ATPperGlyNoCIbypass)) = [];
ATPperProteinNoCIbypass(isnan(ATPperProteinNoCIbypass)) = [];
ATPperGlyNoCIbypass(end+1) = 0;
ATPperProteinNoCIbypass(end+1) = ATPperProteinNoCIbypass(end);

area(ATPperProtein, ATPperGly, 'FaceColor', [17 115 187]/256, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
fill([ATPperProteinNoCIbypass flip(ATPperProtein)], [ATPperGlyNoCIbypass flip(ATPperGly)], [17 115 187]/256, 'FaceAlpha', 0.5,  'EdgeColor', 'none')
plot(ATPperProtein, ATPperGly, 'linewidth', 2, 'color', [17 115 187]/256)

ylim([0 7])
xlim([0 701])
ylabel('ATP/Cmol')
xlabel('mmol ATP/g protein/h')
scatter(massYield, glyYields, 50, 'filled', 'markerfacecolor', [17 115 187]/256)
for i = 1:length(glyYields)
    text(massYield(i), glyYields(i), solutionNames{i})
end

%%
figure()
hold all
X = 1./ATPperProtein;
Y =  1./ATPperGly;
X(end) = [];
Y(end) = [];
Ymax = round(max(Y)*1.3,1);

area([X min(X)*0.999 0], [Y Ymax Ymax], 'FaceColor', [150 150 150]/256, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(1./ATPperProtein, 1./ATPperGly, 'linewidth', 2)
plot(1./ATPperProteinNoCIbypass, 1./ATPperGlyNoCIbypass, '--', 'linewidth', 2, 'color', [17 115 187]/256)


scatter(1./massYield, 1./glyYields, 50, 'filled', 'markerfacecolor', [17 115 187]/256)
for i = 1:length(glyYields)
    text(1./massYield(i), 1./glyYields(i), solutionNames{i})
end
xlim([0 10^-2])
ylim([0 Ymax])

