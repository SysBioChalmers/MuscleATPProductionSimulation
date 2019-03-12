close all
hold all
addpath('src1')
proteinContent = 0.03;
yeastSA = 500;

load('model/reducedModel')
yeastModel = model;

%Make reference model
model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = addSpecificActivityConstraint(model, 0.5, proteinContent, 60);
model = addReversedReactions(model);
[ATPperCmol, ATPperProtein, cMol, fullSolution] = getParetoCurve(model);
GLYPHOS = findIndex(model.rxns, 'HMR_0483_back');

%Add yeast NDI
glycogenPhos = createRXNStuct(yeastModel, 'yeastNDI', 'H+[m] + NADH[m] + ubiquinone[m] => NAD+[m] + ubiquinol[m]', 0, 1000, 'Starch and sucrose metabolism');
yeastModel=addRxns(yeastModel,glycogenPhos,3,'m',false);
yeastRxn = findIndex(yeastModel.rxns, 'yeastNDI');

yeastModel = mapDataToRxns(yeastModel, 'data/RxnAndSA.txt');
yeastModel.specificActivity(yeastRxn) = yeastSA;

yeastModel = addSpecificActivityConstraint(yeastModel, 0.5, proteinContent, 60);
yeastModel = addReversedReactions(yeastModel);
[ATPperCmolYeast, ATPperProteinYeast, cMolYeast, fullSolutionYeast] = getParetoCurve(yeastModel);


%%
hold all
subplot(1,3,1)
hold all
plot(ATPperProtein, ATPperCmol, 'linewidth', 2)
plot(ATPperProteinYeast, ATPperCmolYeast, 'linewidth', 2)
legend('w/o NDI', 'with NDI')
legend boxoff
ylim([0 7])
xlim([0 701])
ylabel('substrate efficency [ATP/Cmol]')
xlabel('catalytic capacity [mmol ATP/g protein/h]')

subplot(1,3,2)
hold all
plot(ATPperProtein, fullSolution(:,GLYPHOS)./proteinContent, 'linewidth', 2)
plot(ATPperProteinYeast, fullSolutionYeast(:,GLYPHOS)./proteinContent, 'linewidth', 2)
ylabel('GLY-PHOS flux [mmol/g protein/h')
xlabel('catalytic capacity [mmol ATP/g protein/h]')
ylim([0 120])

subplot(1,3,3)
hold all
plot(ATPperProtein, zeros(length(ATPperProtein),1), 'linewidth', 2)
plot(ATPperProteinYeast, fullSolutionYeast(:,yeastRxn)./proteinContent, 'linewidth', 2)
ylabel('NDI flux [mmol/g protein/h')
xlabel('catalytic capacity [mmol ATP/g protein/h]')
ylim([0 120])

