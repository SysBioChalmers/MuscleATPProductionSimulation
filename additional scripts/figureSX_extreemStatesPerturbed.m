load('../model/reducedModel')
addpath('../src1')
model = mapDataToRxns(model, '../data/RxnAndSA.txt');
model = addSpecificActivityConstraint(model, 0.5, 0.03, 60);
model = addReversedReactions(model);
massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightRow = full(model.S(massConstraintRow,:));

primColor = [0    0.4470    0.7410];
secColor = [0.8500    0.3250    0.0980];

[ATPperGly,ATPperProtein] = getParetoCurve(model);
ATPperGly(isnan(ATPperProtein)) = [];
ATPperProtein(isnan(ATPperProtein)) = [];
ATPperGly(end+1) = 0; 
ATPperProtein(end+1) = ATPperProtein(end); 

[glyYields, massYield, solutionNames, solutions] = calculateSubSolutions(model);

nrPerturbations = 10;
assumedError = 0.2;
massPerturbed= zeros(length(massYield),nrPerturbations);


for i = 1:nrPerturbations
    weights = (1 + assumedError*2*(rand(1,length(weightRow))-0.5)) .* weightRow;
    massPerturbed(:,i) = 1./(solutions*weights');
end
%%
maxPerturb = 1;
nrPerturbations = 5000;
errorValues = linspace(0,maxPerturb,50);
cplx1Fraction = zeros(length(errorValues),1);
cplx1FractionPercent = zeros(length(errorValues),2);

fermentationVSCplx1 = zeros(length(errorValues),1);
fermentationVSCplx1Percent = zeros(length(errorValues),2);

fermentationFraction = zeros(length(errorValues),1);
fermentationFractionPercent = zeros(length(errorValues),2);

curPerturbed = zeros(3,nrPerturbations);

for i = 1:length(errorValues)
    curError = errorValues(i);
    for j = 1:nrPerturbations
        perturbation = (1 + curError*2*(rand(1,length(weightRow))-0.5));
        weights =  perturbation.* weightRow;
        curPerturbed(:,j) = (solutions(1:3,:)*weights');
    end
    
    cplx1Fraction(i) = sum(curPerturbed(1,:)>curPerturbed(2,:))/nrPerturbations;
    cplx1FractionPercent(i,1) = mean(curPerturbed(1,:)./curPerturbed(2,:));
    cplx1FractionPercent(i,2) = std(curPerturbed(1,:)./curPerturbed(2,:));   

    fermentationVSCplx1(i) = sum(curPerturbed(2,:)>curPerturbed(3,:))/nrPerturbations;
    fermentationVSCplx1Percent(i,1) = mean(curPerturbed(2,:)./curPerturbed(3,:));
    fermentationVSCplx1Percent(i,2) = std(curPerturbed(2,:)./curPerturbed(3,:));     
    
    fermentationFraction(i) = sum(curPerturbed(1,:)>curPerturbed(3,:))/nrPerturbations; 
    fermentationFractionPercent(i,1) = mean(curPerturbed(1,:)./curPerturbed(3,:));
    fermentationFractionPercent(i,2) = std(curPerturbed(1,:)./curPerturbed(3,:));       
    
    
end

subplot(2,2,1)
hold all
plot(errorValues, cplx1Fraction, 'linewidth', 2)
plot(errorValues, fermentationVSCplx1, 'linewidth', 2)
plot(errorValues, fermentationFraction, 'linewidth', 2)

legend('CI+ vs CI-', 'CI- vs Ferm', 'CI+ vs Ferm', 'location', 'SW')
legend boxoff
xlabel('Perturbation factor')
ylabel('fraction of simulations')
ylim([0 1])

%set(gca,'xscale','log')
subplot(2,2,2)
errorbar(errorValues,cplx1FractionPercent(:,1), cplx1FractionPercent(:,2))
xlabel('Perturbation factor')
ylabel('factor more mass')
ylim([1 2])
subplot(2,2,3)
errorbar(errorValues,fermentationVSCplx1Percent(:,1), fermentationVSCplx1Percent(:,2))
xlabel('Perturbation factor')
ylabel('factor more mass')
ylim([1 5])
subplot(2,2,4)
errorbar(errorValues,fermentationFractionPercent(:,1), fermentationFractionPercent(:,2))
ylim([1 6])
% legend('Complex I bypass', 'Fermentation', 'location', 'NW')
% legend boxoff
xlabel('Perturbation factor')
ylabel('factor more mass')

%%
nrPerturbations = 10;
paretoYields = linspace(0, 850, 200);
catalyticYield = zeros(length(paretoYields), nrPerturbations);

for i = 1:nrPerturbations
    weights = (1 + assumedError*2*(rand(1,length(weightRow))-0.5)) .* weightRow;
    model.S(end,:) = weights;
    [AG,AP] = getParetoCurve(model);
    %standardize results by interpolating to fixed X values
    AG(isnan(AP)) = [];
    AP(isnan(AP)) = [];
    [a, b] =unique(AP, 'stable');
    AG = AG(b);
    AP = AP(b);
    catalyticYield(:,i) = interp1(AP,AG,paretoYields);
end

%%
figure()
hold all
area(ATPperProtein, ATPperGly, 'FaceColor', [17 115 187]/256, 'FaceAlpha', 0.2, 'EdgeColor', 'none')

ylim([0 36])
xlim([0 800])
ylabel('ATP/glycogen')
xlabel('mmol ATP/g protein/h')


%errorbar(mean(massPerturbed(4,:)), glyYields(4), std(massPerturbed(4,:)),'horizontal', 'kx') 

plot(ATPperProtein, ATPperGly, 'linewidth', 2, 'color', primColor)

scatter(massYield, glyYields, 50, 'filled', 'MarkerFaceColor', secColor)

for i = 1:length(glyYields)
    text(massYield(i), glyYields(i), solutionNames{i})
end

for i = 1:size(catalyticYield,2)
    plot(paretoYields, catalyticYield(:,i), 'color', [0 0 0 0.1])
end


for i = 1:size(massPerturbed,2)
    plot(massPerturbed(4,i), glyYields(4), '.', 'color', [0 0 0 0.1])
end



