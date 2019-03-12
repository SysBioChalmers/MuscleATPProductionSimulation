function [ATPperCmol,ATPperProtein, cMol, fullSolution] = getParetoCurve(model)
massConstraintRow = findIndex(model.mets, 'MassConstraint');
weightValue = model.b(massConstraintRow,2);
weightRow = full(model.S(massConstraintRow,:));

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

ATPRates = linspace(0, 25, 1000);

reactionNumbers = getBounds(model, minimalMedia);

fatNr = reactionNumbers(2);
glucoseNr = reactionNumbers(3);
stochiometry = [8 6];

model = setParam(model, 'obj', [fatNr glucoseNr], [8 6]);
model = configureModel(model, minimalMedia, minimalFlux);

fullSolution = runChemostatExperiment(model, ATPRates, 'human_ATPMaintainance');
fullMass = fullSolution .* repmat(weightRow, size(fullSolution,1),1);
massOfFlux2 = sum(fullMass,2);

cMol = [-stochiometry(1)*fullSolution(:,fatNr) -stochiometry(2)*fullSolution(:,glucoseNr)];
%o2Flux = -fullSolution(:,O2Nr)';

ATPperCmol = ATPRates./sum(cMol,2)';
ATPperProtein = ATPRates./massOfFlux2';

%Specify the 0 case
ATPperProtein(1) = 0;
ATPperCmol(1) = ATPperCmol(2);

ATPperProtein(ATPperCmol == -inf) = [];
cMol(ATPperCmol == -inf,:) = [];
fullSolution(ATPperCmol == -inf,:) = [];
ATPperCmol(ATPperCmol == -inf) = [];
end

