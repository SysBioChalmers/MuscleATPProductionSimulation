function saveFluxes(filename, model, fluxes, rates, resolution)
tresh = 10^-6;
fileID = fopen(filename,'w');

includedFluxes = find(sum(abs(fluxes),2)>tresh);
width = length(rates);
height = length(includedFluxes);
samples = linspace(1, width,resolution);
samples = round(samples);

textFields = {'rxn', 'eqn', 'sub'};

fprintf(fileID,'%s\t%s\t%s', textFields{1}, textFields{2}, textFields{3});
for i = 1:length(samples)
    curRate = samples(i);
    fprintf(fileID,'\t %2.2f', rates(curRate));
end
fprintf(fileID,'\n');

allEq = constructEquations(model, includedFluxes);
allRxn = model.rxns(includedFluxes);
allSub = model.subSystems(includedFluxes);

for i = 1:height
    curRxn = includedFluxes(i);
    fprintf(fileID,'%s\t%s\t%s', allRxn{i}, allEq{i}, allSub{i});
    for j = 1:length(samples)
        curRate = samples(j);
        fprintf(fileID,'\t %2.2f', fluxes(curRxn,curRate));
    end
    fprintf(fileID,'\n');
end

fclose(fileID);

end

