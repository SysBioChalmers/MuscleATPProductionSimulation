function model = mapProteomToRxns(model, rxnFile)
proteinData = importdata(rxnFile);
rxns = proteinData.textdata;
proteinAbundance = proteinData.data;

results = -ones(length(model.eccodes), 1);

for i = 1:length(rxns)
    results(ismember(model.rxns, rxns{i})) = proteinAbundance(i);
end

model.proteinMass = results;
end

