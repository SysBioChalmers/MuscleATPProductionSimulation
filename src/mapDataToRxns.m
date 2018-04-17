function model = mapDataToRxns(model, rxnFile)
SAData = importdata(rxnFile);
rxns = SAData.textdata;
specificActivity = SAData.data;

results = -ones(length(model.eccodes), 1);
medianRate = median(specificActivity);

for i = 1:length(rxns)
    currentRxn = ismember(model.rxns, rxns{i});

    if specificActivity(i) == -1
        results(currentRxn) = medianRate;
    else
        results(currentRxn) = specificActivity(i);
    end
end

model.specificActivity = results;
end

