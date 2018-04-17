function model = changeBiomassEquation(model, input, output)
reactionName =  'human_biomass';
biomassRxn = findIndex(model.rxns, reactionName);

%reset Equation
model.S(:,biomassRxn) = 0;

for i = 1:length(input)
    model = configureGrowthRelatedMaintainance(model, 1, reactionName, input{i});
end

for i = 1:length(output)
    model = configureGrowthRelatedMaintainance(model, -1, reactionName, output{i});
end

model = configureGrowthRelatedMaintainance(model, -1, reactionName, 'human_biomass[c]');


end

