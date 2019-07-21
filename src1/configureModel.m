function model = configureModel(model, minimalMedia, minimalFlux)
maintenance = 0;


[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

%Minimal media simulation

reactionNumbers= getBounds(model, minimalMedia);

model = setParam(model, 'lb', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', exchangeRxnsIndexes, 1000);


model = setParam(model, 'lb', reactionNumbers, minimalFlux);
model = setParam(model, 'ub', reactionNumbers, 1000);


model = setParam(model, 'lb', 'human_ATPMaintainance', maintenance);
model = setParam(model, 'ub', 'human_ATPMaintainance', maintenance);
end

