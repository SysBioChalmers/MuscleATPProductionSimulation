function model = addExchangeMedum(model)
%Minimal media simulation
inMedia = {
    'glycogen[s]'
    'palmitate[s]'
    'H2O[s]'
    'CO2[s]'
    'O2[s]'
    %'glutamate[s]'
};
constraints =  [
    -1000
    -1000
    -1000
    -1000
    -0
    %-10000
];

outMedia = {
    'H2O[s]'
    'stearate[s]'
    'glycogen[s]'
    'O2[s]'
    'CO2[s]'
    %'alanine[s]'
    'ATPwork[s]'
    'ATPmaintainance[s]'
};


outConstraints = [
    1000
    1000
    1000
    1000
    0
    %1000
    1000
    1000
    ];


[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

reactionIn= getBounds(model, inMedia);
reactionOut = getBounds(model, outMedia);

model = setParam(model, 'lb', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', exchangeRxnsIndexes, 0);

model = setParam(model, 'lb', reactionIn, constraints);
model = setParam(model, 'ub', reactionOut, outConstraints);




end

