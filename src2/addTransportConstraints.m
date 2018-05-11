function model = addTransportConstraints(model, tissue, mets, lb, ub)
%Internal constraints
mExchange = getTransportReactions(model, tissue, 'sb');
model.lb(mExchange) = 0;
model.ub(mExchange) = 0;

%TransportMedia
transpMedia = {
     'glycogen'
     'glucose'
     'O2'
     'CO2'
     'HCO3-'
     'H+'
     'L-lactate'
     'palmitate'
     'H2O'
};

%Open exchange reactions
mTransp = getTransport(model, transpMedia, tissue, 'sb');
model.lb(mTransp) = -1000;
model.ub(mTransp) = 1000;

%Turn of glycogen production
mTransp = getTransport(model, {'glycogen', 'glucose'}, tissue, 'sb');
model.ub(mTransp) = 0;


for i = 1:length(mets)
    mTransp = getTransport(model, mets(i), tissue, 'sb');
    model.lb(mTransp) = lb(i);
    model.ub(mTransp) = ub(i);
end

end

