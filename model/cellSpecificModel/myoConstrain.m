function model = myoConstrain(model)
    kockRxns = importdata('myoKnock.txt');

    affectedRxns = ismember(model.rxns, kockRxns);
    model.lb(affectedRxns) = 0;
    model.ub(affectedRxns) = 0;
    
end