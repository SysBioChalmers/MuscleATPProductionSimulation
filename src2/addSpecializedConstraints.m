function model = addSpecializedConstraints(model, ventilation, maintainance, energy)
%Allow transport of objective metabolite:
objRxn = findIndex(model.rxns, 'JoinMuscleATP1');
model.ub(objRxn) = energy;
model.lb(objRxn) = 0;

%Allow maintainance
objRxn = findIndex(model.rxns, 'Maintainance');
model.lb(objRxn) = maintainance;
model.ub(objRxn) = 1000;

%Allow ventilation
objRxn = findIndex(model.rxns, 'ventilation');
model.ub(objRxn) = ventilation;
model.lb(objRxn) = 0;
end

