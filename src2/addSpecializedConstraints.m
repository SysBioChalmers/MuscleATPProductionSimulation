function model = addSpecializedConstraints(model, ventilation, maintainance, internalWork, m1Ratio)
%Allow transport of objective metabolite:
% objRxn = findIndex(model.rxns, 'JoinMuscleATP1');
% model.ub(objRxn) = energy;
% model.lb(objRxn) = 0;
atpDisipationM1 = findIndex(model.rxns, 'human_ATPMaintainance_m1');
atpDisipationM2 = findIndex(model.rxns, 'human_ATPMaintainance_m2');
model.lb(atpDisipationM1) = internalWork * m1Ratio;
model.lb(atpDisipationM2) = internalWork * (1-m1Ratio);

%Allow maintainance
objRxn = findIndex(model.rxns, 'Maintainance');
model.lb(objRxn) = maintainance;
model.ub(objRxn) = 1000;

%Allow ventilation
objRxn = findIndex(model.rxns, 'ventilation');
model.ub(objRxn) = ventilation;
model.lb(objRxn) = 0;
end

