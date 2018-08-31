function model = addMaintenanceConstraints(model, ventilation, maintainance, internalWork)
%assume the internal work is preformed by type 1 muscle
atpDisipationM1 = findIndex(model.rxns, 'human_ATPMaintainance_m1');
model.lb(atpDisipationM1) = internalWork; 

%Allow maintainance
objRxn = findIndex(model.rxns, 'Maintainance');
model.lb(objRxn) = maintainance;
model.ub(objRxn) = 1000;

%Allow ventilation
objRxn = findIndex(model.rxns, 'ventilation');
model.ub(objRxn) = ventilation;
model.lb(objRxn) = 0;
end

