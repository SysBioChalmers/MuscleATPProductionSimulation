function model = addSupportReactions(model, m1Ratio)
    m2Ratio = 1 - m1Ratio;
    reactionString = sprintf('%f ATPwork[sm1] + %f ATPwork[sm2] => ATPwork[s]', m1Ratio, m2Ratio);
    equationStruct = createRXNStuct(model, 'JoinMuscleATP1', reactionString, 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);
    equationStruct = createRXNStuct(model, 'MuscleATPOut', 'ATPwork[s] =>', 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);
    equationStruct = createRXNStuct(model, 'ventilation', 'CO2[s] => O2[s]', 0, 1000, 'Ventilation');
    model=addRxns(model,equationStruct,3,'s',true);
    equationStruct = createRXNStuct(model, 'Maintainance', 'ATPwork[sm3] => ATPmaintainance[s]', 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);
    equationStruct = createRXNStuct(model, 'MaintainanceOut', 'ATPmaintainance[s] =>', 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);    
end

