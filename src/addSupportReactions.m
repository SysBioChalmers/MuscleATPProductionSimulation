function model = addSupportReactions(model, m1Ratio)
    m2Ratio = 1 - m1Ratio;
    reactionString = sprintf('%f human_biomass[sm1] + %f human_biomass[sm2] => obectiveMetabolite[s]', m1Ratio, m2Ratio);
    equationStruct = createRXNStuct(model, 'JoinMuscleATP1', reactionString, 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);
    equationStruct = createRXNStuct(model, 'MuscleATP', 'obectiveMetabolite[s] =>', 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);
    equationStruct = createRXNStuct(model, 'ventilation', 'CO2[s] => O2[s]', 0, 1000, 'Ventilation');
    model=addRxns(model,equationStruct,3,'s',true);
    equationStruct = createRXNStuct(model, 'Maintainance', 'human_biomass[sm3] =>', 0, 1000, 'Objective Funcion');
    model=addRxns(model,equationStruct,3,'sb',true);
end

