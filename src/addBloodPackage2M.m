function model = addBloodPackage2M(model, bloodConc, bloodName, cellComp, bloodComp)
    inMets = {
        ['O2[' bloodComp ']']
        ['glucose[' bloodComp ']']
        ['L-lactate[' cellComp ']']
        ['CO2[' cellComp ']']
        };
    
    outMets = {
        ['O2[' cellComp ']']
        ['glucose[' cellComp ']']
        ['L-lactate[' bloodComp ']']
        ['CO2[' bloodComp ']']
        [bloodName '[' bloodComp ']']        
        };

    bloodEquation = sprintf('%2.2f %s', bloodConc(1), inMets{1});
    for i = 2:length(inMets)
        bloodEquation = sprintf('%s + %2.2f %s', bloodEquation, bloodConc(i), inMets{i});
    end

    bloodEquation = sprintf('%s => %2.2f %s', bloodEquation, bloodConc(1), outMets{1});
    
    for i = 2:length(outMets)
        bloodEquation = sprintf('%s + %2.2f %s', bloodEquation, bloodConc(i), outMets{i});
    end
    
    equationStruct = createRXNStuct(model, bloodName, bloodEquation, 0, 1000, 'Blood Transport');
    model=addRxns(model,equationStruct,3,'sb',true);
    
    for i = 1:length(inMets)
         rxnName = sprintf('BloodBackflux%s%i', bloodName, i);
         rxnEquation = sprintf('%s => %s', outMets{i}, inMets{i});
         equationStruct = createRXNStuct(model, rxnName, rxnEquation, 0, 1000, 'Blood Transport');
         model=addRxns(model,equationStruct,3,'sb',false);       
    end
end

