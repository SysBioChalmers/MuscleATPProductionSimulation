function model = addBloodPackage(model, bloodConc, comp2, comp1)
    inMets = {
        ['O2[' comp1 ']']
        ['glucose[' comp1 ']']
        ['L-lactate[' comp2 ']']
        ['CO2[' comp2 ']']
        ['human_biomass[' comp1 ']']
        };
    
    outMets = {
        ['O2[' comp2 ']']
        ['glucose[' comp2 ']']
        ['L-lactate[' comp1 ']']
        ['CO2[' comp1 ']']
        };
    inMets
    outMets
    bloodEquation = sprintf('%2.2f %s', bloodConc(1), inMets{1});
    for i = 2:length(inMets)
        bloodEquation = sprintf('%s + %2.2f %s', bloodEquation, bloodConc(i), inMets{i});
    end

    bloodEquation = sprintf('%s => %2.2f %s', bloodEquation, bloodConc(1), outMets{1});
    
    for i = 2:length(outMets)
        bloodEquation = sprintf('%s + %2.2f %s', bloodEquation, bloodConc(i), outMets{i});
    end
    
    equationStruct = createRXNStuct(model, 'BloodPackage', bloodEquation, 0, 1000, 'Blood Transport');
    model=addRxns(model,equationStruct,3,'sb',false);
    
    for i = 1:length(outMets)
         rxnName = sprintf('BloodBackflux%i', i);
         rxnEquation = sprintf('%s => %s', outMets{i}, inMets{i});
         equationStruct = createRXNStuct(model, rxnName, rxnEquation, 0, 1000, 'Blood Transport');
         model=addRxns(model,equationStruct,3,'sb',false);       
    end
end

