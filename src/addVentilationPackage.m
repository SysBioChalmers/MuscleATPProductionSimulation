function model = addVentilationPackage(model, bloodConc, comp2, comp1)
    inMets = {
        ['O2[' comp1 ']']
        ['CO2[' comp2 ']']
        };
    
    outMets = {
        ['O2[' comp2 ']']
        ['CO2[' comp1 ']']
        };

    ventEquation = sprintf('%2.2f %s', bloodConc(1), inMets{1});
    for i = 2:length(inMets)
        ventEquation = sprintf('%s + %2.2f %s', ventEquation, bloodConc(i), inMets{i});
    end

    ventEquation = sprintf('%s => %2.2f %s', ventEquation, bloodConc(1), outMets{1});
    
    for i = 2:length(outMets)
        ventEquation = sprintf('%s + %2.2f %s', ventEquation, bloodConc(i), outMets{i});
    end
    
    equationStruct = createRXNStuct(model, 'VentPackage', ventEquation, 0, 1000, 'Lung Transport');
    model=addRxns(model,equationStruct,3,'sb',false);
    
    for i = 1:length(outMets)
         rxnName = sprintf('VentBackflux%i', i);
         rxnEquation = sprintf('%s => %s', outMets{i}, inMets{i});
         equationStruct = createRXNStuct(model, rxnName, rxnEquation, 0, 1000, 'Lung Transport');
         model=addRxns(model,equationStruct,3,'sb',false);       
    end
end

