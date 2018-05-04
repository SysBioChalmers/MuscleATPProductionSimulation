%This loads the model, makes some changes (adds rxns, makes the rxn
%directionality uniform etc.) and saves the result as a matlab struct for
%future use.

%load Raven Model
addpath('../src1')
model = importExcelModel('HMRdatabase2_00.xlsx');
addpath('cellSpecificModel')

%Standardize the direction of exchange fluxes, negative flux = uptake
model = createConsistentReactionDirection(model);

%Make sure that NADPH is used in anabolism and NADH produced in catabolism
model = fixNADHDirectionality(model);

%Fix P/O ratio
model = configureSMatrix(model, 3, 'HMR_6916', 'H+[c]');
model = configureSMatrix(model, -3, 'HMR_6916', 'H+[m]');

%Add fat reactions
model = addFatExchange(model, 's', 'AddedFatExchange');

%Allow reversed IDH flux in the cytoplasm
model.lb(findIndex(model.rxns, 'HMR_0710')) = -1000;

%Add LDHA to cytosolic gene list.
model.grRules{findIndex(model.rxns, 'HMR_4388')} = '(ENSG00000111716 or ENSG00000151116 or ENSG00000166796 or ENSG00000166800 or ENSG00000171989 or ENSG00000134333)';


%Free lactate Transport
lactRxn = createRXNStuct(model, 'FreeLactateTransport', 'L-lactate[s] <=> L-lactate[c]', -1000, 1000, 'Transport, extracellular');
model=addRxns(model,lactRxn,3,'c',false);

%Allow export of alanine from mitochondria
model = setParam(model, 'lb', 'HMR_5113', -1000);





model = myoConstrain(model)

blockedRxns = and(model.ub==0, model.lb==0);

model=removeRxns(model,blockedRxns,true,true,true);

%remove lactate shunt
model.lb(findIndex(model.rxns, 'HMR_4280')) = 0;
model.ub(findIndex(model.rxns, 'HMR_4280')) = 0;

%Allow Glycerol phosphate shuttle
model.lb(findIndex(model.rxns, 'HMR_0483')) = -1000;

%Prevent Proline->ubiqinol cycle
model.ub(findIndex(model.rxns, 'HMR_3837')) = 0;
model.lb(findIndex(model.rxns, 'HMR_3838')) = 0;
model.ub(findIndex(model.rxns, 'HMR_3838')) = 0;
model.ub(findIndex(model.rxns, 'HMR_8611')) = 0;

%Remove reversed Succinate fumarate loop
model = setParam(model, 'ub', 'HMR_8743', 0);

%Remove ubiqinol FAD miss anotations 
model.lb(findIndex(model.rxns, 'HMR_3783')) = 0; %2-methylbutyryl-CoA[m] + ubiquinone[m] <=> tiglyl-CoA[m] + ubiquinol[m]
model.ub(findIndex(model.rxns, 'HMR_3783')) = 0;
model.lb(findIndex(model.rxns, 'HMR_3751')) = 0; %isobutyryl-CoA[m] + ubiquinone[m] => methacrylyl-CoA[m] + ubiquinol[m]
model.ub(findIndex(model.rxns, 'HMR_3751')) = 0;
model.lb(findIndex(model.rxns, 'HMR_3769')) = 0; %3-methylcrotonyl-CoA[m] + ubiquinol[m] <=> isovaleryl-CoA[m] + ubiquinone[m]
model.ub(findIndex(model.rxns, 'HMR_3769')) = 0;
model.lb(findIndex(model.rxns, 'HMR_4242')) = 0; %glutaryl-CoA[m] + ubiquinone[m] => CO2[m] + crotonyl-CoA[m] + ubiquinol[m]
model.ub(findIndex(model.rxns, 'HMR_4242')) = 0;
model.lb(findIndex(model.rxns, 'HMR_3212')) = 0; %propanoyl-CoA[m] + ubiquinone[m] => acrylyl-CoA[m] + ubiquinol[m]
model.ub(findIndex(model.rxns, 'HMR_3212')) = 0;

%remove hypothetical lactade dehydrogenase to ferrocytochrome
model = setParam(model, 'lb', 'HMR_8514', 0);
model = setParam(model, 'ub', 'HMR_8514', 0);
model = setParam(model, 'lb', 'HMR_3859', 0);
model = setParam(model, 'ub', 'HMR_3859', 0);

%Add glycogen consumption
glycogenPhos = createRXNStuct(model, 'GlycogenPhosphorylase', 'glycogen[c] + Pi[c] => glucose-1-phosphate[c]', 0, 1000, 'Starch and sucrose metabolism');
model=addRxns(model,glycogenPhos,3,'c',false);
model.grRules{findIndex(model.rxns, 'GlycogenPhosphorylase')} = 'ENSG00000068976';


%Remove cytosolic oxidation
O2mets = ismember(model.metNames, 'O2');
O2rxns = sum(abs(model.S(O2mets,:)))>0;

%Ignore fatty acid oxidation
O2rxns(contains(model.subSystems,'fatty acid')) = 0;

%Ignore transport and oxphos
O2rxns(ismember(model.rxns, {'HMR_4896', 'HMR_4898',  'HMR_6914', 'HMR_9048'})) = 0;
model.lb(O2rxns) = 0;
model.ub(O2rxns) = 0;


%Overwrite raven model
save('muscleModel', 'model')
