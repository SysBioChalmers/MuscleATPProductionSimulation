load('../model/muscleModel.mat')
A = importdata('proteomicLean.txt');
C = importdata('ID_conversion_key.txt');
uniprotMap = containers.Map(C.textdata(:,4), C.textdata(:,1));



%Map uniprot to ensembl genes
geneList = A.textdata(:,1);
for i = 2:length(geneList)
    geneList{i} = strrep(geneList{i}, ';', '');
    if isKey(uniprotMap, geneList{i})
       geneList{i} =  uniprotMap(geneList{i});
    else
       disp(geneList{i}) 
    end
end
geneList(1) = [];

%Normalize proteomics by mass and sort
totalMass = A.data(:,1) .* A.data(:,2);
massSum = sum(totalMass);
normalizedMass = totalMass./massSum;
geneMap = containers.Map(geneList, normalizedMass);

[normalizedMass, indx]= sort(normalizedMass, 'descend');
geneList = geneList(indx);
plot(cumsum(normalizedMass), 'o-')
ylim([0 1])

%%
%Plot proteomics by subsystem
%plotAllSubs(model, geneMap, 0.02);
%inspectSubsystem(model, 'Glycine, serine and threonine metabolism', geneMap, 0.001);
%inspectReactions(model, {'HMR_4652', 'HMR_6911', 'HMR_6914', 'HMR_6916', 'HMR_6918', 'HMR_6921'},geneMap, 0.01)

%%
%Plot mass distribution
mapOfSystems = containers.Map;    
mapOfSystems('Glycolysis') = getInvolvedGenes(model, ismember(model.subSystems, 'Glycolysis / Gluconeogenesis'));
mapOfSystems('OXPHOS') = getInvolvedGenes(model, ismember(model.subSystems, 'Oxidative phosphorylation'));
mapOfSystems('TCA') = getInvolvedGenes(model, ismember(model.subSystems, 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism'));

mapOfSystems('Fatty Acids') = getInvolvedGenes(model, contains(model.subSystems,'fatty acid'));

mapOfSystems('Myosin') = {'ENSG00000092054' 'ENSG00000125414' 'ENSG00000180209' 'ENSG00000168530' 'ENSG00000111245' 'ENSG00000109061' 'ENSG00000160808' 'ENSG00000086967' 'ENSG00000198467' 'ENSG00000196465' 'ENSG00000143549' 'ENSG00000101306' 'ENSG00000144821', 'ENSG00000101470', 'ENSG00000175084', 'ENSG00000177791', 'ENSG00000130595','ENSG00000114854','ENSG00000130598','ENSG00000105048'};
mapOfSystems('Glycogen phosphorylase') = {'ENSG00000068976'};
mapOfSystems('Creatine') = {'ENSG00000104879' 'ENSG00000131730' 'ENSG00000166165'};
mapOfSystems('Albumin') = {'ENSG00000163631'};
mapOfSystems('Myoglobin') = {'ENSG00000198125'};
%mapOfSystems('Aktin') = {'ENSG00000075624', 'ENSG00000184009', 'ENSG00000077522', 'ENSG00000143632', 'ENSG00000248746', 'ENSG00000154553', 'ENSG00000120729', 'ENSG00000035403'};
%mapOfSystems('Heatshock proteins') = {'ENSG00000109971', 'ENSG00000126803', 'ENSG00000115541', 'ENSG00000170606', 'ENSG00000144381', 'ENSG00000126602'};

%Add Aktin to myosin
mapOfSystems('Myosin') = [mapOfSystems('Myosin') 'ENSG00000075624', 'ENSG00000184009', 'ENSG00000077522', 'ENSG00000143632', 'ENSG00000248746', 'ENSG00000154553', 'ENSG00000120729', 'ENSG00000035403'];

%add malate aspartate shuttle and glyceron phosphate shuttle to TCA 
mapOfSystems('TCA') = [mapOfSystems('TCA') 'ENSG00000120053', 'ENSG00000125166', 'ENSG00000108528', 'ENSG00000167588'];

%add cytochrome C to OXPHOS
%mapOfSystems('OXPHOS') = [mapOfSystems('OXPHOS') 'ENSG00000172115'];


%ENSG00000151729 %transporter of ATP from mitochondria
%ENSG00000112992 %NNT

allKeys = mapOfSystems.keys;

allGenes = [];
results = zeros(length(allKeys),1);
for i = 1:length(allKeys)
    curGenes = mapOfSystems(allKeys{i});
    uniqueGenes = setdiff(curGenes, allGenes);
    allGenes = [allGenes uniqueGenes];
    results(i) = sumOfGenes(geneMap, curGenes);   
end
[results, indx] = sort(results, 'descend');
allKeys = allKeys(indx);

other = (1- sum(results));
results = [results; other];
allKeys = [allKeys, 'Other'];
plotPie(results, allKeys)

uniqueGenes = setdiff(geneList, allGenes);
clc
printGeneList(geneMap, uniqueGenes)

%%
ATPgenes = getInvolvedGenes(model, ismember(model.rxns, 'HMR_6916'));
revMap = containers.Map(C.textdata(:,1), C.textdata(:,4));
ATPUniprot = ATPgenes;
for i = 1:length(ATPgenes)
    if isKey(revMap, ATPgenes{i})
        if strcmp(revMap(ATPgenes{i}), '')
            disp(ATPgenes{i})
        else
            ATPUniprot{i} = revMap(ATPgenes{i});
        end
    else
        disp(ATPgenes{i})
    end
end


for i = 1:18
    ATPUniprot{i} = uniprotMap(ATPUniprot{i})
end

%%
load('../model/reducedModel.mat')
result = zeros(length(model.rxns),1);
eqns = constructEquations(model);

detectionLimit = min(A.data(:,2));

%Asign genes from the NADPH version of NADH rxns
model.grRules{findIndex(model.rxns, 'HMR_3957')} = 'ENSG00000182054';
model.grRules{findIndex(model.rxns, 'HMR_6510')} = 'ENSG00000090013';


for i = 1:length(model.rxns)
    curGenes = getInvolvedGenes(model, ismember(model.rxns, model.rxns{i}));
    results(i) = sumOfGenes(geneMap, curGenes);   
end

%Did not detect HMR_4391	DHAP[c] <=> GAP[c]
results(findIndex(model.rxns, 'HMR_4391')) = (detectionLimit * 30791)/massSum;

clc
for i = 1:length(model.rxns)
    fprintf('%s\t%s\t%2.6f\n', model.rxns{i}, eqns{i}, results(i))
end

