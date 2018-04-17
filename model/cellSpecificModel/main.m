addpath('C:\Program Files\SBML\libSBML-5.10.2-libxml2-x64\bindings\matlab\matlab')
modelMyo=importModel('iMyocyte2419.xml');

load('../genericHuman2.mat')

%ignore knocks that have no gene association
myoRxns = modelMyo.rxns;
genRxns = model.rxns;
myoKnock = setdiff(genRxns, myoRxns);

%remove aquaporin mitochondrial anotation:
model.grRules{findIndex(model.rxns,'HMR_4888')} = '';
model.grRules{findIndex(model.rxns,'HMR_4951')} = '';

noRules = model.rxns(ismember(model.grRules,''));
myoKnock = setdiff(myoKnock, noRules);

fid = fopen('myoKnock.txt','w');
fprintf(fid,'%s\n', myoKnock{:});
fclose(fid);

rxns = constructEquations(model, myoKnock);

fid = fopen('myoKnockRxns.txt','w');
fprintf(fid,'%s\n', rxns{:});
fclose(fid);
