function metNames = modifyMetNames(modelMod)
% modifyMetNames
% Adds the compartment to the metabolite names e.g. 'glucose' in
% compartment s becomes 'glucose[s]'.
%
%   modelMod         a model structure
%   metNames         metabolite names with compartment indications
%
%   Avlant Nilsson, 2016-05-16
%
    metNames = modelMod.metNames;
    for i=1:length(metNames)
        metNames{i} = [metNames{i} '[' modelMod.compNames{modelMod.metComps(i)} ']'];
    end
end

