function model = configureSMatrix(model, value, reaction, metabolite)
% configureSMatrix
% Adds a value to the S matrix using rxn and metabolite names as reference
%
%   model             a model structure
%   value             the stocheometry of the metabolite substrates are
%                     given positive values and products negative
%   reaction          the name of the rxn e.g. 'human_biomass'
%   metabolite        the name of the metabolite e.g. 'ATP[c]'
%   model             an updated model structure
%
%   Avlant Nilsson, 2016-05-16
%
    reaction = findIndex(model.rxns, reaction);
    metaboliteNames = modifyMetNames(model);
    metabolite = findIndex(metaboliteNames, metabolite);
    model.S(metabolite, reaction) = -value;
end

