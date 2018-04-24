function printExchangeFluxes(model, fullSolution)
% printExchangeFluxes
% Prints the fluxes for each reaction together with the
% name of the subsystem and the equation of the reaction.
%
%   model               a model struct
%   fullSolution        a vector with fluxes
%
%   Avlant Nilsson, 2016-05-17
%
fluxThresh = 10^-6;
    for i=1:length(model.rxns)
       if sum(abs(fullSolution(i))) > fluxThresh
           eq = constructEquations(model, i);
           fprintf('%s\t%s\t%s', model.rxns{i}, model.subSystems{i}, eq{1})
           fprintf('\t%f', fullSolution(i));
           fprintf('\n');
       end
    end
end

