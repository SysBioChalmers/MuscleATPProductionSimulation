function [fullSolution, ATPrate]  = optimizeFluxes(model, settings, optimization)
if strcmp(optimization, 'quad')
    settings.minRxns = {'glycogen_Exchange_'};
    settings.minVal = [1];
    settings.QminRxns = {'ventilation', 'HMR_9135_m2'};
    settings.QminVal = [0.05 0.05];
    [fullSolution, ATPrate] = runFullModel(model, settings);
elseif strcmp(optimization, 'uncoupling')
    settings.minRxns = {'glycogen_Exchange_', 'ventilation', 'HMR_9135_m2'};
    settings.minVal = [1, -0.5, -1];
    [fullSolution, ATPrate] = runFullModel(model, settings);
elseif strcmp(optimization, 'RER')
    settings.minRxns = {'glycogen_Exchange_', 'ventilation', 'HMR_9135_m2'};
    settings.minVal = [1, -0.5, -0.01];
    [fullSolution, ATPrate] = runFullModel(model, settings);
elseif strcmp(optimization, 'AandB')
    settings.minRxns = {'glycogen_Exchange_', 'ventilation', 'HMR_9135_m2', 'HMR_9135_m3'};
    settings.minVal = [1, -1, -0.4, 0.01];
    [fullSolution, ATPrate] = runFullModel(model, settings);
end

end

