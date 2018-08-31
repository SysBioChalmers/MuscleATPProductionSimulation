function settings = setSimulationSettings(timeSteps, pFBA)
    settings = [];
    settings.timeSteps = timeSteps;
    settings.pfba = pFBA;
    %objective function
    settings.primaryObjective = 'MuscleATPOut';    
%     settings.minRxns = {'glycogenUsage', 'O2_Exchange_b',  'HMR_9135_m1' 'HMR_9135_m3', 'L-lactate_Exchange_'};
%     settings.minVal = [0, 1, 0, 0, 5];
    settings.minRxns = {'O2_Exchange_', 'glycogen_Exchange_', 'L-lactate_Exchange_', 'ATPworkOut_m1', 'HMR_9135_m1', 'HMR_9135_m3'};
    settings.minVal = [1, 1, -4, 0.01, 1, 3];
end

