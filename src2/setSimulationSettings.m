function settings = setSimulationSettings(timeSteps, pFBA)
    settings = [];
    settings.timeSteps = timeSteps;
    settings.pfba = pFBA;
    %objective function
    settings.primaryObjective = 'MuscleATPOut';    
    settings.minRxns = {'glycogen_Exchange_', 'O2_Exchange_b',  'HMR_9135_m1', 'HMR_9135_m2', 'HMR_9135_m3'};
    settings.minVal = [1, 1, 1 -1, 1];

end

