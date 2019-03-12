function settings = setSimulationSettings(timeSteps, pFBA)
    settings = [];
    settings.timeSteps = timeSteps;
    settings.pfba = pFBA;
    %objective function
    settings.primaryObjective = 'MuscleATPOut';    
%   settings.minRxns = {'glycogen_Exchange_b', 'O2_Exchange_b',  'HMR_9135_m1' 'HMR_9135_m3', 'L-lactate_Exchange_'};
%   settings.minVal = [0, 1, 0, 0, 5];
%   settings.minRxns = {'CO2_Exchange_b', 'ATPworkOut_m1'};
%   settings.minVal = [-1, 0.01];    

    settings.minRxns = {'O2_Exchange_b', 'glycogen_Exchange_b', 'glucose_Exchange_b', 'lactateBuffering', 'ATPworkOut_m1', 'HMR_9135_m1', 'HMR_9135_m3'};
    settings.minVal = [1, 1, 1, -3, 0.01, 0.9, 3];
    
end

