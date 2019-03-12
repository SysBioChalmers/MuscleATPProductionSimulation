function plotPathwayMap(model, ATPrate, fullSolution, settings)



if settings.pfba == true
    threshold = 10^-6;
    
    absFlux = abs(fullSolution)';
    roundedWatt = round(1000*molToW(ATPrate));
    roundedWatt = cellstr(num2str(roundedWatt'))';
    
    for j = 2:(length(roundedWatt)-1)
        if mod(j,2) ==0
            roundedWatt{j} = '';
        end
    end    
    
    %RemoveTissueSpecificNess
    for i = 1:length(model.subSystems)
        %tmp = split(model.subSystems{i}, '_');
        tmp = strrep(model.subSystems{i}, '_', '-');
        model.subSystems{i} = tmp; 
    end
    
    
    allSubsystems = unique(model.subSystems);
    subsystemFlux = zeros(length(allSubsystems), length(ATPrate));

    for i = 1:length(allSubsystems)
        curRxns = ismember(model.subSystems, allSubsystems{i});
        subsystemFlux(i,:) = sum(absFlux(curRxns, :), 1);
    end


    emptySubsystems = sum(subsystemFlux,2) <threshold;
    
    subsystemFlux(emptySubsystems,:) = [];
    allSubsystems(emptySubsystems) = [];
    
    %Remove muscle 3, artifical, etc
    emptySubsystems = contains(allSubsystems, {'Muscle3', 'blood', 'artificial', 'Exchange', 'Transport', 'Artificial', 'Ventilation', 'Objective Funcion'});
    subsystemFlux(emptySubsystems,:) = [];
    allSubsystems(emptySubsystems) = [];   

    normalizationMatrix = repmat(max(subsystemFlux')', 1,length(roundedWatt));

    CGObject = clustergram(subsystemFlux./normalizationMatrix, 'RowLabels', allSubsystems, 'ColumnLabels', roundedWatt, 'Colormap', redbluecmap, 'Cluster', 1, 'ColumnLabelsRotate', 0, 'Symmetric', false);
    plot(CGObject)
        
else
    disp('simulations have to be run with pFBA for this to be meaningful')
end


end

