function model = limitOxygenDelivery(model, maxLevel)
     %mTransp = getTransport(model, {'O2', 'L-lactate'}, 'sm2', 'sb');
     if size(model.b,2) == 1
         model.b = [model.b model.b];
     end

    model.mets = [model.mets; 'oxygenTradeof'];
    model.metNames = [model.metNames; 'oxygenTradeof'];
    model.metComps = [model.metComps; model.metComps(end)]; 
    model.b = [model.b; [0 maxLevel]];
    model.S = [model.S; zeros(1,length(model.rxns))];    
    
    model.S(end,findIndex(model.rxns,'HMR_9048_m1')) = -1; 
    model.S(end,findIndex(model.rxns,'HMR_9048_m2')) = -1;
%     

end

