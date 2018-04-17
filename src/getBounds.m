function result = getBounds(model, foodLabels)
    metaboliteNames = modifyMetNames(model);

    metaboliteNumbers = getIndexFromText(metaboliteNames, foodLabels);
    
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
    exchangeS = full(model.S(:, exchangeRxnsIndexes)); 
    result = zeros(length(metaboliteNumbers),1);
    
    for i = 1:length(metaboliteNumbers)
        result(i) = -1;    
        
        if metaboliteNumbers(i) ~= -1
            
            currentReaction = find(exchangeS(metaboliteNumbers(i), :));
            
            if sum(currentReaction) ~= 0
                result(i) = exchangeRxnsIndexes(currentReaction);
            end
        end
        if result(i) == -1
           fprintf('Error with %s\n', foodLabels{i});
            
        end
    end
end

