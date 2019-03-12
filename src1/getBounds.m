function result = getBounds(model, foodLabels)
% getBounds
% returns the indexes to the exchange reactions corresponding to the
% metabolites in foodLabels
%
%   model           a model struct
%   foodLabels      a cell array with metabolite names, e.g. 'glucose[s]'
%   result          a array with indexes to exchange rxns
%
%   Avlant Nilsson, 2016-05-16
%

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