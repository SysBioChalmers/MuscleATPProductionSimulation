function [model, rxnList] = addProteinConstrant(model, saturation)

rxnList= zeros(length(model.rxns),1);

for i = 1:length(model.rxns)
   curProtein = model.proteinMass(i);
   if curProtein>0
       
       curSA = model.specificActivity(i);
       if curSA>0
           capacity = saturation*curSA*60*curProtein;
           if model.lb(i)>0
               model.lb(i) = -capacity;
           end
           if model.ub(i)>0
               model.ub(i) = capacity;
           end
           rxnList(i) = 1;
       end
   end
end
rxnList = rxnList == 1;

end

