function model = modify_model(model, rxn_names, rxns, rxns_to_use, add)
   if add
     for i = 1:rxns_to_use  
       model = addReaction(model,rxn_names{i},rxns{i}); %% Usage =>  model = addReaction(model,'newRxn1','A -> B + 2 C')    
     end  
   else
       model = removeRxns(model,rxn_names(1:rxns_to_use));
   end

end