function [rxnNames, rxn_mets, rxns] = find_rxn_for_gene(model1, model2, Gene, add_or_del)
  
 if add_or_del
  GeneId = strmatch(Gene, model2.genes_name);
  rxn_ids2 = get_critical_rxns(model2, GeneId); 
  rxnName = model2.rxnNames(rxn_ids2);
  n_reactions = length(rxn_ids2);
  
  k = 0;
  
  for i = 1:n_reactions
    if isempty(find(model1.recon2RxnMap == rxn_ids2(i), 1)) == 1  
     rxnS = recon2.S(:,rxn_ids2(i));
     if isempty(rxnS)
         continue;
     end    
     mets = find(rxnS);
     mets1 = arrayfun(@(x) find(model1.recon2MetMap == x, 1), mets);
     if sum(isempty(mets1)) > 0
      fprintf('Need to add metabolites too\n');   
      continue;
     end 
     k = k+1;
     rxnNames(K) = rxnName(i);
     rxn_mets{k} = {model1.metNames(mets1)'};
     rxns{k} = rxnS(mets)';      
    end  
  end  
  
  if K == 0
      rxns =[];
      rxn_mets =[];
      rxnNames =[];
  end    
 else
    GeneId = strmatch(Gene, model1.genes_unique_names);
    rxn_ids1 = get_critical_rxns(model1, GeneId); 
    rxnNames = model2.rxnNames(rxn_ids1); 
    rxns =[];
    rxn_mets =[];
 end    
end


function crxn_ids = get_critical_rxns(model, GeneId)
  crxn_ids = [];
  ngenes = length(model.genes);
  x = ones(ngenes,1);
  x(GeneId) = 0;
  rxn_ids = find(~cellfun(@isempty,strfind(model.rules,sprintf('x(%d)',GeneId)))>0);
  sz = length(rxn_ids);
  for i = 1:sz
     if eval(model.rules{rxn_ids(i)}) < 1
         crxn_ids = [crxn_ids, rxn_ids(i)];
     end    
  end    
end