function [rxnNames, rxns] = find_rxn_for_gene(model1, model2, Gene)
  
  GeneId = strmatch(Gene, model1.genes_unique_names);
  ngenes = length(model1.genes);
  x = ones(ngenes,1);
  x(GeneId) = 0;
  r1 = find(~cellfun(@isempty,strfind(model.rules,sprintf('x(%d)',GeneId)))>0);
  rxn1 = model.rules(r1);
  GeneId = strmatch(Gene, model2.genes_unique_names);
  ngenes = length(model2.genes);
  x = ones(ngenes,1);
  x(GeneId) = 0;
  r2 = find(~cellfun(@isempty,strfind(model.rules,sprintf('x(%d)',GeneId)))>0);
  rxn2 = model.rules(r2);
  
  [rxns] = setDiff(rxn2,rxn1);
   
end