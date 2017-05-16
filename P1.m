load('Class_files/Achiles.mat');
load('Class_files/CCLE.mat');
load('Class_files/recon1.mat');
load('Class_files/recon2.mat');
changeCobraSolver('tomlab_cplex','QP');
changeCobraSolver('tomlab_cplex','LP');

% Declarations
ge_threshold = 9;
celline_id = 32;
cellline = Achiles.cellines(celline_id);

id_A = find(~cellfun(@isempty,strfind(Achiles.cellines,'SKIN'))>0);

%% Initial Model
%for i=4:4 %length(id_A)
 i = 4; 
 eGenes{i} = essGenes(Achiles,CCLE,Achiles.cellines(id_A(i)),ge_threshold);
 recon1_m = defineHumanMediaRPMI(recon1);  
 EG_r1m_FVA_id = essFVA(recon1_m);
 EG_r1m_FVA{i} = EG_r1m_FVA_id;
 EG_r1m{i} = recon1_m.genes_unique_names(unique(recon1_m.genes_unique_map(EG_r1m_FVA_id>0)));
 per_r1m{i} = 100*length(intersect(EG_r1m{i},eGenes{i}))/length(eGenes{i});

%end
%% Initial accuracy
accuracy = evaluateModel(recon1_m, celline, celline_id);

%% Improving model

isstable = 0;
initial_accuracy = accuracy;
best_accuracy = initial_accuracy;

previous_accuracy = initial_accuracy;
count_nochanges = 0;
model = recon1_m;
add_if_true = 1;
g_a = {};
g_d = {};
a = 0;
d = 0;
iter = 1;
best_iter = iter;

while ~isstable && iter < 100
    
    %% option 1
   found_reactions = 0;
   
   while found_reactions == 0
      [gene, g_a, g_d] = get_next_gene(model, cellineId, g_a, g_d, add_if_true); 
      [rxnNames, meta, rxn] = find_rxn_for_gene(model,recon2, gene);
      rxns_to_use = 1; 
      rxn_used{iter} = rxns;
      operation_used{iter} = add_if_true;
      if ~isempty(rxn_used{iter})
          found_reactions = 1;
      end    
   end 
   
   model = modify_model(model, rxnNames, meta, rxns, rxns_to_use,  add_if_true); 
   EG_r1m_FVA_id = essFVA(model);
   accuracy_new = evaluateModel(model, celline);
   
   if accuracy_new > best_accuracy
       best_accuracy = accuracy_new;
       best_iter = iter;
       count_nochanges = 0;
   else 
       if accuracy_new < precious_acuracy
           add_if_true = ~add_if_true;
           count_nochanges = 0;
       else
           if accuracy_new == precious_acuracy
               count_nochanges = count_nochanges + 1;
           else
               count_nochanges = 0;
           end    
       end  
   end
   
   if count_nochanges >= 4
       isstable = 1;
   end    
   
   iter = iter + 1;
end  
    
    
    
    
    
    



  