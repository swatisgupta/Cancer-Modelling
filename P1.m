load('../Class_files/Achiles.mat');
load('../Class_files/CCLE.mat');
load('../Class_files/recon1.mat');
load('../Class_files/recon2.mat');
changeCobraSolver('tomlab_cplex','QP');
changeCobraSolver('tomlab_cplex','LP');

%% Celline & Essential Genes

ge_threshold = 9;
celline_id = 32;
celline = Achiles.cellines(celline_id);
id_A = find(~cellfun(@isempty,strfind(Achiles.cellines,'SKIN'))>0);

essG = essGenes(Achiles,CCLE,celline,ge_threshold);

%% Initial Model

recon1_m = defineHumanMediaRPMI(recon1);

%for i=4:4 %length(id_A)
% i = 4; 
% eGenes{i} = essGenes(Achiles,CCLE,Achiles.cellines(id_A(i)),ge_threshold);
% recon1_m = defineHumanMediaRPMI(recon1);  
% EG_r1m_FVA_id = essFBA(recon1_m);
% EG_r1m_FVA{i} = EG_r1m_FVA_id;
% EG_r1m{i} = recon1_m.genes_unique_names(unique(recon1_m.genes_unique_map(EG_r1m_FVA_id>0)));
% per_r1m{i} = 100*length(intersect(EG_r1m{i},eGenes{i}))/length(eGenes{i});
%
%end

%% Initial model

[acc_best, essGM_best] = evaluateModel(recon1_m, essG);

%% Improving model

isstable = 0;
acc_init = acc_best;
acc_prev = acc_best;
count_nochanges = 0;
model = recon1_m;
add_if_true = 1;
g_a = [];
g_d = [];
a = 0;
d = 0;
iter = 1;
essGM = [];

while ~isstable && iter < 4

    found_reactions = 0;
	while found_reactions == 0
        [gene,g_a,g_d] = get_next_gene(essG,essGM,g_a,g_d,add_if_true);
        [rxnNames,meta,rxn] = find_rxn_for_gene(model,recon2,gene,add_if_true);
        rxns_to_use = 1;
        %rxn_used{iter} = rxns;
        %operation_used{iter} = add_if_true;
        if ~isempty(rxn)
            found_reactions = 1;
        end
    end
   
    model = modify_model(model,rxnNames,meta,rxns,rxns_to_use,add_if_true);
	[acc, essGM] = evaluateModel(model, essG);
   
   if acc > acc_best
       acc_best = acc;
       essGM_best = essGM;
       count_nochanges = 0;
   else 
       if acc < acc_prev
           add_if_true = ~add_if_true;
           count_nochanges = 0;
       else
           if acc == acc_prev
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

%% Select a gene
function [gene, g_a, g_d] = get_next_gene(essG,essGM,g_a,g_d,add_if_true)

    if add_if_true
        from = setdiff(essG,essGM);
    else
        from = setdiff(essGM,essG);
    end
    
    gene = datasample(from,1);
    
    if add_if_true
        g_a = [g_a gene];
    else
        g_d = [g_d gene];
    end

end