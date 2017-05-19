load('../Class_files/Achiles.mat');
load('../Class_files/CCLE.mat');
load('../Class_files/recon1.mat');
load('../Class_files/recon2.mat');
load('../Class_files/HMR2.mat');
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
model = recon1_m;
rxnValid = 1:length(HMR2.rxnNames);
[acc_best, essGM_best] = evaluateModel(model, essG);

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
essGM = essGM_best;

while ~isstable && iter < 10

%     found_reactions = 0;
% 	while found_reactions < 5  %#ok<ALIGN>
%        [gene,g_a,g_d] = get_next_gene(essG,essGM,g_a,g_d,add_if_true);
%        fprintf('Gene %s (Add %d)\n',gene{1},add_if_true);
%        [rxnNames,meta,rxns] = find_rxn_for_gene(model,recon2,gene,add_if_true);
%         rxns_to_use = 1;
%         %rxn_used{iter} = rxns;
%         %operation_used{iter} = add_if_true;
%         if ~isempty(rxnNames)
%             found_reactions = found_reactions + 1;
%             model = modify_model(model,rxnNames,meta,rxns,rxns_to_use,add_if_true);
%             add_if_true = ~add_if_true;
%         end
%     end

    [model,rxnValid] = increaseModel(model,HMR2,rxnValid);
    model = decreaseModel(model);
	[acc, essGM] = evaluateModel(model, essG);
    fprintf('Iteration %d (%f)\n',iter,acc);
   
    if acc > acc_best
        acc_best = acc;
        essGM_best = essGM;
        count_nochanges = 0;
%     else 
%        if acc < acc_prev
%            add_if_true = ~add_if_true;
%            count_nochanges = 0;
%        else
%            if acc == acc_prev
%                count_nochanges = count_nochanges + 1;
%            else
%                count_nochanges = 0;
%            end
%        end
    end
   
%    if count_nochanges >= 4
%        isstable = 1;
%    end
   
   iter = iter + 1;
end

%% Add reactions to the model
function [model_new,rxnValid] = increaseModel(model1,model2,rxnValid)
    added = 0;
    model_new = model1;
    while added < 5 && ~isempty(rxnValid)
        id = datasample(rxnValid,1);
        fprintf('Reaction %d (Added %d)\n',id,added);
        rxnValid = setdiff(rxnValid,id);
        rxn_name = model2.rxnNames(id);
        rxn_name = rxn_name{1};
        met_names = model2.metNames(model2.S(:,id)~=0).';
        if isequal(sort(model2.metNames(model2.S(:,id)==1)),sort(model2.metNames(model2.S(:,id)==-1)))
            continue
        end
        S_new_sparse = full(model2.S(:,id));
        S_new = S_new_sparse(S_new_sparse~=0).';
        try
            [model_new,~] = addReaction(model_new,rxn_name,met_names,S_new);
        catch
        end
        added = length(model_new.rxnNames) - length(model1.rxnNames);
    end
end

%% Remove reactions to the model

function model = decreaseModel(model)
    for i=1:5
        model = removeRxns(model, datasample(model.rxnNames,1));
    end
end