function [accuracy, essGM, essGM_id] = evaluateModel(model, essG)

    essGM_id = essFBA(model);
    essGM = model.genes_unique_names(unique(model.genes_unique_map(essGM_id>0)));
    
    TP = intersect(essGM,essG);
    FP = setdiff(essGM,essG);
    Neg = setdiff(model.genes_unique_names,essG);
    
    accuracy = length(TP)/length(essG) - length(FP)/length(Neg);

end