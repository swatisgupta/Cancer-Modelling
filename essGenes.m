function egenes = essGenes(Achiles, CCLE, celline, ge_threshold)
    
    id_A = find(~cellfun(@isempty,strfind(Achiles.cellines,celline))>0);
    genes_A = Achiles.genes(Achiles.Mat(:,id_A) < -0.5);
    
    id_C = find(~cellfun(@isempty,strfind(CCLE.cell_line,celline))>0);
    genes_C = CCLE.genes(CCLE.GE(:,id_C) > ge_threshold);
    
    egenes = intersect(genes_A,genes_C);
end