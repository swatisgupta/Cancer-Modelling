function egenes = essMOMA(model)

 biomass = model.biomass_rxn;
 ngenes = length(model.genes);
 egenes = zeros(ngenes,1);
 for i = 1 : ngenes
     new_model = knockdown(model,find(model.rxnGeneMat(:,i)>0));
     [sM,sO,~,~] = MOMA(model,new_model);
     if (sO.x(biomass)*0.1 > sM.x(biomass))
         egenes(i) = 1;
         %fprintf('%d\n',i);
     end
 end
end

function new_model = knockdown(model,reaction)
    new_model = model;
    for i = 1 : length(reaction)
        new_model.lb(reaction(i)) = 0;
        new_model.ub(reaction(i)) = 0;
    end
end