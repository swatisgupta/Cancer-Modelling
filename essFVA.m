function egenes = essFVA(model)

 biomass = model.biomass_rxn;
 ngenes = length(model.genes);
 egenes = zeros(ngenes,1);
 sO = optimizeCbModel(model);
 for i = 1 : ngenes
     x = ones(ngenes,1);
     x(i) = 0;
     r = find(~cellfun(@isempty,strfind(model.rules,sprintf('x(%d)',i)))>0);
     new_model = knockdown(model,r,x);
     sM = optimizeCbModel(new_model);
     if (sO.x(biomass)*0.1 > sM.x(biomass))
         egenes(i) = 1;
         %fprintf('%d\n',i);
     end
 end
 fprintf('FBA: found%d essential genes\n', length(egenes));
end

function new_model = knockdown(model,reaction)
    new_model = model;
    for i = 1 : length(reaction)
        if (eval(model.rules{reaction(i)})<1)
            new_model.lb(reaction(i)) = 0;
            new_model.ub(reaction(i)) = 0;
        end
    end
end