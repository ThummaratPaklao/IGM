function [modelIGM, solIGM] = IGML1(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod)
% IGML1 - Multi-Condition Integration of Gene Expression into Genome-Scale Metabolic Models 
%         with L1-Norm Regularization
%
% USAGE:
%    [modelIGM, solIGM] = IGML1(model, uptakeRatesTable, geneexpressionTable)
%    [modelIGM, solIGM] = IGML1(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod)
%
% INPUTS:
%    model:                 COBRA model structure (genome-scale metabolic model).
%
%    uptakeRatesTable:      Table of uptake rates where each column represents a condition 
%                           and each row corresponds to a reaction flux.
%                           The reaction list must match the reactions in the metabolic model.
%
%    geneexpressionTable:   Table of gene expression profiles where each column represents a 
%                           condition and each row corresponds to a gene in the metabolic model.
%
% OPTIONAL INPUTS:
%    condition:             Row vector specifying the condition indices corresponding to 
%                           the columns in geneexpressionTable.
%                           Examples:
%                               1:3       → conditions 1 to 3
%                               [1:5,7:8] → conditions 1 to 8 excluding condition 6
%                           Default: all conditions.
%
%    normalizemethod:       Method for normalizing gene expression data for each gene (row) 
%                           to relative values. Options:
%                               'max'    → use the maximum value in each row as the reference
%                               'maxmin' → use min-max scaling
%                               'mean'   → use the mean value in each row as the reference midpoint
%                           Default: 'mean'
%
% OUTPUTS:
%    modelIGM:              IGM optimization model incorporating relative gene expression data.
%
%    solIGM:                Structure containing the optimization results:
%                               solIGM.v   → reaction flux distribution for each condition
%                               solIGM.x   → solution values for all variables in each condition
%                               solIGM.f   → optimal objective value for each condition
%
%                           L1-norm regularization results:
%                               solIGM.vL1 → reaction flux distribution for each condition using L1 norm
%                               solIGM.xL1 → solution values for all variables in each condition using L1 norm
%                               solIGM.fL1 → optimal objective value for each condition using L1 norm
%
% Authors: 
%    -Thummarat Paklao, 09/08/2025, Department of Mathematics and Computer Science, Faculty of Science, Chulalongkorn University, Thailand. 
%    -Apichat Suratanee, 09/08/2025, Department of Mathematics, Faculty of Applied Science, King Mongkut's University of Technology North Bangkok. 
%    -Kitiporn Plaimas, 09/08/2025, Department of Mathematics and Computer Science, Faculty of Science, Chulalongkorn University, Thailand.    

if (nargin < 5 || isempty(condition))
      condition = 1 : size(geneexpressionTable,2) - 1;
end

if (nargin < 6 || isempty(normalizemethod))
      normalizemethod = 'mean';
end


% Prepare Data
exp = table2array(geneexpressionTable(:, 2:size(geneexpressionTable, 2)));    % extract gene expression value
geneName = table2cell(geneexpressionTable(:, 1));                            % extract gene name from gene expression table
uptakeRxnList = table2cell(uptakeRatesTable(:, 1));                          % extract lis of reactions in uptake reaction table

% remove missing value
[exp,TF] = rmmissing(exp, 1);
geneName = geneName(~TF, 1);

% Creation of Template Metabolic Model
model.lb(model.lb < 0) = -1000;
model.lb(model.lb >= 0) = 0;
model.ub(model.ub > 0) = 1000;
model.ub(model.ub <= 0) = 0;

% Perform FBA and FVA for each condition
minFlux1 = table();
maxFlux1 = table();
solbiomass = zeros(length(condition), 1);
for ch = condition
    % integrate uptake rates into metabolic model
    [rate, TF1] = rmmissing(uptakeRatesTable{:,ch+1}, 1);
    id_rxn = findRxnIDs(model, uptakeRxnList(~TF1));
    modelFBA = changeRxnBounds(model, model.rxns(id_rxn), rate, 'l');
    
    % convert FBA model to irreversible model
    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(modelFBA);
    modelIrrev.grRules = modelIrrev.rules;

    % solve FBA problem
    solFBA = optimizeCbModel(modelIrrev);
    solbiomass(ch) = solFBA.f;

    % solve FVA programming
    optPercentage = 99;
    [minFlux, maxFlux] = fluxVariability(modelIrrev, optPercentage, 'max', modelIrrev.rxns, 0, true);
    minFlux1(:, ch) = table(minFlux);
    maxFlux1(:, ch) = table(maxFlux);
end

% construct model and calculate flux distribution in each condition
solIGM.x = table();
solIGM.v = table(model.rxns);
solIGM.f = zeros(size(condition, 2), 1);
solIGM.vL1 = table(model.rxns);;
n = 0;
for ch = condition
    n = n + 1;
    % integrate uptake rates to metabolic model
    [rate, TF1] = rmmissing(uptakeRatesTable{:, condition(ch)+1}, 1);
    id_rxn = findRxnIDs(model, uptakeRxnList(~TF1));
    model4 = changeRxnBounds(model, model.rxns(id_rxn), rate, 'l');

    % convert to irreversible model
    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model4);
    
    %Set Matrix for GPR constraints
    M = 1;
    k = 1;
    p = 1;
    nn = 1;
    qi = sparse(2, 2);                         % new variables for GPR 
    yi = sparse(2, 2);                         % new binary variables for GPR
    G = sparse(2, size(geneName, 1));          % Gene matrix which correspond to variable g
    sumYi = sparse(2, 2);                      % Coefficint Matrix of variable yi
    BsumYi = sparse(2, 1);                     % RHS of matrix sumYi 
    gv = {};
    for i = 1 : length(modelIrrev.rxns)
        % extract GPR rule for each reaction
        rulegene = modelIrrev.rules{i};
        nrule = 0;
        nnrule = 0;
        gg=0;
        if ~isempty(rulegene)                %check empty GPR rule
            newrule = split(rulegene, "|");  % split GPR rule between 'OR'
            nrule = size(newrule, 1);        % number of splited rule between 'OR'
            gene_or = [];
            for j = 1 : nrule
                newrule1 = newrule{j};
                if ~ismember('&', newrule1)                                  % if each splitted rule do not contain 'AND'
                    gene_or(j) = string(regexp(newrule1, '\d*', 'Match'));    % extract number in gene name
                else                                                         % if each splitted rule contain 'AND'
                    newrule_and = split(newrule1, "&");                      % split GPR rule between 'AND'
                    nnrule = length(newrule_and);                            % number of splited rule between 'AND'
                    gene_and = [];
                    for m = 1 : nnrule
                        newrule_and1 = newrule_and{m};
                        gene_and(m) = string(regexp(newrule_and1, '\d*', 'Match'));
                        gg = find(string(geneName) == string(modelIrrev.genes(gene_and(m))));
                        G(k, gg) = -1;
                        G(nnrule + k, gg) = 1;
                        qi(k, nn) = 1;
                        qi(nnrule + k, nn) = -1;
                        yi(nnrule + k, p) = -M;
                        sumYi(nn, p) = 1;
                        p = p + 1;
                        k = k + 1;
                    end
                    gene_or(j) = 10^7 + nn;
                    k = k + nnrule;
                    b_sumYi(nn) = nnrule - 1;
                    nn = nn + 1;
                end
            end 
            if nrule > 1
                for c = 1 : nrule
                    qi(k, nn) = -1;
                    qi(nrule + k, nn) = 1;
                    if gene_or(c) > 10^7
                        qi(k, gene_or(c) - 10^7) = 1;
                        qi(nrule + k, gene_or(c) - 10^7) = -1;
                    else
                        gg = find(string(geneName) == string(modelIrrev.genes(gene_or(c))));
                        G(k, gg) = 1;
                        G(nrule + k, gg) = -1;
                    end
                    yi(nrule + k, p) = -M;
                    sumYi(nn, p) = 1;
                    p = p + 1;
                    k = k + 1;
                end
                gv{i} = 'x' + string(nn);
                BsumYi(nn) = nrule - 1;
                nn = nn + 1;
                k = k + nrule;
            elseif nrule == 1 && nnrule < 1 
                gv{i} = 'g' + string(regexp(newrule1, '\d*', 'Match'));
            elseif nrule == 1 && nnrule > 1
                gv{i} = 'x' + string(nn - 1);
            end
        else
            gv{i} = [];
        end 
    end

    M62 = sparse(G);
    M63 = sparse(qi);
    M66 = sparse(yi);

    M11 = sparse(modelIrrev.S);
    M21 = sparse(eye(length(modelIrrev.rxns)));
    M27 = sparse(eye(length(modelIrrev.rxns)));
    M41 = sparse(2, size(modelIrrev.rxns, 1));
    M42 = sparse(2, length(geneName));
    M43 = sparse(2, size(qi, 2));
    M45 = sparse(2, 2*length(modelIrrev.rxns));
    b4 = sparse(2, 1);
    gg = 0;
    k = 1;
    for i = 1 : length(modelIrrev.rxns)
        M27(i, i) = -modelIrrev.ub(i);
        if ~isempty(gv{i})
            a = gv{i}{1};   
            if max(table2array(maxFlux1(i, :))) - min(table2array(minFlux1(i, :))) > 0.0001 && a(1) == 'g'
                M41(k, i) = 1/(max(table2array(maxFlux1(i, :))) - min(table2array(minFlux1(i, :))));
                gg = find(string(geneName) == string(modelIrrev.genes(str2num(replace(gv{i}{1}, 'g', '')))));
                M42(k, gg) = -1;
                M43(k, :) = 0;
                M45(k, i) = 1;
                M45(k, length(modelIrrev.rxns) + i) = -1;
                b4(k, 1) = min(table2array(minFlux1(i, :)))/(max(table2array(maxFlux1(i, :))) - min(table2array(minFlux1(i, :))));

            elseif max(table2array(maxFlux1(i, :))) - min(table2array(minFlux1(i, :))) > 0.0001 && a(1) == 'x'
                M41(k, i) = 1/(max(table2array(maxFlux1(i, :))) - min(table2array(minFlux1(i, :))));
                M43(k, str2num(replace(gv{i}{1}, 'x', ''))) = -1;
                M42(k, :) = 0;
                M45(k, i) = 1;
                M45(k, length(modelIrrev.rxns) + i) = -1;
                b4(k, 1) = min(table2array(minFlux1(i, :)))/(max(table2array(maxFlux1(i, :))) - min(table2array(minFlux1(i, :))));
            end
            k = k + 1;
        end
    end


    M37 = sparse(2,size(modelIrrev.rxns, 1));
    k = 1;
    for i = 1:length(model.rxns)
        if length(rev2irrev{i, 1}) == 2
            M37(k, rev2irrev{i, 1}(1)) = 1;
            M37(k, rev2irrev{i, 1}(2)) = 1;
            k = k + 1;
        end
    end


    if strcmp(normalizemethod, 'max')
        k = 1;
        b5 = sparse(2, 1);
        M52 = sparse(2, length(geneName));
        M54 = sparse(2, 2);
        for i = 1 : length(geneName(:, 1)) - 1
            if max(exp(i, :)) - min(exp(i, :)) > 0.1 
                M52(k, k) = 1;
                M54(k, k) = 1;
                M54(k, length(geneName(:, 1)) - 1 + k) = -1;
                b5(k, 1) = (exp(i, ch))/ (max(exp(i, :)));
                k = k + 1;
            end
        end
    end

    if strcmp(normalizemethod, 'maxmin')
        k = 1;
        b5 = sparse(2, 1);
        M2 = sparse(2, length(geneName));
        M54 = sparse(2, 2);
        for i = 1:length(geneName(:, 1)) - 1
            if max(exp(i, :)) - min(exp(i, :)) > 0.1 
                M52(k, k) = 1;
                M54(k, k) = 1;
                M54(k, length(geneName(:, 1)) - 1 + k) = -1;
                b5(k, 1) = (exp(i, ch) - min(exp(i, :)))/ (max(exp(i, :)) - min(exp(i, :)));
                k = k + 1;
            end
        end
    end


    if strcmp(normalizemethod, 'mean')
        k = 1;
        b5 = sparse(2, 1);
        M52 = sparse(2, length(geneName));
        M54 = sparse(2, 2);
        for i = 1 : length(geneName(:, 1)) - 1
            if max(exp(i, :)) - min(exp(i, :)) > 0.1 
                meang = mean(exp(i, :));
                M52(k, k) = 1;
                M54(k, k) = 1;
                M54(k, length(geneName(:, 1)) - 1 + k) = -1;
                if exp(i, ch) >= meang
                    b5(k, 1) = (((exp(i, ch) - meang)/ (max(exp(i, :)) - meang)) + 1)/ 2;
                else 
                    b5(k, 1) = ((meang - exp(i, ch))/ (meang - min(exp(i, :))))/ 2;
                end    
                k = k + 1;
            end
        end
    end

    % Prepare zero Matrix for Optimization
    M12 = sparse(size(M11, 1), length(geneName));
    M13 = sparse(size(M11, 1), size(qi, 2));
    M14 = sparse(size(M11, 1), size(M54, 2));
    M15 = sparse(size(M11, 1), size(M45, 2));
    M16 = sparse(size(M11, 1), size(yi, 2));
    M17 = sparse(size(M11, 1), length(modelIrrev.rxns));

    M22 = sparse(size(M21, 1), length(geneName));
    M23 = sparse(size(M21, 1), size(qi, 2));
    M24 = sparse(size(M21, 1), size(M54, 2));
    M25 = sparse(size(M21, 1), size(M45, 2));
    M26 = sparse(size(M21, 1), size(yi, 2));

    M31 = sparse(size(M37, 1), size(M11, 2));
    M32 = sparse(size(M37, 1), size(M12, 2));
    M33 = sparse(size(M37, 1), size(M13, 2));
    M34 = sparse(size(M37, 1), size(M14, 2));
    M35 = sparse(size(M37, 1), size(M15, 2));
    M36 = sparse(size(M37, 1), size(M16, 2));

    M44 = sparse(size(M41, 1), size(M14, 2));
    M46 = sparse(size(M41, 1), size(M16, 2));
    M47 = sparse(size(M41, 1), size(M17, 2));

    M51 = sparse(size(M52, 1), size(M11, 2));
    M53 = sparse(size(M52, 1), size(M13, 2));
    M55 = sparse(size(M52, 1), size(M15, 2));
    M56 = sparse(size(M52, 1), size(M16, 2));
    M57 = sparse(size(M52, 1),size(M17, 2));

    M61 = sparse(size(G, 1), size(M11, 2));
    M64 = sparse(size(M62, 1), size(M14, 2));
    M65 = sparse(size(M62, 1), size(M15, 2));
    M67 = sparse(size(M62, 1), size(M17, 2));

    M76 = sparse(sumYi);
    M71 = sparse(size(M76, 1), size(M11, 2));
    M72 = sparse(size(M76, 1), size(M12, 2));
    M73 = sparse(size(M76, 1), size(M13, 2));
    M74 = sparse(size(M76, 1), size(M14, 2));
    M75 = sparse(size(M76, 1), size(M15, 2));
    M77 = sparse(size(M76, 1), size(M17, 2));

    b1 = sparse(zeros(size(M11, 1), 1));
    b2 = sparse(zeros(size(M21, 1), 1));
    b3 = sparse(ones(size(M31, 1), 1));
    b4 = sparse(b4);
    b5 = sparse(b5);
    b6 = sparse(zeros(size(M61, 1), 1));
    b7 = sparse(BsumYi);
    
    % Matrix of Optimization Problem
    Aeq = sparse([M11, M12, M13, M14, M15, M16, M17; M31, M32, M33, M34, M35, M36, M37; M41, M42, M43, M44, M45, M46, M47; M51, M52, M53, M54, M55, M56, M57; M71, M72, M73, M74, M75, M76, M77]);
    beq = sparse([b1; b3; b4; b5; b7]);
    A = sparse([M21, M22, M23, M24, M25, M26, M27; M61, M62, M63, M64, M65, M66, M67]);
    b = sparse([b2; b6]);
    f = sparse([-(((size(M14, 2) + size(M15, 2)))/ (solbiomass(ch)))*modelIrrev.c; zeros(size(M12, 2), 1); zeros(size(M13, 2), 1); ones(size(M14, 2), 1); ones(size(M15, 2), 1); zeros(size(M16, 2), 1); zeros(size(M17, 2), 1)]);
    lb = [modelIrrev.lb; zeros(size(M12, 2) + size(M13, 2) + size(M14, 2) + size(M15, 2) + size(M16, 2) + size(M17, 2), 1)];
    ub = [modelIrrev.ub; ones(size(M12, 2), 1); ones(size(M13, 2), 1); inf(size(M14, 2), 1); inf(size(M15, 2), 1); ones(size(M16, 2), 1); ones(size(M17, 2), 1)];


    intcon = (length(modelIrrev.c) + size(M12, 2) + size(M13, 2) + size(M14, 2) + size(M15, 2) + 1) : (length(modelIrrev.c) + size(M12, 2) + size(M13, 2) + size(M14, 2) + size(M15, 2) + size(M16, 2) + size(M17, 2));
    
    % Model IGM construction
    modelIGM = struct();
    modelIGM.A = sparse([Aeq; A]);
    modelIGM.obj = full(f);
    modelIGM.rhs = full([beq; b]);
    modelIGM.sense = [char('='*ones(size(Aeq, 1), 1)); char('<'*ones(size(A, 1), 1))];
    modelIGM.vtype = [char('C'*ones(intcon(1) - 1, 1)); char('B'*ones(max(intcon) - intcon(1) + 1, 1)); char('C'*ones(size(f, 1) - max(intcon), 1))];
    modelIGM.modelsense = 'min';
    modelIGM.lb = lb;
    modelIGM.ub = ub; 

    params.outputflag = 1;

    result = gurobi(modelIGM, params);

    sol_flux=[];
    for i = 1 : length(model.rxns)
        if length(rev2irrev{i, 1}) == 1
            sol_flux(i, 1) = result.x(i, 1);
        else
            sol_flux(i, 1) = (result.x(rev2irrev{i, 1}(1, 1)) - result.x(rev2irrev{i, 1}(1, 2)));
        end
    end
    solIGM.v(:, n + 1) = table(sol_flux);
    solIGM.x(:, n) = table(result.x);
    solIGM.f(n, 1) = result.objval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimization with L1 norm
    f1 = sparse(zeros(size(f)));
    A=sparse([M21, M22, M23, M24, M25, M26, M27; M71, M72, M73, M74, M75, M76, M77; f']);
    b=sparse([b2;b7; solIGM.f(ch, 1) + 0.001*abs(solIGM.f(ch, 1))]);
    for i = 1 : length(modelIrrev.rxns)
        f1(i) = 1;
    end
    modelL1 = struct();
    modelL1.A = sparse([Aeq; A]);
    modelL1.obj = full(f1);
    modelL1.rhs = full([beq; b]);
    modelL1.sense = [char('='*ones(size(Aeq, 1), 1)); char('<'*ones(size(A, 1), 1))];
    modelL1.vtype = [char('C'*ones(intcon(1)-1, 1)); char('B'*ones(max(intcon) - intcon(1) + 1, 1)); char('C'*ones(size(f, 1) - max(intcon), 1))];
    modelL1.modelsense = 'min';
% modelg1.Q=sparse(Q);
    modelL1.lb=lb;
    modelL1.ub=ub; 

    params.outputflag = 1;

    result = gurobi(modelL1, params);

    sol_flux=[];
    for i = 1 : length(model.rxns)
        if length(rev2irrev{i,1}) == 1
            sol_flux(i, 1) = result.x(i, 1);
        else
            sol_flux(i, 1) = (result.x(rev2irrev{i, 1}(1, 1)) - result.x(rev2irrev{i, 1}(1, 2)));
        end
    end
    solIGM.vL1(:, n + 1) = table(sol_flux);
    solIGM.xL1(:, n) = table(result.x);
    solIGM.fL1(n, 1) = result.objval;
end
end
