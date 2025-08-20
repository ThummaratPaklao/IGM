function [modelIGM solIGM] = IGMRUN(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod, method)
% IGMRUN - Multi-Condition Integration of Gene Expression into Genome-Scale Metabolic Models 
%          with Optional L1 or L2 Norm Regularization
%
% USAGE:
%    [modelIGM, solIGM] = IGMRUN(model, uptakeRatesTable, geneexpressionTable)
%    [modelIGM, solIGM] = IGMRUN(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod, method)
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
%    method:                Method for solving the optimization problem:
%                               'IGM' → IGM without additional regularization
%                               'L1'  → IGM with L1 norm regularization
%                               'L2'  → IGM with L2 norm regularization
%
% OUTPUTS:
%    modelIGM:              IGM optimization model incorporating relative gene expression data.
%
%    solIGM:                Structure containing the optimization results:
%                               solIGM.v   → reaction flux distribution for each condition
%                               solIGM.x   → solution values for all variables in each condition
%                               solIGM.f   → optimal objective value for each condition
%
%                           L1-norm regularization results (if method = 'L1'):
%                               solIGM.vL1 → reaction flux distribution for each condition using L1 norm
%                               solIGM.xL1 → solution values for all variables in each condition using L1 norm
%                               solIGM.fL1 → optimal objective value for each condition using L1 norm
%
%                           L2-norm regularization results (if method = 'L2'):
%                               solIGM.vL2 → reaction flux distribution for each condition using L2 norm
%                               solIGM.xL2 → solution values for all variables in each condition using L2 norm
%                               solIGM.fL2 → optimal objective value for each condition using L2 norm
%
% EXAMPLES:
%    % Run with default settings (IGM without regularization):
%    [modelIGM, solIGM] = IGMRUN(model, uptakeRatesTable, geneexpressionTable)
%
%    % Run with L1-norm regularization:
%    [modelIGM, solIGM] = IGMRUN(model, uptakeRatesTable, geneexpressionTable, 1:3, 'maxmin', 'L1')
%
%    % Run with L2-norm regularization:
%    [modelIGM, solIGM] = IGMRUN(model, uptakeRatesTable, geneexpressionTable, [1:5,7:8], 'mean', 'L2')
%
%..Author: 
%    -Thummarat Paklao, 07/08/2025, Department of Mathematics and Computer Science, Faculty of Science, Chulalongkorn University, Thailand. 
%    -Apichat Suratanee, 07/08/2025, Department of Mathematics, Faculty of Applied Science, King Mongkut's University of Technology North Bangkok. 
%    -Kitiporn Plaimas, 07/08/2025, Department of Mathematics and Computer Science, Faculty of Science, Chulalongkorn University, Thailand.    

if (nargin < 4 || isempty(condition))
      condition=1:size(geneexpressionTable, 2) - 1;
end
if (nargin < 5 || isempty(method))
      method = 'IGM';
end
if (nargin < 6 || isempty(normalizemethod))
      normalizemethod = 'mean';
end

switch method
    case 'IGM' 
        [modelIGM solIGM] = IGM_function(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod);

    case 'L1'
        [modelIGM solIGM] = IGML1(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod);

    case 'L2'         
        [modelIGM solIGM] = IGML2(model, uptakeRatesTable, geneexpressionTable, condition, normalizemethod);

    otherwise
            error('Unknown method: %s', method);
end
    solIGM.status = 'complete';
end