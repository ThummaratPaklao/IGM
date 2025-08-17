function scatterplotcompareflux(solIGM, conditionIdx, filtervalue)
% SCATTERPLOTCOMPAREFLUX - Generate pairwise scatter plots to compare flux values
%
%   scatterplotcompareflux(solIGM, conditionIdx, filtervalue)
%
%   INPUTS:
%    solIGM:                Structure containing the optimization results from IGM model:
%                               solIGM.v   → reaction flux distribution for each condition
%                               solIGM.x   → solution values for all variables in each condition
%                               solIGM.f   → optimal objective value for each condition
%       conditionIdx  : vector of indices (relative to fluxTable columns, excluding first column)
%                       specifying which conditions to compare
%       filtervalue   : threshold to filter out fluxes within range [-filtervalue, +filtervalue]
%
%   This function:
%       1. Extracts flux names and flux values for selected conditions.
%       2. Removes fluxes that are zero across all selected conditions.
%       3. Filters out fluxes with absolute values less than or equal to `filtervalue`.
%       4. Creates scatter plots for every pair of selected conditions.
%       5. Applies log2-transformation (with +1 offset) to handle small values and visualize scaling.
%       6. Adds a reference line (y = x) to indicate equality between two conditions.
%
%   Example:
%       scatterplotcompareflux(solIGM, [1 2 3], 100)

    if (nargin < 2 || isempty(conditionIdx))
        conditionIdx = 1 : size(solIGM.v,2)-1;
    end

    if (nargin < 3 || isempty(filtervalue))
        filtervalue = 100;
    end

    % Extract flux table from solIGM
    T = solIGM.v;

    % Get flux names (first column of the table)
    fluxNames = T{:,1};

    % Extract flux values for selected conditions
    % +1 adjustment because the first column contains flux names
    fluxValues = T{:, conditionIdx+1};

    % Remove rows where flux values are all zero across selected conditions
    rowsAllZero = all(fluxValues == 0, 2);
    fluxValues(rowsAllZero, :) = [];
    fluxNames(rowsAllZero) = [];

    % Filter out fluxes with absolute values within [-filtervalue, +filtervalue]
    rowsValid = all(fluxValues >= -filtervalue & fluxValues <= filtervalue, 2);
    fluxValues = fluxValues(rowsValid, :);
    fluxNames = fluxNames(rowsValid);

    % Number of selected conditions
    numConditions = length(conditionIdx);

    % Total number of pairwise scatter plots
    numPlots = nchoosek(numConditions, 2);

    % Arrange subplot grid size (square layout)
    nRows = ceil(sqrt(numPlots));
    nCols = ceil(numPlots / nRows);

    plotCount = 1;

    % Loop through all pairs of conditions
    for i = 1:numConditions-1
        for j = i+1:numConditions
            subplot(nRows, nCols, plotCount)

            % Scatter plot with log2(abs(value)+1) transformation
            scatter(log2(abs(fluxValues(:,i))+1), log2(abs(fluxValues(:,j))+1), ...
                5, 'filled')
            hold on

            % Reference line y = x
            lims = [min([log2(abs(fluxValues(:,i))+1); log2(abs(fluxValues(:,j))+1)]), ...
                    max([log2(abs(fluxValues(:,i))+1); log2(abs(fluxValues(:,j))+1)])];
            plot(lims, lims, 'r--', 'LineWidth', 1.5)

            % Set axis limits
            xlim(lims)
            ylim(lims)

            % Labels and title
            xlabel(['Condition ' num2str(conditionIdx(i))])
            ylabel(['Condition ' num2str(conditionIdx(j))])
            title(['Cond ' num2str(conditionIdx(i)) ' vs ' num2str(conditionIdx(j))])

            % Grid for better readability
            grid on
            hold off

            plotCount = plotCount + 1;
        end
    end
end
