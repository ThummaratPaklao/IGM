function plotcomparetwocond(solIGM, conditionIdx, filtervalue, top)
% PLOTCOMPARETWOCOND - Compare flux distributions between two conditions
% and identify top reaction fluxes change
%
%   plotcomparetwocond(solIGM, conditionIdx, filtervalue, top)
%
%   INPUTS:
%       solIGM        : solution structure containing flux values (table format)
%                       - first column = flux names
%                       - other columns = flux values across conditions
%       conditionIdx  : vector of 2 indices (relative to fluxTable columns, excluding first column)
%                       specifying which two conditions to compare
%       filtervalue   : threshold to filter out fluxes with values in [-filtervalue, +filtervalue]
%       top           : number of top up- and down-regulated fluxes to display
%
%   STEPS:
%       1. Extract flux names and values for the two selected conditions.
%       2. Remove fluxes that are zero across both conditions.
%       3. Filter fluxes within [-filtervalue, +filtervalue].
%       4. Compute relative flux change (log2 fold change).
%       5. Select top "up-regulated" and "down-regulated" fluxes.
%       6. Plot:
%          - Horizontal bar plot of relative flux changes.
%          - Scatter plot comparing the two conditions, highlighting top changes.
%
%   Example:
%       plotcomparetwocond(solIGM, [1 2], 1e-6, 15)

    %% Step 1: Extract flux data
    T = solIGM.v;
    fluxNames = T{:,1};                       % Flux names
    fluxValues = T{:, conditionIdx+1};        % Flux values for the two selected conditions

    %% Step 2: Remove all-zero flux rows
    rowsAllZero = all(fluxValues == 0, 2);
    fluxValues(rowsAllZero,:) = [];
    fluxNames(rowsAllZero,:) = [];

    %% Step 3: Filter fluxes within valid range
    rowsValid = all(fluxValues >= -filtervalue & fluxValues <= filtervalue, 2);
    fluxValues = fluxValues(rowsValid,:);
    fluxNames = fluxNames(rowsValid,:);

    %% Step 4: Compute relative flux change (log2 fold change)
    fluxValues1 = abs(fluxValues);
    epsilon = 0.1; % to avoid log(0)
    Var1 = fluxValues1(:,1);
    Var2 = fluxValues1(:,2);
    RelativeFluxChange = log2((Var2 + epsilon) ./ (Var1 + epsilon));

    % Build result table
    resultTable = table(fluxNames, fluxValues1(:,1), fluxValues1(:,2), RelativeFluxChange, ...
        'VariableNames', {'FluxName', 'Cond1', 'Cond2', 'RelativeFluxChange'});

    %% Step 5: Select top up- and down-regulated fluxes
    T_sorted = sortrows(resultTable, 'RelativeFluxChange');
    top_down = T_sorted(1:top, :);                  % most down-regulated
    top_up   = T_sorted(end-top+1:end, :);          % most up-regulated
    plot_T   = sortrows([top_down; top_up], 'RelativeFluxChange');

    fluxNames1 = plot_T.FluxName;
    relFluxChange = plot_T.RelativeFluxChange;

    %% Step 6a: Horizontal bar plot of relative flux changes
    figure();
    h = barh(relFluxChange, 'FaceColor', 'flat');

    % Color bars: green = up, red = down
    for i = 1:length(relFluxChange)
        if relFluxChange(i) >= 0
            h.CData(i,:) = [0 1 0];
        else
            h.CData(i,:) = [1 0 0];
        end
    end

    set(gca, 'ytick', 1:length(fluxNames1), 'yticklabel', string(fluxNames1));
    grid on;
    xlabel('Relative Flux Change (log2 fold change)');
    title('Top Up/Down-Regulated Fluxes');
    xlim([min(relFluxChange)*1.1, max(relFluxChange)*1.1]);

    %% Step 6b: Scatter plot of two conditions
    % Identify fluxes from top list
    specialList0 = plot_T.FluxName(plot_T.RelativeFluxChange < 0); % down
    specialList1 = plot_T.FluxName(plot_T.RelativeFluxChange > 0); % up
    highlightMask0 = ismember(fluxNames, specialList0);
    highlightMask1 = ismember(fluxNames, specialList1);

    % Log2 transform
    x = log2(abs(fluxValues(:,1))+1);
    y = log2(abs(fluxValues(:,2))+1);

    figure();
    % Default points (blue)
    scatter(x(~highlightMask0 & ~highlightMask1), y(~highlightMask0 & ~highlightMask1), ...
        20, 'b', 'filled'); 
    hold on;

    % Highlighted points (red = down, green = up)
    scatter(x(highlightMask0), y(highlightMask0), 200, 'r', 'filled');
    scatter(x(highlightMask1), y(highlightMask1), 200, 'g', 'filled');

    % Reference line (y = x)
    lims = [min([x; y]), max([x; y])];
    plot(lims, lims, 'k--', 'LineWidth', 2);

    % Label highlighted points
    for i = find(highlightMask0)' % down
        text(x(i)+0.05, y(i)+0.1, fluxNames{i}, 'FontSize', 10, 'Color', 'r', 'Interpreter', 'none');
    end
    for i = find(highlightMask1)' % up
        text(x(i)+0.05, y(i)+0.1, fluxNames{i}, 'FontSize', 10, 'Color', 'g', 'Interpreter', 'none');
    end

    xlabel(sprintf('Condition %d Flux (log2 scale)', conditionIdx(1)));
    ylabel(sprintf('Condition %d Flux (log2 scale)', conditionIdx(2)));
    title(sprintf('Condition %d vs Condition %d Flux Comparison', conditionIdx(1), conditionIdx(2)));
    axis equal;
    grid on;
    hold off;
end
