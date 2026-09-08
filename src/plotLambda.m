% fig4_sim_mastercurve.m
% Plot either e or \tilde{e} vs. Lambda_measured, analogous to Fig. 4.

clear; clc;

% --- choose what to plot on the y-axis: "e" or "etilde" ---
% scalPlotMode = "e";   % options: "e" or "etilde"
scalPlotMode = "etilde";   % options: "e" or "etilde"

% --- parameters ---
NArr           = 3:30;
vecMassHat     = [20];      % ball-to-chain mass ratios
vecVImpactHat  = [0.1];     % impact velocities
scalGravityHat = 0;

% --- spring–damping pairs: each row is [k_hat, gamma_hat] ---
matSpringDampHat = [ ...
    % 0.7,  0.03;    % set 1
    1.0,  0.006;     % set 2
    % 2.5,  0.0025;  % set 3
    % 3.5,  0.0001;  % set 4
];
% ------------------------------------------------------------

% --- plotting style ---
scalDottedBy = "mass";    % "mass" or "v"

figure; hold on;

numPairs = size(matSpringDampHat, 1);
t        = linspace(1, 0, numPairs)';
colors   = flip([t, zeros(numPairs,1), 1-t]);

% marker sizes per damping (larger damping -> larger markers)
vecDampAll    = matSpringDampHat(:,2);
scalDampMin   = min(vecDampAll);
scalDampMax   = max(vecDampAll);

if scalDampMax == scalDampMin
    % single pair (or identical dampings): use a fixed marker size
    vecMarkerSize = 8 * ones(numPairs, 1);
else
    vecMarkerSize = 4 + 8 * (vecDampAll - scalDampMin) / (scalDampMax - scalDampMin);  % 4–12
end

% line styles keyed by mass and by impact velocity
vecLineStyleMass = repmat({"-"}, size(vecMassHat));     % default solid
vecLineStyleMass(2:end) = {":"};

vecLineStyleV    = repmat({"-"}, size(vecVImpactHat));
vecLineStyleV(2:end)    = {":"};

legendEntries = cell(numPairs, 1);

% ===== LOOP OVER SPRING–DAMPING PAIRS =====
for idxPair = 1:numPairs

    scalSpringHat = matSpringDampHat(idxPair, 1);
    scalDampHat   = matSpringDampHat(idxPair, 2);

    for idxMass = 1:length(vecMassHat)
        scalMassHat = vecMassHat(idxMass);

        for idxV = 1:length(vecVImpactHat)
            scalVImpactHat = vecVImpactHat(idxV);

            if scalDottedBy == "mass"
                scalLineStyle = vecLineStyleMass(idxMass);
            else
                scalLineStyle = vecLineStyleV(idxV);
            end

            ratios           = nan(size(NArr));
            lambdas_theory   = nan(size(NArr)); %#ok<NASGU>
            lambdas_measured = nan(size(NArr));

            % --- run sim1d over N ---
            for idxN = 1:length(NArr)
                [ratios(idxN), lambdas_theory(idxN), lambdas_measured(idxN)] = ...
                    sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
                          scalVImpactHat=scalVImpactHat, scalGravityHat=scalGravityHat);
            end

            % restitution e(N): e^2 = ratios
            e = sqrt(max(ratios, 0));

            % basic validity mask (finite e and Lambda)
            valid = ~isnan(e) & isfinite(e) & ...
                    ~isnan(lambdas_measured) & isfinite(lambdas_measured);

            if nnz(valid) < 4
                continue;
            end

            % --- optionally compute tilde{e} ---
            if scalPlotMode == "etilde"
                % use e_max over valid points
                e_valid = e(valid);
                e_max   = max(e_valid);

                % plateau: mean of last few valid points in N (tail of valid indices)
                validIdxAll  = find(valid);
                lastValidIdx = validIdxAll(end);
                firstTailIdx = validIdxAll(max(1, end-3));
                tailIdxN     = validIdxAll(validIdxAll >= firstTailIdx);  % up to lastValidIdx

                e_inf = mean(e(tailIdxN), 'omitnan');

                % guard against degenerate normalization
                if abs(e_max - e_inf) < 1e-8
                    % skip this curve if e_max ~ e_inf
                    continue;
                end

                % normalized restitution
                e_tilde = (e - e_inf) / (e_max - e_inf);

                y_raw = e_tilde;
            else
                % plain restitution
                y_raw = e;
            end

            % restrict to finite y values and positive Lambda
            goodIdx = valid & isfinite(y_raw) & ...
                      isfinite(lambdas_measured) & (lambdas_measured > 0);

            if nnz(goodIdx) < 2
                continue;
            end

            lambda_plot = lambdas_measured(goodIdx);
            y_plot_raw  = y_raw(goodIdx);

            % sort by Lambda for smooth curves
            [lambda_sorted, idxSort] = sort(lambda_plot);
            y_sorted = y_plot_raw(idxSort);

            % final filter after sorting
            posIdx = isfinite(lambda_sorted) & (lambda_sorted > 0) & isfinite(y_sorted);
            if nnz(posIdx) < 2
                continue;
            end
            lambda_sorted = lambda_sorted(posIdx);
            y_sorted      = y_sorted(posIdx);

            % --- plot: simple curve, no inset ---
            plot(lambda_sorted, y_sorted, ...
                'LineWidth', 1.5, ...
                'Color', colors(idxPair,:), ...
                'Marker', 'o', ...
                'MarkerSize', vecMarkerSize(idxPair), ...
                'LineStyle', scalLineStyle{1});

            legendEntries{idxPair} = sprintf('$\\hat{k}=%.2f,\\ \\hat{\\gamma}=%.3f$', ...
                                             scalSpringHat, scalDampHat);
        end
    end
end
% =============================================================

xlabel('$\\Lambda_\\text{measured}$', 'Interpreter', 'latex', 'FontSize', 20);

if scalPlotMode == "etilde"
    ylabel('$\\tilde{e} = (e - e_{\\infty}) / (e_{\\max} - e_{\\infty})$', ...
           'Interpreter', 'latex', 'FontSize', 18);
else
    ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
end

legend(legendEntries, ...
       'Location', 'best', ...
       'Interpreter', 'latex', ...
       'FontSize', 13, ...
       'NumColumns', 1);

grid on;
box on;
set(gca, 'XScale', 'linear', 'YScale', 'linear');

theme(gcf, 'light');

