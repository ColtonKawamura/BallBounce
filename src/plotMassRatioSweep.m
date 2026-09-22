function plotMassRatioSweep(normalizeByPredPeak)
% plotMassRatioSweep: Task 1 bar-theory baseline
%
% Optional argument:
%   normalizeByPredPeak (logical, default false)
%     false: x-axis is N, red stars at predicted peaks are shown
%     true:  x-axis is N / N_pred, NO red stars on main axes,
%            but a bottom-left INSET shows non-normalized e(N) with stars.

    if nargin < 1
        normalizeByPredPeak = false;
    end

    %% fixed parameters (same as plotSingle.m / plotPaper.m)
    NArr           = 3:30;
    vecMassHat     = [10, 12, 16];   % ball-to-chain mass ratios to sweep
    scalSpringHat  = 4;              % k_b/k_c
    scalDampHat    = 0.0025;         % gamma_hat
    scalVImpactHat = 0.1;            % v_hat
    scalGravityHat = 0;              % g_hat

    % marker-size scaling based on mass ratio
    baseMarkerSize = 6;                      % size for the smallest mass ratio
    minMassHat     = min(vecMassHat);

    %% figure + explicit MAIN AXES
    figure;
    mainAx = axes('Units', 'normalized', ...
                  'Position', [0.12 0.15 0.80 0.78]);  % big main panel
    hold(mainAx, 'on');
    set(mainAx, 'Color', 'none');   % transparent so inset can sit on top

    numMass   = numel(vecMassHat);
    matNPeak  = nan(numMass, 1);
    matNDelta = nan(numMass, 1);

    % store e(N) and predicted-peak info for inset
    eAll           = nan(numMass, numel(NArr));
    vecPredNPeak   = nan(numMass, 1);
    vecIdxPredPeak = nan(numMass, 1);

    %% MAIN AXES: sweep mass ratio, then chain length
    for idxMass = 1:numMass
        scalMassHat = vecMassHat(idxMass);

        ratios = nan(size(NArr));
        for idxN = 1:numel(NArr)
            [ratios(idxN), ~, ~] = ...
                sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
                      scalVImpactHat = scalVImpactHat, scalGravityHat = scalGravityHat);
        end

        % restitution e(N): e^2 = KE ratio
        e = sqrt(max(ratios, 0));
        eAll(idxMass, :) = e;  % store for inset

        % measured peak: first interior local maximum (skip N = 3 boundary)
        idxPeak = NaN;
        for i = 2:numel(NArr)-1
            if e(i) > e(i-1) && e(i) > e(i+1)
                idxPeak = i;
                break;
            end
        end
        if ~isnan(idxPeak)
            matNPeak(idxMass)  = NArr(idxPeak);
            matNDelta(idxMass) = NArr(idxPeak) - 1;
        end

        % predicted peak location from bar theory: N_peak - 1 ~ 0.6 * mHat
        predNPeak = 1 + 0.6 * scalMassHat;
        [~, idxPredPeak] = min(abs(NArr - predNPeak));  % nearest integer N
        vecPredNPeak(idxMass)   = predNPeak;
        vecIdxPredPeak(idxMass) = idxPredPeak;

        % choose x-values: either raw N, or normalized N / N_pred
        if normalizeByPredPeak
            xVals = NArr / predNPeak;
            % no star on main axes in normalized plot
        else
            xVals = NArr;
            xStar = NArr(idxPredPeak);
        end

        % marker size for this mass ratio (larger mass -> much larger squares)
        scaleFactor     = scalMassHat / minMassHat;
        markerSizeCurve = baseMarkerSize * (scaleFactor^2);   % quadratic scaling

        % red line with square markers (MAIN AXES)
        semilogx(mainAx, xVals, e, ...
            'LineWidth', 1.5, ...
            'Color', [0.65 0.00 0.30], ...
            'Marker', 's', ...
            'MarkerSize', markerSizeCurve, ...
            'LineStyle', '-', ...
            'DisplayName', sprintf('$\\hat{m} = %.0f$', scalMassHat));

        % predicted peak star ONLY in non-normalized plot
        if ~normalizeByPredPeak && ~isnan(idxPredPeak)
            semilogx(mainAx, xStar, e(idxPredPeak), 'r*', ...
                'MarkerSize', 12, 'LineWidth', 1.5);
        end
    end

    %% main-axis labels, legend, title
    axes(mainAx);  % set current axes explicitly

    if normalizeByPredPeak
        xlabel('$N / N_{\mathrm{pred}}$', 'Interpreter', 'latex', 'FontSize', 20);
    else
        xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
    end
    ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
    legend(mainAx, 'show', 'Location', 'best', ...
           'Interpreter', 'latex', 'FontSize', 13);
    grid(mainAx, 'on');
    set(mainAx, 'XScale', 'log');
    box(mainAx, 'on');
    title(mainAx, ...
          sprintf('$\\hat{g}=0,\\ \\hat{k}=%.1f,\\ \\hat{\\gamma}=%.4f,\\ \\hat{v}=%.1f$', ...
                  scalSpringHat, scalDampHat, scalVImpactHat), ...
          'Interpreter', 'latex', 'FontSize', 16);

    %% INSET: non-normalized e(N) with stars, only when normalizeByPredPeak == true
    if normalizeByPredPeak
        % bottom-left inset: [left bottom width height] in figure-normalized units
        insetAx = axes('Units', 'normalized', ...
                       'Position', [0.18 0.22 0.30 0.30]);  % small panel bottom-left
        hold(insetAx, 'on');
        box(insetAx, 'on');
        set(insetAx, 'Color', 'none');  % transparent background

        % plot non-normalized curves on inset from stored eAll
        for idxMass = 1:numMass
            scalMassHat = vecMassHat(idxMass);
            e           = eAll(idxMass, :);
            idxPredPeak = vecIdxPredPeak(idxMass);

            % marker size for this mass ratio (same scaling as main axes)
            scaleFactor     = scalMassHat / minMassHat;
            markerSizeCurve = baseMarkerSize * (scaleFactor^2);

            % non-normalized curve: x = NArr
            semilogx(insetAx, NArr, e, ...
                'LineWidth', 1.5, ...
                'Color', [0.65 0.00 0.30], ...
                'Marker', 's', ...
                'MarkerSize', markerSizeCurve, ...
                'LineStyle', '-');

            % red star at predicted peak on inset
            if ~isnan(idxPredPeak)
                xStarInset = NArr(idxPredPeak);
                semilogx(insetAx, xStarInset, e(idxPredPeak), 'r*', ...
                    'MarkerSize', 10, 'LineWidth', 1.5);
            end
        end

        % log x-scale and N-range
        set(insetAx, 'XScale', 'log');
        xlim(insetAx, [NArr(1), NArr(end)]);

        % auto y-limits from stored eAll
        allE = eAll(:);
        allE = allE(~isnan(allE));
        if ~isempty(allE)
            ymin = min(allE);
            ymax = max(allE);
            pad  = 0.02 * (ymax - ymin);
            ylim(insetAx, [ymin - pad, ymax + pad]);
        end

        % simplify inset (no labels, smaller font)
        set(insetAx, 'XTickLabel', [], 'YTickLabel', []);
        set(insetAx, 'FontSize', 10);
        grid(insetAx, 'on');
        box(insetAx, 'on');

        % make sure inset is drawn on top of main axes
        uistack(insetAx, 'top');

        % restore current axes to main (for any theming)
        axes(mainAx);
    end

    %% peak table: mHat, N_peak, N_peak-1, 0.6*mHat
    fprintf('\n=== Mass-ratio sweep: restitution peak (bar-theory baseline) ===\n');
    fprintf('Parameters: k_hat=%.2f, gamma_hat=%.4f, v_hat=%.2f, g_hat=%.1f, N=%d:%d\n', ...
            scalSpringHat, scalDampHat, scalVImpactHat, scalGravityHat, ...
            NArr(1), NArr(end));
    fprintf('Peak = first interior local maximum of e(N) (N=%d boundary excluded).\n', ...
            NArr(1));
    fprintf('%8s %10s %12s %12s\n', 'mHat', 'N_peak', 'N_peak-1', '0.6*mHat');
    for idxMass = 1:numMass
        fprintf('%8.0f %10.0f %12.0f %12.2f\n', ...
            vecMassHat(idxMass), matNPeak(idxMass), matNDelta(idxMass), ...
            0.6 * vecMassHat(idxMass));
    end
    fprintf('===\n\n');

    %% save figure into <repo>/figures/ (this file lives in <repo>/src/)
    thisDir = fileparts(mfilename('fullpath'));
    figDir  = fullfile(fileparts(thisDir), 'figures');
    figPath = fullfile(figDir, 'massRatioSweep.png');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
    print(figPath, '-dpng', '-r150');
    fprintf('[plotMassRatioSweep] Saved figure to %s\n', figPath);
end

