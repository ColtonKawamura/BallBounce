function plotMassRatioSweep(normalizeByPredPeak)
% plotMassRatioSweep: Task 1 bar-theory baseline
%
% Show how the restitution peak in e(N) shifts with the ball-to-chain mass
% ratio mHat = m_b/m_c, using the existing 1D DEM (sim1d.m).
%
% Sweeps the chain length NArr and the mass ratios vecMassHat at fixed
% spring ratio, damping, impact velocity and zero gravity (the same
% "good" parameters used by plotSingle.m / plotPaper.m). Restitution
% e(N) = sqrt(KE_after/KE_before) is computed from the kinetic-energy
% ratio returned by sim1d, exactly as in the existing plotting scripts.
%
% Peak definition: the first INTERIOR local maximum of e(N) (the first
% index i with 1 < i < numel(NArr) and e(i) > e(i-1) and e(i) > e(i+1),
% i.e. after the trivial low-N rise). The global maximum of e(N) sits at
% the N = 3 boundary for every mass ratio (a very short stack behaves
% like a rigid wall, so e is highest there), which would make N_peak-1 a
% constant and meaningless for the bar-theory comparison; the interior
% peak is the feature the theory addresses.
%
% For each mHat the peak location is printed in a small table
% (mHat, N_peak, N_peak-1, 0.6*mHat) so it can be compared against the
% bar-theory scaling N_peak - 1 ~ 0.6 * mHat.
%
% All three e(N) curves are drawn on a single figure with:
%   - red lines
%   - square markers on the lines
%   - square marker size strongly increasing with mass ratio
%   - red stars at the predicted peak positions ONLY on the non-normalized plot
% and saved to <repo>/figures/massRatioSweep.png (created if needed).
%
% Optional argument:
%   normalizeByPredPeak (logical, default false)
%     false: x-axis is N, red stars at predicted peaks are shown
%     true:  x-axis is N / N_pred, NO red stars are plotted

    if nargin < 1
        normalizeByPredPeak = false;
    end

    %% fixed parameters (same as plotSingle.m / plotPaper.m)
    NArr           = 3:30;
    vecMassHat     = [10, 12, 16];   % ball-to-chain mass ratios to sweep, [10, 12, 16] was good
    scalSpringHat  = 4;              % k_b/k_c [4] was good
    scalDampHat    = 0.0025;         % gamma_hat [.0025] was good
    scalVImpactHat = 0.1;            % v_hat
    scalGravityHat = 0;              % g_hat

    % marker-size scaling based on mass ratio
    baseMarkerSize = 6;                      % size for the smallest mass ratio
    minMassHat     = min(vecMassHat);

    figure; hold on;

    matNPeak  = nan(numel(vecMassHat), 1);
    matNDelta = nan(numel(vecMassHat), 1);

    %% sweep mass ratio, then chain length
    for idxMass = 1:numel(vecMassHat)
        scalMassHat = vecMassHat(idxMass);

        ratios = nan(size(NArr));
        for idxN = 1:numel(NArr)
            [ratios(idxN), ~, ~] = ...
                sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
                      scalVImpactHat = scalVImpactHat, scalGravityHat = scalGravityHat);
        end

        % restitution e(N): e^2 = KE ratio (same as existing plotting code)
        e = sqrt(max(ratios, 0));

        % measured peak: first interior local maximum (skip the N = 3 boundary)
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
        [~, idxPredPeak] = min(abs(NArr - predNPeak));  % nearest integer N in NArr

        % choose x-values: either raw N, or normalized N / N_pred
        if normalizeByPredPeak
            xVals = NArr / predNPeak;
            xStar = NArr(idxPredPeak) / predNPeak;  % should be ~1 (not used for star in normalized plot)
        else
            xVals = NArr;
            xStar = NArr(idxPredPeak);
        end

        % marker size for this mass ratio (larger mass -> much larger squares)
        scaleFactor     = scalMassHat / minMassHat;
        markerSizeCurve = baseMarkerSize * (scaleFactor^2);   % quadratic scaling for larger differences

        % red line with square markers
        semilogx(xVals, e, ...
            'LineWidth', 1.5, ...
            'Color', [0.65 0.00 0.30], ...
            'Marker', 's', ...
            'MarkerSize', markerSizeCurve, ...
            'LineStyle', '-', ...
            'DisplayName', sprintf('$\\hat{m} = %.0f$', scalMassHat));

        % mark the PREDICTED peak (red star at theoretical location) ONLY in non-normalized plot
        if ~normalizeByPredPeak && ~isnan(idxPredPeak)
            semilogx(xStar, e(idxPredPeak), 'r*', ...
                'MarkerSize', 12, 'LineWidth', 1.5);
        end
    end

    if normalizeByPredPeak
        xlabel('$N / N_{\mathrm{pred}}$', 'Interpreter', 'latex', 'FontSize', 20);
    else
        xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
    end
    ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
    legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    xscale(gca, 'log');
    box on;
    title(sprintf('$\\hat{g}=0,\\ \\hat{k}=%.1f,\\ \\hat{\\gamma}=%.4f,\\ \\hat{v}=%.1f$', ...
                  scalSpringHat, scalDampHat, scalVImpactHat), ...
          'Interpreter', 'latex', 'FontSize', 16);

    %% peak table: mHat, N_peak, N_peak-1, 0.6*mHat
    fprintf('\n=== Mass-ratio sweep: restitution peak (bar-theory baseline) ===\n');
    fprintf('Parameters: k_hat=%.2f, gamma_hat=%.4f, v_hat=%.2f, g_hat=%.1f, N=%d:%d\n', ...
            scalSpringHat, scalDampHat, scalVImpactHat, scalGravityHat, ...
            NArr(1), NArr(end));
    fprintf('Peak = first interior local maximum of e(N) (N=%d boundary excluded).\n', ...
            NArr(1));
    fprintf('%8s %10s %12s %12s\n', 'mHat', 'N_peak', 'N_peak-1', '0.6*mHat');
    for idxMass = 1:numel(vecMassHat)
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

