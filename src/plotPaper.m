mainTwo = true;
% Full script: sweep over (spring, damping) pairs and plot e(N)

if mainTwo == true
    clear; clc;

    NArr           = 3:30;
    vecMassHat     = [20];      % ball-to-chain mass ratios
    vecVImpactHat  = [0.1];     % impact velocities
    scalGravityHat = 0;

    % --- spring–damping pairs: each row is [k_hat, gamma_hat] ---
    matSpringDampHat = [ ...
        0.7, 0.03;     % line 1
        1.0, 0.006;    % line 2
        2.5, 0.0025;   % line 3
        3.5, 0.0001;   % line 4
        % add more rows as needed: [k_hat, gamma_hat]
    ];
    % ------------------------------------------------------------

    % --- single switch: choose which parameter controls dotted vs solid ---
    scalDottedBy = "mass";    % set to "mass" or "v"
    % ----------------------------------------------------------------------

    figure; hold on;

    % colors per spring–damping pair
    numPairs = size(matSpringDampHat, 1);
    t        = linspace(1, 0, numPairs)';
    colors   = flip([t, zeros(numPairs,1), 1-t]);

    % marker sizes per damping (larger damping -> larger markers)
    vecDampAll    = matSpringDampHat(:,2);
    scalDampMin   = min(vecDampAll);
    scalDampMax   = max(vecDampAll);
    vecMarkerSize = 4 + 8 * (vecDampAll - scalDampMin) / (scalDampMax - scalDampMin);  % 4–12

    % line styles keyed by mass and by impact velocity
    vecLineStyleMass = repmat({"-"}, size(vecMassHat));     % default solid
    vecLineStyleMass(2:end) = {":"};                        % higher masses dotted

    vecLineStyleV    = repmat({"-"}, size(vecVImpactHat));  % default solid
    vecLineStyleV(2:end)    = {":"};                        % second/others dotted

    for idxPair = 1:numPairs
        scalSpringHat = matSpringDampHat(idxPair, 1);
        scalDampHat   = matSpringDampHat(idxPair, 2);

        for idxMass = 1:length(vecMassHat)
            scalMassHat = vecMassHat(idxMass);

            for idxV = 1:length(vecVImpactHat)
                scalVImpactHat = vecVImpactHat(idxV);

                % choose line style based on single switch
                if scalDottedBy == "mass"
                    scalLineStyle = vecLineStyleMass(idxMass);
                else
                    scalLineStyle = vecLineStyleV(idxV);
                end

                ratios           = nan(size(NArr));
                lambdas_theory   = nan(size(NArr));
                lambdas_measured = nan(size(NArr));

                for idxN = 1:length(NArr)
                    [ratios(idxN), lambdas_theory(idxN), lambdas_measured(idxN)] = ...
                        sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
                              scalVImpactHat=scalVImpactHat, scalGravityHat=scalGravityHat);
                end

                e = sqrt(max(ratios, 0));
                e(ratios <= .2) = NaN;

                % DATA LINES ONLY (tagged)
                semilogx(NArr, e, ...
                    'LineWidth', 1.5, ...
                    'Color', colors(idxPair,:), ...
                    'Marker', 'o', ...
                    'MarkerSize', vecMarkerSize(idxPair), ...
                    'LineStyle', scalLineStyle{1}, ...
                    'DisplayName', sprintf('$\\hat{k}=%.2f,\\ \\hat{\\gamma}=%.3f$', ...
                                           scalSpringHat, scalDampHat), ...
                    'Tag', 'dataLine');   % <-- tag so inset can select only data
            end
        end
    end

    xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);

    legend('show', ...
           'Location', 'best', ...
           'Interpreter', 'latex', ...
           'FontSize', 13, ...
           'NumColumns', 1);

    xlim([6, 25]);
    ylim([.76, .92]);
    grid on;
    xscale(gca, 'log');
    box on;

    % --- HARD-CODED VERTICAL DOTTED LINES AND ARROW ON MAIN AXES ---
    % vertical dotted lines at x = 10 and x = 12
    % xline(10, 'LineStyle', ':', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
    % xline(12, 'LineStyle', ':', 'Color', [0.2 0.2 0.8], 'LineWidth', 1.5);
    xline(10, 'LineStyle', ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);  % dark gray
    xline(10, 'LineStyle', ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);  % light gray

    % fixed y-position for the horizontal arrow line (fits in [0.76, 0.92])
    yArrow = 0.90;

    % horizontal line between x = 10 and x = 12
    % line([10 12], [yArrow yArrow], ...
    %      'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
    % arrow from x = 12 to x = 10 (pointing left)
    % arrow from x = 12 to x = 10 (pointing left), with a bigger head
    % --- custom left-pointing arrow from x = 12 to x = 10 ---
xTail  = 12;           % where shaft starts (right)
xTip   = 10;           % arrow tip (left)
yArrow = 0.90;         % already defined above

shaftEnd   = xTip + 0.3;   % where shaft meets head (leave 0.3 in N for head)
shaftColor = [0.3 0.3 0.3];

% shaft (horizontal line)
line([xTail shaftEnd], [yArrow yArrow], ...
     'Color', shaftColor, 'LineStyle', '-', 'LineWidth', 1.5);

% triangular head
headWidth  = 0.1;      % extent in x
headHeight = 0.003;    % extent in y

xHead = [shaftEnd, shaftEnd, xTip];
yHead = [yArrow - headHeight, yArrow + headHeight, yArrow];

patch(xHead, yHead, shaftColor, ...
      'EdgeColor', shaftColor);
% --- end custom arrow ---

    % center of the arrow text in log-x coordinates (geometric mean)
    xCenter = sqrt(10 * 12);

    % text with left-pointing arrow
    text(xCenter, yArrow + 0.01, ...
         '$N_c$ shift', ...
         'Interpreter', 'latex', ...
         'HorizontalAlignment', 'center', ...
         'Color', [0.2 0.2 0.2], ...
         'FontSize', 14);
    % --- END ANNOTATION BLOCK ---

    % --- inset zoom: N in [8,16], e in [0.8,0.92] ---
    mainAx  = gca;

    % position = [left bottom width height] in figure-normalized units
    insetAx = axes('Position', [0.58 0.61 0.28 0.28]);
    box(insetAx, 'on'); hold(insetAx, 'on');

    % copy ONLY data lines (tagged 'dataLine') from the main axes into the inset
    hLines = findobj(mainAx, 'Type', 'line', '-and', 'Tag', 'dataLine');
    copyobj(hLines, insetAx);

    % keep semilog x-scale and zoom into desired region
    set(insetAx, 'XScale', 'log');
    xlim(insetAx, [8 16]);
    ylim(insetAx, [0.8 0.92]);

    % optional: simplify inset (no labels, smaller font)
    set(insetAx, 'XTickLabel', [], 'YTickLabel', []);
    set(insetAx, 'FontSize', 10);
    grid on;
    box on;
    % -----------------------------------------------
    theme(gcf, 'light');
else
    clear; clc;

    NArr           = 3:30;
    vecMassHat     = [20];      % ball-to-chain mass ratios
    vecVImpactHat  = [0.1];     % impact velocities
    scalGravityHat = 0;

    % --- spring–damping pairs: each row is [k_hat, gamma_hat] ---
    matSpringDampHat = [ ...
        0.7, 0.03;     % line 1
        1.0, 0.006;    % line 2
        2.5, 0.0025;   % line 3
        3.5, 0.0001;   % line 4
        % add more rows as needed: [k_hat, gamma_hat]
    ];
    % ------------------------------------------------------------

    % --- single switch: choose which parameter controls dotted vs solid ---
    scalDottedBy = "mass";    % set to "mass" or "v"
    % ----------------------------------------------------------------------

    figure; hold on;

    % colors per spring–damping pair
    numPairs = size(matSpringDampHat, 1);
    t        = linspace(1, 0, numPairs)';
    colors   = flip([t, zeros(numPairs,1), 1-t]);

    % marker sizes per damping (larger damping -> larger markers)
    vecDampAll    = matSpringDampHat(:,2);
    scalDampMin   = min(vecDampAll);
    scalDampMax   = max(vecDampAll);
    vecMarkerSize = 4 + 8 * (vecDampAll - scalDampMin) / (scalDampMax - scalDampMin);  % 4–12

    % line styles keyed by mass and by impact velocity
    vecLineStyleMass = repmat({"-"}, size(vecMassHat));     % default solid
    vecLineStyleMass(2:end) = {":"};                        % higher masses dotted

    vecLineStyleV    = repmat({"-"}, size(vecVImpactHat));  % default solid
    vecLineStyleV(2:end)    = {":"};                        % second/others dotted

    for idxPair = 1:numPairs
        scalSpringHat = matSpringDampHat(idxPair, 1);
        scalDampHat   = matSpringDampHat(idxPair, 2);

        for idxMass = 1:length(vecMassHat)
            scalMassHat = vecMassHat(idxMass);

            for idxV = 1:length(vecVImpactHat)
                scalVImpactHat = vecVImpactHat(idxV);

                % choose line style based on single switch
                if scalDottedBy == "mass"
                    scalLineStyle = vecLineStyleMass(idxMass);
                else
                    scalLineStyle = vecLineStyleV(idxV);
                end

                ratios           = nan(size(NArr));
                lambdas_theory   = nan(size(NArr));
                lambdas_measured = nan(size(NArr));

                for idxN = 1:length(NArr)
                    [ratios(idxN), lambdas_theory(idxN), lambdas_measured(idxN)] = ...
                        sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
                              scalVImpactHat=scalVImpactHat, scalGravityHat=scalGravityHat);
                end

                e = sqrt(max(ratios, 0));
                e(ratios <= .2) = NaN;

                semilogx(NArr, e, ...
                    'LineWidth', 1.5, ...
                    'Color', colors(idxPair,:), ...
                    'Marker', 'o', ...
                    'MarkerSize', vecMarkerSize(idxPair), ...
                    'LineStyle', scalLineStyle{1}, ...
                    'DisplayName', sprintf('$\\hat{k}=%.2f,\\ \\hat{\\gamma}=%.3f$', ...
                                           scalSpringHat, scalDampHat), ...
                    'Tag', 'dataLine');   % <--- added tag
            end
        end
    end

    xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);

    legend('show', ...
           'Location', 'best', ...
           'Interpreter', 'latex', ...
           'FontSize', 13, ...
           'NumColumns', 1);

    % xlim([6, 25]);
    % ylim([.76, .92]);
    grid on;
    xscale(gca, 'log');
    box on;
    % --- inset zoom: N in [8,16], e in [0.8,0.92] ---
    mainAx  = gca;

    % position = [left bottom width height] in figure-normalized units
    insetAx = axes('Position', [0.58 0.61 0.28 0.28]);
    box(insetAx, 'on'); hold(insetAx, 'on');

    % copy all line objects from the main axes into the inset
    % copy only data lines (tagged 'dataLine') from the main axes into the inset
    hLines = findobj(mainAx, 'Type', 'line', '-and', 'Tag', 'dataLine');
    copyobj(hLines, insetAx);

    % keep semilog x-scale and zoom into desired region
    set(insetAx, 'XScale', 'log');
    xlim(insetAx, [6 25]);
    ylim(insetAx, [0.76 0.92]);

    % optional: simplify inset (no labels, smaller font)
    set(insetAx, 'XTickLabel', [], 'YTickLabel', []);
    set(insetAx, 'FontSize', 10);
    grid on;
    box on;
    % -----------------------------------------------
    theme(gcf, 'light');
end
