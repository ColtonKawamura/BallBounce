NArr           = 3:60;
vecDampHat     = [0.0006, .0023];
vecMassHat     = [20];               % ball-to-chain mass ratios
vecVImpactHat  = [0.1];           % impact velocities
% vecSpringHat   = [1, 1.25, 1.5, 2.5, 3, 3.5, 3.7];
vecSpringHat   = [.5,1, 2.5];
scalGravityHat = 0;
% vecDampHat = .005./(.5*sqrt(vecSpringHat * 1))


% --- single switch: choose which parameter controls dotted vs solid ---
scalDottedBy = "mass";    % set to "mass" or "v"
% ----------------------------------------------------------------------

figure; hold on;

% colors per spring value
t      = linspace(1, 0, length(vecSpringHat))';
colors = flip([t, zeros(length(t),1), 1-t]);

% marker sizes per damping (larger damping -> larger markers)
scalDampMin   = min(vecDampHat);
scalDampMax   = max(vecDampHat);
vecMarkerSize = 4 + 8 * (vecDampHat - scalDampMin) / (scalDampMax - scalDampMin);  % 4–12

% line styles keyed by mass and by impact velocity
vecLineStyleMass = repmat({"-"}, size(vecMassHat));     % default solid
vecLineStyleMass(2:end) = {":"};                        % higher masses dotted

vecLineStyleV    = repmat({"-"}, size(vecVImpactHat));  % default solid
vecLineStyleV(2:end)    = {":"};                        % second/others dotted

for idxSpring = 1:length(vecSpringHat)
    scalSpringHat = vecSpringHat(idxSpring);

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

            for idxDamp = 1:length(vecDampHat)
                scalDampHat = vecDampHat(idxDamp);

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
                    'Color', colors(idxSpring,:), ...
                    'Marker', 'o', ...
                    'MarkerSize', vecMarkerSize(idxDamp), ...
                    'LineStyle', scalLineStyle{1}, ...
                    'DisplayName', sprintf('$\\hat{k}=%.2f,\\ \\hat{\\gamma}=%.3f,\\ \\hat{m}=%.2f,\\ \\hat{v}=%.2f$', ...
                                           scalSpringHat, scalDampHat, scalMassHat, scalVImpactHat));
            end
        end
    end
end

xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
xscale(gca, 'log');
title('$\hat{g} = 0$', 'Interpreter', 'latex', 'FontSize', 20);


%%

% --- second figure with only two selected curves ---
figure; hold on;

% indices for requested parameter combinations
idxSpringMax = find(vecSpringHat == max(vecSpringHat), 1);  % largest k
idxSpringMin = find(vecSpringHat == min(vecSpringHat), 1);  % smallest k
idxDampMin   = find(vecDampHat   == min(vecDampHat),   1);  % smallest damping
idxDampMax   = find(vecDampHat   == max(vecDampHat),   1);  % largest damping

scalMassHat    = vecMassHat(1);
scalVImpactHat = vecVImpactHat(1);

% --- curve 1: largest k, smallest damping ---
scalSpringHat = vecSpringHat(idxSpringMax);
scalDampHat   = vecDampHat(idxDampMin);

ratios = nan(size(NArr));
for idxN = 1:numel(NArr)
    [ratios(idxN), lambdas_theory, lambdas_measured] = ...
        sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
              scalVImpactHat=scalVImpactHat, scalGravityHat=scalGravityHat);
end

e1 = sqrt(max(ratios, 0));
e1(ratios <= .2) = NaN;

semilogx(NArr, e1, ...
    'LineWidth', 1.8, ...
    'Color', colors(idxSpringMax,:), ...
    'Marker', 'o', ...
    'MarkerSize', vecMarkerSize(idxDampMin), ...
    'LineStyle', '-', ...
    'DisplayName', sprintf('$\\hat{k}=%.2f,\\ \\hat{\\gamma}=%.4f,\\ \\hat{m}=%.2f,\\ \\hat{v}=%.2f$', ...
                           scalSpringHat, scalDampHat, scalMassHat, scalVImpactHat));

% --- curve 2: smallest k, largest damping ---
scalSpringHat = vecSpringHat(idxSpringMin);
scalDampHat   = vecDampHat(idxDampMax);

ratios = nan(size(NArr));
for idxN = 1:numel(NArr)
    [ratios(idxN), lambdas_theory, lambdas_measured] = ...
        sim1d(scalDampHat, NArr(idxN), scalMassHat, scalSpringHat, ...
              scalVImpactHat=scalVImpactHat, scalGravityHat=scalGravityHat);
end

e2 = sqrt(max(ratios, 0));
e2(ratios <= .2) = NaN;

semilogx(NArr, e2, ...
    'LineWidth', 1.8, ...
    'Color', colors(idxSpringMin,:), ...
    'Marker', 'o', ...
    'MarkerSize', vecMarkerSize(idxDampMax), ...
    'LineStyle', '--', ...
    'DisplayName', sprintf('$\\hat{k}=%.2f,\\ \\hat{\\gamma}=%.4f,\\ \\hat{m}=%.2f,\\ \\hat{v}=%.2f$', ...
                           scalSpringHat, scalDampHat, scalMassHat, scalVImpactHat));

xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
xscale(gca, 'log');
title('$\hat{g} = 0$ (two selected curves)', 'Interpreter', 'latex', 'FontSize', 20);
