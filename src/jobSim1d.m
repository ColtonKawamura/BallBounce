NArr           = 3:30;
vecDampHat     = [0.016, .015];
vecMassHat     = [1];               % ball-to-chain mass ratios
vecVImpactHat  = [0.2];           % impact velocities
vecSpringHat   = [1];

scalGravityHat = 0;

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

