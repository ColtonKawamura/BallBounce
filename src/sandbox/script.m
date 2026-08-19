% All defaults addpath("~/repos/BallBounce/src/");
NArr = 2:15;
scalDampHat    = 0.001;
% scalH = 0.22 / (2 * sqrt(.0001));  % = 11
scalH = 5;  % = 11
scalMassRatio= 6; % ratio of ball to chain
scalSpringRatio = .0718; % ratio of ball to chain (increases contact time)
scalSpringRatio = .1; % ratio of ball to chain (increases contact time)
% scalSpringRatio = scalSpringRatio * (30/1000); % 30 mbar / 1000 mbar

ratios  = nan(size(NArr));
lambdas_theory = nan(size(NArr));
lambdas_measured = nan(size(NArr));

for i = 1:length(NArr)
    fprintf('N=%d\n', NArr(i));
    [ratios(i), lambdas_theory(i), lambdas_measured(i)] = sim1d(scalH,...
        scalDampHat,...
        NArr(i), ...
        scalMassRatio, ...
        scalSpringRatio ...
        );
end

e = sqrt(max(ratios, 0)); % i'm  outputing KE
% e_inf = mean(e(end-5:end), 'omitnan'); % average the last 5 values
e_inf = 0; % average the last 5 values
e_max = max(e, [], 'omitnan');
e_tilde = (e - e_inf) ./ (e_max - e_inf);

figure;
semilogx(lambdas_theory, e_tilde, 'o-', 'LineWidth', 2, 'DisplayName', 'Theory');
hold on;
semilogx(lambdas_measured, e_tilde, 'o-', 'LineWidth', 2, 'DisplayName', 'Measured');
xlabel('$\Lambda = 2Nd/c\tau$', Interpreter='latex', FontSize=20);
ylabel('$\tilde{e}$', Interpreter='latex', FontSize=20);
legend('show', 'Location', 'NorthEastOutside', 'Interpreter', 'latex', FontSize=15);
title(sprintf('$\\hat{\\gamma} = %.3f$', scalDampHat), 'Interpreter', 'latex', 'FontSize', 20);
grid on;

figure;
semilogx(NArr, e_tilde, 'o', LineWidth=2);
yline(1, '--r');
xline(1, '--k');
xlabel('$N$', Interpreter='latex', FontSize=20);
ylabel('$\tilde{e}$', Interpreter='latex', FontSize=20);
grid on;

sim1d(scalH,scalDampHat, 5, scalMassRatio, scalSpringRatio , visSim=true);


% ----
NArr = 2:15;
scalDampHat   = 0.0001;
scalH         = 5;
scalMassRatio = 6;
scalSpringRatioArr = [0.2, 0.15, 0.11, 0.1, 0.09, 0.08, 0.072];  % loop over these

figure; hold on;
colors = lines(length(scalSpringRatioArr));

for j = 1:length(scalSpringRatioArr)
    scalSpringRatio = scalSpringRatioArr(j);

    ratios         = nan(size(NArr));
    lambdas_theory = nan(size(NArr));
    lambdas_measured = nan(size(NArr));

    for i = 1:length(NArr)
        [ratios(i), lambdas_theory(i), lambdas_measured(i)] = sim1d(scalH, ...
            scalDampHat, NArr(i), scalMassRatio, scalSpringRatio);
    end

    e       = sqrt(max(ratios, 0));
    e(ratios <= .2) = NaN;
    e_inf   = 0;
    e_max   = max(e, [], 'omitnan');
    e_tilde = (e - e_inf) ./ (e_max - e_inf);

    t = linspace(1, 0, length(scalSpringRatioArr))';  % 1=red (high), 0=blue (low)
    colors = [t, zeros(length(t),1), 1-t];             % [R, G, B]

    semilogx(lambdas_theory, e_tilde, 'o-', ...
        'LineWidth', 2, ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('$k_r = %.3f$', scalSpringRatio));
end

xlabel('$\Lambda = 2Nd/c\tau$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\tilde{e}$',           'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
title(sprintf('$\\hat{\\gamma} = %.4f$', scalDampHat), 'Interpreter', 'latex', 'FontSize', 20);
grid on;
xscale(gca, 'log');



% ---
NArr = 2:15;
scalH         = 5;
scalMassRatio = 6;
scalSpringRatio = 0.1;  % fixed
scalDampHatArr = [0.01, 0.003, 0.001, 0.0003, 0.0001];  % high=red, low=blue

figure; hold on;
t = linspace(1, 0, length(scalDampHatArr))';
colors = [t, zeros(length(t),1), 1-t];

for j = 1:length(scalDampHatArr)
    scalDampHat = scalDampHatArr(j);

    ratios           = nan(size(NArr));
    lambdas_theory   = nan(size(NArr));
    lambdas_measured = nan(size(NArr));

    for i = 1:length(NArr)
        [ratios(i), lambdas_theory(i), lambdas_measured(i)] = sim1d(scalH, ...
            scalDampHat, NArr(i), scalMassRatio, scalSpringRatio);
    end

    e = sqrt(max(ratios, 0));
    e(ratios <= .2) = NaN;
    e_max   = max(e, [], 'omitnan');
    e_tilde = (e - 0) ./ (e_max - 0);

    markerSizes = linspace(15, 4, length(scalDampHatArr));  % large=high damp, small=low damp

    % then inside the loop, replace the semilogx call with:
    semilogx(lambdas_theory, e_tilde, '-', ...
        'LineWidth', 2, ...
        'Color', colors(j,:), ...
        'HandleVisibility', 'off');
    semilogx(lambdas_theory, e_tilde, 'o', ...
        'MarkerSize', markerSizes(j), ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('$\\hat{\\gamma} = %.4f$', scalDampHat));
end

xlabel('$\Lambda = 2Nd/c\tau$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\tilde{e}$',           'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
title(sprintf('$k_r = %.3f$', scalSpringRatio), 'Interpreter', 'latex', 'FontSize', 20);
grid on;
xscale(gca, 'log');
