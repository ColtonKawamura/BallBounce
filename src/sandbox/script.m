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

sim1d(scalH,scalDampHat, 5, scalMassRatio, scalSpringRatio ,10, visSim=true);


%% ---- Sweep over spring ratios
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



%% --- Swee pver damping
NArr = 2:15;
scalH         = 5;
scalMassRatio = 6;
scalSpringRatio = 0.1;  % fixed
scalDampHatArr = [0.01, 0.003, 0.001, 0.0003, 0.0001];  % high=red, low=blue

figure; hold on;
t_color = scalSpringRatio / 0.15;  % 1=red, 0=blue
lineColor = [t_color, 0, 1-t_color];
markerSizes = linspace(15, 4, length(scalDampHatArr));  % large=high damp, small=low damp

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

    semilogx(lambdas_theory, e_tilde, '-', ...
        'LineWidth', 1, ...
        'Color', lineColor, ...
        'HandleVisibility', 'off');
    semilogx(lambdas_theory, e_tilde, 'o', ...
        'MarkerSize', markerSizes(j), ...
        'Color', lineColor, ...
        'DisplayName', sprintf('$\\hat{\\gamma} = %.4f$', scalDampHat));
end

xlabel('$\Lambda = 2Nd/c\tau$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\tilde{e}$',           'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
title(sprintf('$k_r = %.3f$', scalSpringRatio), 'Interpreter', 'latex', 'FontSize', 20);
grid on;
xscale(gca, 'log');


%% --- Sweep over Chain Stiffness vs Lambda ---
NArr = 2:40;
scalDampHat   = 0.0001;
scalH         = 5;
scalMassRatio = 6;
scalSpringRatio = 0.2;  % fixed
scalSpringConstArr = [10, 7, 5, 3, 1, 0.5, 0.1];  % high=red (high pressure), low=blue (low pressure)

figure; hold on;
t = linspace(1, 0, length(scalSpringConstArr))';
colors = [t, zeros(length(t),1), 1-t];

for j = 1:length(scalSpringConstArr)
    scalSpringConst = scalSpringConstArr(j);

    ratios           = nan(size(NArr));
    lambdas_theory   = nan(size(NArr));
    lambdas_measured = nan(size(NArr));

    for i = 1:length(NArr)
        [ratios(i), lambdas_theory(i), lambdas_measured(i)] = sim1d(scalH, ...
            scalDampHat, NArr(i), scalMassRatio, scalSpringRatio, scalSpringConst);
    end

    e = sqrt(max(ratios, 0));
    e(ratios <= .2) = NaN;
    e_max   = max(e, [], 'omitnan');
    e_tilde = e ./ e_max;

    semilogx(lambdas_theory, e_tilde, 'o-', ...
        'LineWidth', 2, ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('$k = %.2f$', scalSpringHat));
end

xlabel('$\Lambda = 2Nd/c\tau$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\tilde{e}$',           'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
title(sprintf('Pressure sweep, $\\hat{\\gamma} = %.4f$, $k_r = %.2f$', scalDampHat, scalSpringRatio), ...
    'Interpreter', 'latex', 'FontSize', 20);
grid on;
xscale(gca, 'log');

%% --- Sweep over Chain Stiffness vs N ---
NArr = 2:40;
scalDampHat   = 0.0004;
scalH         = 5;
scalMassRatio = 6;
scalSpringRatio = 0.2;  % fixed
% scalSpringConstArr = [10, 7, 5, 3, 1, 0.5, 0.1];  % high=red (high pressure), low=blue (low pressure)
scalSpringConstArr = logspace(.6, 1.14, 10);  % high=red (high pressure), low=blue (low pressure)
% scalSpringConstArr = logspace(.6,1.5 , 10);  % high=red (high pressure), low=blue (low pressure)

figure; hold on;
t = linspace(1, 0, length(scalSpringConstArr))';
colors = [t, zeros(length(t),1), 1-t];
colors = flip(colors);

for j = 1:length(scalSpringConstArr)
    scalSpringConst = scalSpringConstArr(j);

    ratios           = nan(size(NArr));
    lambdas_theory   = nan(size(NArr));
    lambdas_measured = nan(size(NArr));

    for i = 1:length(NArr)
        [ratios(i), lambdas_theory(i), lambdas_measured(i)] = sim1d(scalH, ...
            scalDampHat, NArr(i), scalMassRatio, scalSpringRatio, scalSpringConst);
    end

    e = sqrt(max(ratios, 0));
    e(ratios <= .2) = NaN;
    e_max   = max(e, [], 'omitnan');
    e_tilde = e ./ e_max;

    % semilogx(lambdas_measured, e_tilde, 'o-', ...
    semilogx(NArr, e_tilde, 'o-', ...
        'LineWidth', 2, ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('$k = %.2f$', scalSpringConst));
end

xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
% xlabel('$\Lambda = 2Nd/c\tau$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\tilde{e}$',           'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
title(sprintf('Pressure sweep, $\\hat{\\gamma} = %.4f$, $k_r = %.2f$', scalDampHat, scalSpringRatio), ...
    'Interpreter', 'latex', 'FontSize', 20);
grid on;
xscale(gca, 'log');
%% ----------Restitution vs SpringConst

NArr_fixed = 8;  % fixed N
scalDampHat   = 0.0001;
scalH         = 5;
scalMassRatio = 6;
scalSpringRatio = 0.2;  % fixed
scalSpringConstArr = logspace(-1, 1, 40);  % sweep from 0.1 to 10

restitution = nan(size(scalSpringConstArr));

for j = 1:length(scalSpringConstArr)
    scalSpringConst = scalSpringConstArr(j);
    [ratio, ~, ~] = sim1d(scalH, scalDampHat, NArr_fixed, scalMassRatio, scalSpringRatio, scalSpringConst);
    e = sqrt(max(ratio, 0));
    restitution(j) = e;
end

figure;
loglog(scalSpringConstArr, restitution, 'o-', 'LineWidth', 2, 'Color', [0.2 0.4 0.8]);
xlabel('$k$ (chain spring constant)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$e$', 'Interpreter', 'latex', 'FontSize', 20);
title(sprintf('$N=%d$, $\\hat{\\gamma}=%.4f$, $k_r=%.2f$', NArr_fixed, scalDampHat, scalSpringRatio), ...
    'Interpreter', 'latex', 'FontSize', 20);
grid on;



%% Single one
scalDampHat = 0.0001;
scalMassHat = 1; % ratio of ball to chain
scalSpringHat = 2; % ratio of ball to chain
sim1d(scalDampHat, 14, scalMassHat, scalSpringHat, visSim=true);



%% double sweep
NArr = 2:20;
scalDampHat = 0.0001;
scalMassHat = 1; % ball to chain
% scalSpringHattArr = logspace(.6, 1.14, 10);  % high=red (high pressure), low=blue (low pressure)
scalSpringHatArr = [.9, .99,1, 1.1, 2, 2.1];

figure; hold on;
t = linspace(1, 0, length(scalSpringHatArr))';
colors = [t, zeros(length(t),1), 1-t];
colors = flip(colors);

for j = 1:length(scalSpringHatArr)
    scalSpringHat = scalSpringHatArr(j);

    ratios           = nan(size(NArr));
    lambdas_theory   = nan(size(NArr));
    lambdas_measured = nan(size(NArr));

    for i = 1:length(NArr)
        [ratios(i), lambdas_theory(i), lambdas_measured(i)] = sim1d(scalDampHat, NArr(i), scalMassHat, scalSpringHat, scalVImpactHat=.5);
    end

    e = sqrt(max(ratios, 0));
    e(ratios <= .2) = NaN;
    e_max   = max(e, [], 'omitnan');
    e_tilde = e ./ e_max;

    % semilogx(lambdas_measured, e_tilde, 'o-', ...
    semilogx(NArr, e_tilde, 'o-', ...
        'LineWidth', 2, ...
        'Color', colors(j,:), ...
        'DisplayName', sprintf('$k = %.2f$', scalSpringHat));
end

xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 20);
% xlabel('$\Lambda = 2Nd/c\tau$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\tilde{e}$',           'Interpreter', 'latex', 'FontSize', 20);
legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
xscale(gca, 'log');
