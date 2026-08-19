% All defaults addpath("~/repos/BallBounce/src/");

scalDampHat    = 0;
% scalH          = 6;
scalH = 0.22 / (2 * sqrt(.0001));  % = 11
scalH = 5;  % = 11
scalMassBall= 12;
scalSpringBall = .4;                          % stiff — no tunneling

ratios  = nan(size(NArr));
lambdas = nan(size(NArr));

% debug
% [ratios(i), lambdas(i)] = sim1d(scalH,...
%         scalDampHat,...
%         6, ...
%         scalSpringBall=scalSpringBall,...
%         visSim=true);

for i = 1:length(NArr)
    fprintf('N=%d\n', NArr(i));
    [ratios(i), lambdas(i)] = sim1d(scalH,...
        scalDampHat,...
        NArr(i), ...
        scalSpringBall=scalSpringBall, ...
        scalMassBall=scalMassBall ...
        );
    fprintf('N=%d  Lambda=%.3f  e=%.4f\n', NArr(i), lambdas(i), sqrt(ratios(i)));
end

e = sqrt(max(ratios, 0)); % i'm  outputing KE
e_inf = mean(e(end-5:end), 'omitnan'); % average the last 5 values
e_max = max(e, [], 'omitnan');
e_tilde = (e - e_inf) ./ (e_max - e_inf);

figure;
semilogx(lambdas, e_tilde, 'o', LineWidth=2);
yline(1, '--r');
xline(1, '--k');
xlabel('$\Lambda = 2Nd/c\tau$', Interpreter='latex', FontSize=20);
ylabel('$\tilde{e}$', Interpreter='latex', FontSize=20);
grid on;

figure;
semilogx(NArr, e_tilde, 'o', LineWidth=2);
yline(1, '--r');
xline(1, '--k');
xlabel('$N$', Interpreter='latex', FontSize=20);
ylabel('$\tilde{e}$', Interpreter='latex', FontSize=20);
grid on;

sim1d(scalH,scalDampHat, 30, scalSpringBall=scalSpringBall, scalMassBall = scalMassBall, visSim=true);
