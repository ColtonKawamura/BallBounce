function [scalRatioKE, scalLambdaTheory, scalLambda_Measured] = sim1d(scalDampHat, scalNumPart, scalMassHat, scalSpringHat, options)
    % sim1d: 1D ball–chain simulation in fully dimensionless units
%
% Inputs (all hats, dimensionless):
%   scalDampHat     = gamma_hat for the chain (normalized damping)
%   scalNumPart     = N (ball + chain particles)
%   scalMassHat     = m_b/m_c (ball-to-chain mass ratio)
%   scalSpringHat   = k_b/k_c (ball-to-chain spring ratio)
%
% Options:
%   options.scalDampHatBall   = ball damping (same normalization)   [default 1e-5]
%   options.visSim            = logical for visualization           [default false]
%   options.scalPressure      = precompression fraction             [default 0]
%   options.scalGravityHat    = g_hat = g m_c/(k_c d_c)             [default 1e-3]
%   options.scalVImpactHat    = v_hat = v0/(d_c sqrt(k_c/m_c))      [default 0.2]

    arguments
        scalDampHat          (1,1) double
        scalNumPart          (1,1) double
        scalMassHat          (1,1) double  % m_b/m_c
        scalSpringHat        (1,1) double  % k_b/k_c
        options.scalDampHatBall   (1,1) double = 1e-5
        options.visSim             (1,1) logical = false
        options.scalPressure       (1,1) double = 0.0
        options.scalGravityHat     (1,1) double = .001
        options.scalVImpactHat     (1,1) double = 0.8
        options.scalDiamHat        (1,1) double = 1.0
    end

    %% GPU setup
    useGPU = canUseGPU();
    if useGPU
        gdev = gpuDevice();
        fprintf('[sim] GPU detected: %s (%.1f GB)\n', gdev.Name, gdev.TotalMemory/1024^3);
    else
        fprintf('[sim] No GPU found — running on CPU.\n');
    end

    savePosVel = true;

    if options.visSim
        figVis = figure('Name', '1D Simulation');
        gifPath = fullfile('~/Desktop/', 'sim1d.gif');
    end

    %% physical constant calculations (dimensionless hats, chain-based)

    % base chain scales: d_c* = 1, m_c* = 1, k_c* = 1
    scalDiamChain   = 1;  % d_c
    scalMassChain   = 1;  % m_c
    scalSpringChain = 1;  % k_c

    % interpret input hats:
    %   scalMassHat   = m_b/m_c
    %   scalSpringHat = k_b/k_c
    %
    % so in these units:
    %   m_b* = scalMassHat,   m_c* = 1
    %   k_b* = scalSpringHat, k_c* = 1
    scalDiamBall   = scalDiamChain* options.scalDiamHat;
    scalMassBall   = scalMassHat   * scalMassChain;
    scalSpringBall = scalSpringHat * scalSpringChain;

    % natural frequencies
    scalNatFreqChain = sqrt(scalSpringChain / scalMassChain);  % chain
    scalNatFreqBall  = sqrt(scalSpringBall      / scalMassBall);    % ball

    % damping (chain and ball)
    scalDamp     = 2 * scalDampHat             * sqrt(scalSpringChain * scalMassChain);
    scalDampBall = 2 * options.scalDampHatBall * sqrt(scalSpringBall       * scalMassBall);

    % gravity: hat{g} = g m_c/(k_c d_c); here m_c = k_c = d_c = 1
    scalGravity  = options.scalGravityHat;

    % total simulation time ~ one round trip along the chain
    scalWavespeed = sqrt(scalSpringChain / scalMassChain) * scalDiamChain;
    scalTimeTotal = 40 * scalNumPart * scalDiamChain / scalWavespeed;

    % time step: fraction of ball period
    scalOmega     = scalNatFreqBall;
    scalPeriod    = 2 * pi / scalOmega;
    scalTimeStep  = scalPeriod * 0.005;

    scalTimeStepHalf        = 0.5 * scalTimeStep;
    scalTimeStepHalfSquared = 0.5 * scalTimeStep^2;
    scalLogInterval         = 1000;

    %% dimensionless impact velocity
    scalVImpact = options.scalVImpactHat;   % v_hat
    fprintf('[params] vImpact_hat = %.3f\n', scalVImpact);

    % displacement per timestep at impact (in diameters)
    scalStepDisp = scalVImpact * scalTimeStep / scalDiamChain;
    fprintf('[params] v_hat * dt / d_c = %.3f (diameters per step)\n', scalStepDisp);

    if scalStepDisp > 0.2
        warning(['[sim] Impact displacement per step is %.2f diameters (v_hat dt / d_c). ' ...
                 'Reduce scalTimeStep or scalVImpact.'], ...
                 scalStepDisp);
    end

    %% pre allocations

    % masses
    vecMass        = ones(scalNumPart,1) * scalMassChain;   % chain
    vecMass(1)     = scalMassBall;                     % ball

    % damping
    vecDamp        = ones(scalNumPart,1) * scalDamp;
    vecDamp(1)     = scalDampBall;

    % positions (precompressed by "pressure" factor)
    vecPosX = ((0:scalNumPart-1) * scalDiamChain * (1 - options.scalPressure))';

    % ball starts two diameter above first chain particle
    vecPosX(1) = vecPosX(2) - 2*scalDiamChain;

    % velocities and accelerations
    vecVelocityX        = zeros(scalNumPart,1);
    vecAccelerationX    = zeros(scalNumPart,1);
    vecAccelerationX_old= zeros(scalNumPart,1);

    vecTime          = 0:scalTimeStep:scalTimeTotal;
    scalMaxTimeSteps = length(vecTime);
    scalVisInterval  = 10;

    if savePosVel
        matPosX   = zeros(scalNumPart, scalMaxTimeSteps);
        matVelX   = zeros(scalNumPart, scalMaxTimeSteps);
        matAccelX = zeros(scalNumPart, scalMaxTimeSteps);
    end

    %% edge list (spring connectivity)

    vecSour = [(1:scalNumPart-1)'; (2:scalNumPart)'];
    vecDest = [(2:scalNumPart)';   (1:scalNumPart-1)'];

    vecDistRest    = repmat(scalDiamChain * (1 - options.scalPressure), 2*(scalNumPart-1), 1);
    vecSpringConst = ones(2*(scalNumPart-1), 1) * scalSpringChain;
    vecSpringConst(1)           = scalSpringBall;   % ball -> first chain
    vecSpringConst(scalNumPart) = scalSpringBall;   % first chain -> ball

    % zero velocities/accelerations, then set impact velocity of ball
    vecVelocityX(:)        = 0;
    vecAccelerationX(:)    = 0;
    vecAccelerationX_old(:)= 0;
    vecVelocityX(1)        = scalVImpact;

    %% GPU transfer

    if useGPU
        vecPosX             = gpuArray(vecPosX);
        vecVelocityX        = gpuArray(vecVelocityX);
        vecAccelerationX_old= gpuArray(vecAccelerationX_old);
    end

    ticLoop = tic;

    if savePosVel
        matPosX(:, 1)   = gather(vecPosX);
        matVelX(:, 1)   = gather(vecVelocityX);
        matAccelX(:, 1) = gather(vecAccelerationX_old);
    end

    %% main time integration

    for step = 2:scalMaxTimeSteps

        if mod(step, scalLogInterval) == 0
            fprintf('[chain1d]   step %d / %d  (%.0f%%, %.1f s)\n', ...
                step, scalMaxTimeSteps, 100*step/scalMaxTimeSteps, toc(ticLoop));
        end

        % update positions and velocities
        vecPosX      = vecPosX + vecVelocityX*scalTimeStep + vecAccelerationX_old.*scalTimeStepHalfSquared;
        vecVelocityX = vecVelocityX + vecAccelerationX_old * scalTimeStepHalf;

        % spring forces
        vecForceEdge = forceLaw(vecPosX, vecSour, vecDest, vecDistRest, vecSpringConst);
        vecForceX    = accumarray(vecSour, vecForceEdge, [scalNumPart, 1]);

        % damping forces
        vecDampForce = -vecDamp .* vecVelocityX;
        inContactNow = (vecPosX(2) - vecPosX(1)) < scalDiamChain;
        vecDampForce(1) = -scalDampBall * vecVelocityX(1) * inContactNow;
        vecForceX = vecForceX + vecDampForce;

        % gravity on all particles
        % vecForceX = vecForceX + vecMass * scalGravity;
        % gravity on ball
        vecForceX(1) = vecForceX(1) + vecMass(1) * scalGravity;

        % bottom particle fixed
        vecForceX(scalNumPart) = 0;

        % update acceleration and complete velocity Verlet
        vecAccelerationX      = vecForceX ./ vecMass;
        vecVelocityX          = vecVelocityX + vecAccelerationX * scalTimeStepHalf;
        vecAccelerationX_old  = vecAccelerationX;

        if savePosVel
            matPosX(:, step)   = gather(vecPosX);
            matVelX(:, step)   = gather(vecVelocityX);
            matAccelX(:, step) = gather(vecAccelerationX);
        end

        if options.visSim && mod(step, scalVisInterval) == 0
            figure(figVis); cla;
            for p = 1:scalNumPart
                theta = linspace(0, 2*pi, 64);
                cx = gather(vecPosX(p));
                plot(cx + (scalDiamChain/2)*cos(theta), (scalDiamChain/2)*sin(theta), 'w');
                hold on;
            end
            axis equal;
            drawnow;

            frame = getframe(figVis);
            img   = frame2im(frame);
            [imgInd, cmap] = rgb2ind(img, 256);
            if step == scalVisInterval
                imwrite(imgInd, cmap, gifPath, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
            else
                imwrite(imgInd, cmap, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
            end
        end
    end

    %% analysis: restitution

    vecBallVel = matVelX(1,:);
    vecKE      = 0.5 * scalMassBall * vecBallVel.^2;

    try
        [vecPeakVals, vecPeakIdx] = findpeaks(vecKE);
        scalKE_before = vecPeakVals(1);
        scalKE_after  = vecPeakVals(2);
    catch
        scalKE_before = vecKE(1);
        scalKE_after  = 0;
    end
    fprintf('[analysis] KE_before: %.4f\n', scalKE_before);
    fprintf('[analysis] KE_after:  %.4f\n', scalKE_after);

    scalRatioKE  = scalKE_after / scalKE_before;
    scalV_before = sqrt(2*scalKE_before / scalMassBall);
    scalV_after  = sqrt(2*scalKE_after  / scalMassBall);
    fprintf('[analysis] e: %.4f\n', scalV_after/scalV_before);

    %% contact time

    vecBallPosX      = matPosX(1,:);
    vecChainLeadPosX = matPosX(2,:);
    inContact        = (vecChainLeadPosX - vecBallPosX) < scalDiamChain;
    idxFirstContact  = find(inContact, 1, 'first');

    vecBallVelWindow = vecBallVel(vecPeakIdx(1):end);
    scalIdxCrossZero = find(diff(sign(vecBallVelWindow)) ~= 0, 1, 'first');
    idxVelZero       = vecPeakIdx(1) + scalIdxCrossZero + 1;

    if isempty(idxVelZero) || isempty(idxFirstContact)
        scalTauContact = NaN;
    else
        scalTauContact = 0.5 * (vecTime(idxVelZero) - vecTime(idxFirstContact));
    end
    fprintf('[analysis] tau_contact: %.4f\n', scalTauContact);

    scalTauContactTheory = 0.5 * pi / sqrt(scalSpringBall / scalMassBall);
    fprintf('[analysis] tau_theory:  %.4f\n', scalTauContactTheory);

    %% measured wavespeed

    velSecondToLast = matVelX(scalNumPart-1, :);
    velBackground   = mean(abs(velSecondToLast(1:idxFirstContact-1)));

    idxWaveArrivalRel = find(abs(velSecondToLast(idxFirstContact+1:end)) ...
                             > velBackground + 1e-4, 1, 'first');
    if isempty(idxWaveArrivalRel)
        warning('idxWaveArrival did not detect wave hitting bottom');
        idxWaveArrival = length(vecTime);
    else
        idxWaveArrival = idxFirstContact + idxWaveArrivalRel;
    end

    scalLeg1 = vecTime(idxWaveArrival) - vecTime(idxFirstContact);
    scalWaveSpeed_Measured = (scalNumPart-1)*scalDiamChain / scalLeg1;
    fprintf('[analysis] c_measured: %.4f\n', scalWaveSpeed_Measured);

    scalWaveSpeed_Theory = sqrt(scalSpringChain / scalMassChain) * scalDiamChain;
    fprintf('[analysis] c_theory:   %.4f\n', scalWaveSpeed_Theory);

%% lambda
    scalLambdaTheory = 2*(scalNumPart-1)*scalDiamChain / (scalNatFreqChain*scalDiamChain*scalTauContactTheory);
    fprintf('[analysis] Lambda_theory:   %.4f\n', scalLambdaTheory);

    scalLambda_Measured = 2*(scalNumPart-1)*scalDiamChain / (scalWaveSpeed_Measured.*scalTauContact);
    fprintf('[analysis] Lambda_measured: %.4f\n', scalLambda_Measured);

    %% plotting

    if options.visSim
        figure;
        plot(vecTime, vecKE); hold on;
        plot(vecTime(vecPeakIdx(1)), vecPeakVals(1), 'rv', 'MarkerFaceColor','r', 'MarkerSize', 10);
        if length(vecPeakIdx) > 1
            plot(vecTime(vecPeakIdx(2)), vecPeakVals(2), 'gv', 'MarkerFaceColor','g', 'MarkerSize', 10);
        end
        legend('KE', '$K_{\mathrm{before}}$', '$K_{\mathrm{after}}$', 'Interpreter','LaTeX');
        xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
        ylabel('$K_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);
        xline(vecTime(idxFirstContact), '--w', 'contact', 'LabelVerticalAlignment','bottom');
        if idxVelZero ~= 0
            xline(vecTime(idxVelZero), '--y', 'vel = 0', 'LabelVerticalAlignment','bottom');
        end
        grid on;
    end

    %% nested force law (contact springs)
    function vecForce = forceLaw(vecPosX_loc, vecSour_loc, vecDest_loc, vecDistRest_loc, vecSpringConst_loc)
        vecContDist = vecPosX_loc(vecDest_loc) - vecPosX_loc(vecSour_loc);
        vecDistCurr = abs(vecContDist);
        vecOverlapX = vecDistRest_loc - vecDistCurr;
        vecForceMag = -vecSpringConst_loc .* (vecDistRest_loc./vecDistCurr - 1) .* (vecOverlapX > 0);
        vecForce    = vecForceMag .* vecContDist;
    end
end

