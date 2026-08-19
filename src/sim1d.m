function [scalRatioKE, scalLambdaTheory, scalLambda_Measured] = sim1d(scalHeightDrop, scalDampHat, scalNumPart, scalMassRatio, scalSpringRatio,scalSpringConst, options)

    arguments
        scalHeightDrop (1,1) double
        scalDampHat (1,1) double
        scalNumPart (1,1) double
        scalMassRatio    (1,1) double 
        scalSpringRatio    (1,1) double
        scalSpringConst   (1,1) double  = 5
        options.scalDampHatBall    (1,1) double = scalDampHat       % defaults to chain valu
        options.visSim (1,1) logical = false
        options.scalPressure (1,1) double = 0.0
        options.scalGravityScale (1,1) double = 0.0001 % should be smaller than compression from impact
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

%%% temp initial parameters
    if options.visSim
        figVis = figure('Name', '1D Simulation');
        gifPath = fullfile('~/Desktop/', 'sim1d.gif');
    end

%% physical constant calculations
    scalDiam = 1;
    scalMass = 2;
    % scalMassBall = scalMassRatio * scalMass;
    scalMassBall = scalMassRatio * scalMass;
    % scalSpringBall = scalSpringRatio * scalSpringConst;
    scalSpringBall = 1;
    scalNatFreq = sqrt(scalSpringConst/scalMass);
    scalDamp = 2 * scalDampHat * sqrt(scalSpringConst * scalMass);
    scalDampBall = 2 * options.scalDampHatBall * sqrt(scalSpringBall * scalMassBall);
    scalGravity = options.scalGravityScale*scalNatFreq^2*2*scalHeightDrop;
    scalTimeTotal = (2*(sqrt(2*scalHeightDrop/scalGravity)+pi/scalNatFreq)); % time to drop + 1/2  period
    % scalFreqCollision = 1/scalTimeTotal;
    scalOmegaMax = sqrt(max(scalSpringConst, scalSpringBall) / min(scalMass, scalMassBall));
    scalTimeStep = pi / scalOmegaMax * 0.005;
    % scalTimeStep = pi*sqrt(scalMass/scalSpringConst)*0.005; % this is 1/1000 of a period
    scalTimeStepHalf = .5 * scalTimeStep;
    scalTimeStepHalfSquared = .5 * scalTimeStep^2;
    scalLogInterval = scalTimeStep/10;

%% pre allocations
    vecMass = ones(scalNumPart,1)*scalMass;
    vecMass(1) = scalMassBall;
    vecDamp = ones(scalNumPart,1) * scalDamp;
    vecDamp(1) = scalDampBall;
    % vecDamp(1) = 0;
    vecPosX = ((0:scalNumPart-1) * scalDiam * (1 - options.scalPressure))';
    vecPosX(1) = vecPosX(2) - scalDiam - scalHeightDrop;
    vecVelocityX = zeros(scalNumPart,1);
    % vecAccelerationX = zeros(scalNumPart,1);
    vecAccelerationX_old = zeros(scalNumPart,1);
    % vecForceX = zeros(scalNumPart,1);
    vecTime = (0:scalTimeStep:scalTimeTotal);
    scalMaxTimeSteps = length( vecTime );
    scalVisInterval = 100; % 

    % depabtable if you want these - more memory for bigger sims
    if savePosVel
        matPosX  = zeros(scalNumPart, scalMaxTimeSteps);
        matVelX  = zeros(scalNumPart, scalMaxTimeSteps);
        matAccelX  = zeros(scalNumPart, scalMaxTimeSteps);
    end

    % hertzian would be  t_contact propto (m^2/E^2 R v_0)^{1/5} or something

%% edge list (using graph theorysyntax)
    % produces a list of inditices of particles that are connected by a spring
    vecSour = [(1:scalNumPart-1)'; (2:scalNumPart)']; % index of particle on one end of spring
    vecDest = [(2:scalNumPart)'; (1:scalNumPart-1)']; % index of particle on other end of spring
    vecDistRest = repmat(scalDiam * (1 - options.scalPressure), 2*(scalNumPart-1), 1);
    vecSpringConst        = ones(2*(scalNumPart-1), 1) * scalSpringConst;
    vecSpringConst(1)     = scalSpringBall;  % sour=1→dest=2
    vecSpringConst(scalNumPart) = scalSpringBall;  % reverse: sour=2→dest=1
    
%% GPU transfer
    if useGPU
        vecPosX = gpuArray(vecPosX);
        vecVelocityX = gpuArray(vecVelocityX);
        % vecAccelerationX = gpuArray(vecAccelerationX);
        vecAccelerationX_old = gpuArray(vecAccelerationX_old);
    end

    % step 1 already set up befor with initial conditions
    ticLoop = tic;

    %
    if savePosVel
        matPosX(:, 1) = gather(vecPosX);
        matVelX(:, 1) = gather(vecVelocityX);
        matAccelX(:, 1) = gather(vecAccelerationX_old);
    end
    for step = 2:length(vecTime)

        if mod(step, scalLogInterval) == 0
            fprintf('[chain1d]   step %d / %d  (%.0f%%, %.1f s)\n', ...
                step, scalMaxTimeSteps, 100*step/scalMaxTimeSteps, toc(ticLoop));
        end

        %% update positions and velocities
        vecPosX = vecPosX + vecVelocityX*scalTimeStep + vecAccelerationX_old.*scalTimeStepHalfSquared;
        vecVelocityX = vecVelocityX + vecAccelerationX_old* scalTimeStepHalf;

        %% update forces and accelerations
        vecForceEdge = forceLaw(vecPosX, vecSour, vecDest, vecDistRest, vecSpringConst);
        % sum values in vecForceEdge by the groups specified in vecSour
        vecForceX = accumarray(vecSour, vecForceEdge, [scalNumPart, 1]); 

        %% damping force
        % vecDampForce = -vecDamp .* vecVelocityX;
        % vecForceX = vecForceX + vecDampForce;
        vecDampForce = -vecDamp .* vecVelocityX;
        inContactNow = (vecPosX(2) - vecPosX(1)) < scalDiam;
        vecDampForce(1) = -scalDampBall * vecVelocityX(1) * inContactNow;
        vecForceX = vecForceX + vecDampForce;

        %% gravity pull all particles
        %% grvaity is much weaker that the compressoin force
        vecForceX(1) = vecForceX(1) + vecMass(1) * scalGravity;

        %% fix endpoints
        vecForceX(scalNumPart) = 0;% the bottom (right) particle is fixed

        %% update acceleration and complete velocity verlet
        vecAccelerationX = vecForceX ./ vecMass;
        vecVelocityX = vecVelocityX + vecAccelerationX * scalTimeStepHalf;
        vecAccelerationX_old = vecAccelerationX;

        %% again depabtable, 
        if savePosVel
            matPosX(:, step) = gather(vecPosX); % gather() pulls from GPU to CPU
            matVelX(:, step) = gather(vecVelocityX);
            matAccelX(:, step) = gather(vecAccelerationX);
        end

        if options.visSim && mod(step, scalVisInterval) == 0
            figure(figVis); cla;
            for p = 1:scalNumPart
                theta = linspace(0, 2*pi, 64);
                cx = gather(vecPosX(p));
                plot(cx + (scalDiam/2)*cos(theta), (scalDiam/2)*sin(theta), 'w');
                hold on;
            end
            axis equal;
            drawnow;

            % capture frame
            frame = getframe(figVis);
            img   = frame2im(frame);
            [imgInd, cmap] = rgb2ind(img, 256);
            if step == scalVisInterval  % first frame
                imwrite(imgInd, cmap, gifPath, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
            else
                imwrite(imgInd, cmap, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
            end
        end
    end

%% plot

%% analysis


%% Find Peaks of KE for restitution
    vecBallVel = matVelX(1,:);
    vecKE = 0.5 * scalMassBall * vecBallVel.^2;

    % if the second peak is not found, set it to 0
    % use first two peaks in chronological order  in case they are about equal
    try 
        [vecPeakVals, vecPeakIdx] = findpeaks(vecKE);
        scalKE_before = vecPeakVals(1);
        scalKE_after  = vecPeakVals(2);
    catch ME
        scalKE_after  = 0;
    end
    fprintf('[analysis] KE_before: %.4f\n', scalKE_before);
    fprintf('[analysis] KE_after: %.4f\n', scalKE_after);

    scalRatioKE = scalKE_after / scalKE_before;
    scalV_before = sqrt(2*scalKE_before / scalMassBall);
    scalV_after  = sqrt(2*scalKE_after  / scalMassBall);
    fprintf('[analysis] e: %.4f\n', scalV_after/scalV_before);


%%  Contact Time

    % find moment bal gets in contact with chain
    vecBallPosX  = matPosX(1,:);
    vecChainLeadPosX = matPosX(2,:);
    inContact = (vecChainLeadPosX - vecBallPosX) < scalDiam; % 1 when touching, 0 when not in contact
    idxFirstContact = find(inContact, 1, 'first');

    % Find when velocity goes to zero after first KE peak
    vecBallVelWindow = vecBallVel(vecPeakIdx(1):end);
    scalIdxCrossZero = find(diff(sign(vecBallVelWindow)) ~= 0, 1, 'first');
    idxVelZero = vecPeakIdx(1) + scalIdxCrossZero + 1;

    % contact time = first contact to velocity zero crossing
    scalTauContact =.5.*( vecTime(idxVelZero) - vecTime(idxFirstContact));
    fprintf('[analysis] tau_contact: %.4f\n', scalTauContact);

    % "τ is half the total contact time for a symmetric impact"
    scalTauContactTheory = .5 * pi / sqrt(scalSpringBall / scalMassBall);
    fprintf('[analysis] tau_thoery: %.4f\n', pi/sqrt(scalSpringBall/scalMassBall));

%% measured wavespeed

    % velocity of second-to-last particle (first to feel the reflected wall)
    velSecondToLast = matVelX(scalNumPart-1, :);

    % baseline velocity before any wave arrives (should be ~0, nonzero if gravity on chain)
    velBackground = mean(abs(velSecondToLast(1:idxFirstContact-1)));

    % wave arrival: first time vel(N-1) exceeds background AFTER first contact
    % find() returns index relative to the subarray starting at idxFirstContact+1
    % so add idxFirstContact to convert back to global time index
    idxWaveArrival = idxFirstContact + ...
        find(abs(velSecondToLast(idxFirstContact+1:end)) > velBackground + 1e-4, 1, 'first');

    % LEG 1: time for wave to travel from particle 1 to particle N-1
    % distance = (N-1) * d, time = idxWaveArrival - idxFirstContact
    scalLeg1 = vecTime(idxWaveArrival) - vecTime(idxFirstContact);

    % wave speed = distance / time
    scalWaveSpeed_Measured = (scalNumPart-1)*scalDiam / scalLeg1;
    fprintf('[analysis] c_measured: %.4f\n', scalWaveSpeed_Measured);

    scalWaveSpeed_Theory = sqrt(scalSpringConst / scalMass) * scalDiam;
    fprintf('[analysis] c_theory: %.4f\n', scalWaveSpeed_Theory);

%% Lambda From Theory

    scalLambdaTheory = 2*(scalNumPart-1)*scalDiam / (scalNatFreq*scalDiam*scalTauContactTheory);
    fprintf('[analysis] Lambda_theory: %.4f\n', scalLambdaTheory);
    % scalLambda = scalRoundTrip / scalTauContact;

%% Lambda from Calculation
    scalLambda_Measured = 2*(scalNumPart-1)*scalDiam / (scalWaveSpeed_Measured.*scalTauContact);
    fprintf('[analysis] Lambda_measured: %.4f\n', scalLambda_Measured);

    % scalLambda = scalLambdaTheory;


%% plotting

    if options.visSim
        % figure;
        % plot(vecTime, matVelX(1,:));
        % xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
        % ylabel('$v_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);
        % grid on;


        figure;
        plot(vecTime, vecKE); hold on;
        plot(vecTime(vecPeakIdx(1)), vecPeakVals(1), 'rv', 'MarkerFaceColor','r', 'MarkerSize', 10);
        if length(vecPeakIdx) > 1
            plot(vecTime(vecPeakIdx(2)), vecPeakVals(2), 'gv', 'MarkerFaceColor','g', 'MarkerSize', 10);
        end
        legend('KE', '$K_{\mathrm{before}}$', '$K_{\mathrm{after}}$', 'Interpreter','LaTeX');
        xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
        ylabel('$K_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);
        xline(vecTime(idxFirstContact), '--w', 'contact',   'LabelVerticalAlignment','bottom');
    xline(vecTime(idxVelZero),      '--y', 'vel = 0',   'LabelVerticalAlignment','bottom');
        grid on;
    end

    function vecForce = forceLaw(vecPosX, vecSour, vecDest, vecDistRest, vecSpringConst)
        vecContDist = vecPosX(vecDest) - vecPosX(vecSour);
        vecDistCurr = abs(vecContDist);
        vecOverlapX = vecDistRest - vecDistCurr;
        vecForceMag = -vecSpringConst .* (vecDistRest./vecDistCurr - 1) .* (vecOverlapX > 0);
        vecForce    = vecForceMag .* vecContDist;
    end
end














