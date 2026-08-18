function [scalRatioKE, scalLambda] = sim1d(scalHeightDrop, scalDampHat, scalNumPart, options)

    arguments
        scalHeightDrop (1,1) double
        scalDampHat (1,1) double
        scalNumPart (1,1) double
        options.scalDampHatBall    (1,1) double = 0        % defaults to chain value
        options.scalMassBall    (1,1) double = 1
        options.scalSpringBall  (1,1) double = 1 % ball spring constant
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
    scalMass = 1;
    scalSpringConst = 1;
    scalNatFreq = sqrt(scalSpringConst/scalMass);
    scalDamp = 2 * scalDampHat * sqrt(scalSpringConst * scalMass);
    scalDampBall = 2 * options.scalDampHatBall * sqrt(options.scalSpringBall * options.scalMassBall);
    scalGravity = options.scalGravityScale*scalNatFreq^2*2*scalHeightDrop;
    scalTimeTotal = 2*(2*(sqrt(2*scalHeightDrop/scalGravity)+pi/scalNatFreq)); % time to drop + 1/2  period
    % scalFreqCollision = 1/scalTimeTotal;
    scalTimeStep = pi*sqrt(scalMass/scalSpringConst)*0.005; % this is 1/1000 of a period
    scalTimeStepHalf = .5 * scalTimeStep;
    scalTimeStepHalfSquared = .5 * scalTimeStep^2;
    scalLogInterval = scalTimeStep/10;

%% pre allocations
    vecMass = ones(scalNumPart,1)*scalMass;
    vecMass(1) = options.scalMassBall;
    vecDamp = ones(scalNumPart,1) * scalDamp;
    vecDamp(1) = scalDampBall;
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
    vecSpringConst(1)     = options.scalSpringBall;  % sour=1→dest=2
    vecSpringConst(scalNumPart) = options.scalSpringBall;  % reverse: sour=2→dest=1
    
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
        vecDampForce = -vecDamp .* vecVelocityX;
        vecForceX = vecForceX + vecDampForce;

        %% gravity pull all particles
        %% grvaity is much weaker that the compressoin force
        vecForceX = vecForceX + vecMass * scalGravity;

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
    if options.visSim
        figure;
        plot(vecTime, matVelX(1,:));
        xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
        ylabel('$v_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);

        figure;
        plot(vecTime, 0.5*options.scalMassBall*matVelX(1,:).^2);
        xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
        ylabel('$K_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);
    end

    %% analysis
    ballVel = matVelX(1,:);

    % find moment ball separates from chain: last time step where ball overlaps particle 2
    ballPos  = matPosX(1,:);
    chainPos = matPosX(2,:);
    inContact = (chainPos - ballPos) < scalDiam; % 1 when touching, 0 when not in contact

    % find first contact with the chain
    idxFirstContact = find(inContact, 1, 'first');

    % find first index (AFTER FIRST CONTACT) when ball is no longer in contact
    % need to add idxFirstContact so it's relative to the original index list
    idxFirstSeparation = idxFirstContact + find(~inContact(idxFirstContact:end), 1, 'first');

    scalV_before = max(abs(ballVel(1:idxFirstContact))); % finds max from start to first contact

    % find the first index of first separation after initial contact
    idxSepRelative = find(~inContact(idxFirstContact:end), 1, 'first');

    % combine two indicies to get the index of first separation relative to the original index list
    idxFirstSeparation = idxFirstContact + idxSepRelative;

    % flip with index of first sepration
    notContact = ~inContact(idxFirstContact:end);

    % avoid false positives by finding wher 3 consecutive time steps are false
    consecFalse = find(conv(double(notContact), ones(1,3), 'valid') == 3, 1, 'first');

    if isempty(consecFalse)
        warning('[sim1d] No clean separation found — using idxFirstSeparation.');
        idxTrueSeparation = idxFirstSeparation;
    else
        idxTrueSeparation = idxFirstContact + consecFalse - 1;
    end

    scalV_after = max(abs(ballVel(idxTrueSeparation:end)));

    if scalV_before < 0.2 || scalV_before > 0.25
        warning('[sim1d] V_before:  %.4f\n', scalV_before);
    else
       fprintf('[analysis] v_before: %.4f\n', scalV_before);
    end
    fprintf('[analysis] v_after:  %.4f\n', scalV_after);

    scalRatioKE  = (scalV_after / scalV_before).^2;  % e^2 for consistency

    % wave travel time:

    % scalTauContactTheory = pi / sqrt(options.scalSpringBall / scalMass); % this is what the paper uses
    scalTauContactTheory = pi / sqrt(options.scalSpringBall / options.scalMassBall);
    scalTauContact = vecTime(idxFirstSeparation) - vecTime(idxFirstContact);
    scalLambdaTheory = 2*(scalNumPart-1)*scalDiam / (scalNatFreq*scalDiam*scalTauContactTheory);

    % to get scalLambda computationally, need to measure wavespeed directly
    % --- LEG 1: time to hit back wall ---
    % find first time gap between particle (N-1) and N < scalDiam
    vecGapLast = matPosX(scalNumPart,:) - matPosX(scalNumPart-1,:); % gap between last two particles

    % find() gives subarray, so need to add idxFirstContact to make relative to original
    % REPLACE the idxWaveArrival block with:
    velSecondToLast = matVelX(scalNumPart-1, :);
    velBackground = mean(abs(velSecondToLast(1:idxFirstContact-1)));
    idxWaveArrival = idxFirstContact + ...
    find(abs(velSecondToLast(idxFirstContact+1:end)) > velBackground + 1e-4, 1, 'first');


    % time it takes for  wave to first hit the "bottom"
    scalLeg1 = vecTime(idxWaveArrival) - vecTime(idxFirstContact);

    % --- LEG 2: reflection: time for particle (N-1) to return to zero velocity ---
    % find first time vel(N-1) crosses zero AFTER wave arrival
    velSecondToLast = matVelX(scalNumPart-1, :);
    idxReturnZero = idxWaveArrival + find(velSecondToLast(idxWaveArrival:end) >= 0, 1, 'first') - 1;
    scalLeg2 = vecTime(idxReturnZero) - vecTime(idxWaveArrival);

    scalRoundTrip= 2 * (scalLeg1 +scalLeg2);  % assuming symmetric return
    scalWaveSpeed = (scalNumPart-1)*scalDiam / (scalLeg1);

    % scalLambda = scalLambdaTheory;
    scalLambda = scalRoundTrip / scalTauContact;



%% Debugging
    fprintf('[analysis] e:        %.4f\n', sqrt(scalV_after/scalV_before));
    % fprintf('[analysis] idxFirstContact:    %d  (t=%.3f)\n',...
    %     idxFirstContact,...
    %     vecTime(idxFirstContact));
    % fprintf('[analysis] idxFirstSeparation: %d  (t=%.3f)\n',...
    %     idxFirstSeparation,...
    %     vecTime(idxFirstSeparation));
    fprintf('[analysis] tau_contact: %.4f\n', scalTauContact);
    fprintf('[analysis] tau_expected: %.4f\n', pi/sqrt(options.scalSpringBall/options.scalMassBall));
    fprintf('[analysis] c_theory: %.4f\n', scalNatFreq * scalDiam);
    fprintf('[analysis] c_measured: %.4f\n', scalWaveSpeed);
    % fprintf('[analysis] Lambda_corrected: %.4f\n',...
        % 2*(scalNumPart-1)*scalDiam / (5.8233 * scalTauContact));
    fprintf('[analysis] Lambda_theory: %.4f\n', scalLambdaTheory);
    fprintf('[analysis] Lambda_measured: %.4f\n', scalLambda);


    function vecForce = forceLaw(vecPosX, vecSour, vecDest, vecDistRest, vecSpringConst)
        vecContDist = vecPosX(vecDest) - vecPosX(vecSour);
        vecDistCurr = abs(vecContDist);
        vecOverlapX = vecDistRest - vecDistCurr;
        vecForceMag = -vecSpringConst .* (vecDistRest./vecDistCurr - 1) .* (vecOverlapX > 0);
        vecForce    = vecForceMag .* vecContDist;
    end
end














