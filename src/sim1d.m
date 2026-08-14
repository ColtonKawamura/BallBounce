function sim1d(scalHeightDrop, scalDamp, scalPressure, scalSpringConst, scalNumPart, scalDiam, scalMass, scalGravity, options)

    arguments
        scalHeightDrop (1,1) double
        scalDamp (1,1) double
        scalPressure (1,1) double
        scalSpringConst (1,1) double
        scalNumPart (1,1) double
        scalDiam (1,1) double
        scalMass (1,1) double
        scalGravity (1,1) double
        options.visSim (1,1) logical = false
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
    scalLogInterval = 100;

%%% temp initial parameters
    % scalNumPart = 5;
    % scalDamp = .1;
    % scalDiam =1;
    % scalMass = 1;
    % scalSpringConst = 1;
    % scalHeightDrop = 10;
    % scalGravity = 3;
    % scalPressure = .00001;
    if options.visSim
        figVis = figure('Name', '1D Simulation');
        gifPath = fullfile('~/Desktop/', 'sim1d.gif');
    end

%% physical constant calculations
     scalNatFreq = sqrt(scalSpringConst/scalMass);
    scalTimeTotal = 2*(sqrt(2*scalHeightDrop/scalGravity)+pi/scalNatFreq); % time to drop + 1/2  period
    scalFreqCollision = 1/scalTimeTotal;
    scalTimeStep = pi*sqrt(scalMass/scalSpringConst)*0.05; % this is 1/1000 of a period
    scalTimeStepHalf = .5 * scalTimeStep;
    scalTimeStepHalfSquared = .5 * scalTimeStep^2;

%% pre allocations
    vecPosX = ((0:scalNumPart-1) * scalDiam * (1 - scalPressure))';
    vecPosX(1) = vecPosX(2) - scalDiam - scalHeightDrop;  % <-- ADD THIS LINE
    vecVelocityX = zeros(scalNumPart,1);
    vecAccelerationX = zeros(scalNumPart,1);
    vecAccelerationX_old = zeros(scalNumPart,1);
    vecForceX = zeros(scalNumPart,1);
    vecTime = (0:scalTimeStep:scalTimeTotal);
    scalMaxTimeSteps = length( vecTime );
    scalVisInterval = .005*scalSpringConst; % plot every 50 steps

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
    vecDistRest = repmat(scalDiam * (1 - scalPressure), 2*(scalNumPart-1), 1);
    
%% GPU transfer
    if useGPU
        vecPosX = gpuArray(vecPosX);
        vecVelocityX = gpuArray(vecVelocityX);
        vecAccelerationX = gpuArray(vecAccelerationX);
        vecAccelerationX_old = gpuArray(vecAccelerationX_old);
    end

    % step 1 already set up befor with initial conditions
    ticLoop = tic;

    %
    if savePosVel
        matPosX(:, 1) = gather(vecPosX);
        matVelX(:, 1) = gather(vecVelocityX);
        matAccelXX(:, 1) = gather(vecAccelerationX_old);
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
    vecContDist = vecPosX(vecDest) - vecPosX(vecSour);        % signed displacement src→dst
    vecDistCurr = abs(vecContDist);                            % scalar distance
    vecOverlapX = vecDistRest - vecDistCurr;                   % positive when compressed
    vecForceMag = -scalSpringConst .* (vecDistRest./vecDistCurr - 1) .* (vecOverlapX > 0);
    vecForceEdge = vecForceMag .* vecContDist;                 % force on source particle
    vecForceX = accumarray(vecSour, vecForceEdge, [scalNumPart, 1]);

    %% damping force
    vecDampForce = -scalDamp * vecVelocityX;
    vecForceX = vecForceX + vecDampForce;

    %% gravity pull particle 1
    vecForceX(1) = vecForceX(1) + scalMass * scalGravity;

    %% fix endpoints
    vecForceX(scalNumPart) = 0;% the bottom (right) particle is fixed

    %% update acceleration and complete velocity verlet
    vecAccelerationX = vecForceX / scalMass;
    vecVelocityX = vecVelocityX + vecAccelerationX * scalTimeStepHalf;
    vecAccelerationX_old = vecAccelerationX;

    %% again depabtable, 
    if savePosVel
        matPosX(:, step) = gather(vecPosX); % gather() pulls from GPU to CPU
        matVelX(:, step) = gather(vecVelocityX);
        matAccelX(:, step) = gather(vecVelocityX);
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
figure;
plot(vecTime, matVelX(1,:));
xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel('$v_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);

figure;
plot(vecTime, 0.5*scalMass*matVelX(1,:).^2);
xlabel('time', 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel('$K_{\mathrm{particle\ 1}}$', 'Interpreter', 'LaTeX', 'FontSize', 20);


    function force = forceLaw(vecOverlap, scalSpringConst)
        % Linear (Hookean) — swap this body for Hertzian etc.
            force = scalSpringConst * max(vecOverlap, 0);
    end
end














