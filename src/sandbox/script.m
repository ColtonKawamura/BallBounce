% All defaults
addpath("~/repos/BallBounce/src/");


chain1d();                              % all defaults
chain1d(20, 1.0, .01, 2.0, 2.0, 500); % heavier mass, softer springs
chain1d(20, 1.0, 1, 1, 0.5, 2e3); % lighter mass, stiffer springs


%% GPU setup
    useGPU = canUseGPU();
    if useGPU
        gdev = gpuDevice();
        fprintf('[sim] GPU detected: %s (%.1f GB)\n', gdev.Name, gdev.TotalMemory/1024^3);
    else
        fprintf('[sim] No GPU found — running on CPU.\n');
    end



scalNumParticles = 5;
scalDiameter =1;
scalMass = 1;
scalSpringConstant = 1;
scalNatFreq = sqrt(scalSpringConstant/scalMass);
scalHeightDrop = 10;
scalGravity = 3;

scalTimeTotal = sqrt(2*scalHeightDrop/scalGravity)+pi/scalNatFreq;

vecPosX = (0:scalNumParticles-1)*scalDiameter;
vecVelocityX = zeros(scalNumParticles,1);
vecAccelerationX = zeros(scalNumParticles,1);
vecForceX = zeros(scalNumParticles,1);
vecAccelX = zeros(scalNumParticles,1);

scalTimeStep = pi*sqrt(scalMass/scalSpringConstant)*0.005;
vecTime = (0:scalTimeTotal)*scalTimeStep;

for step = 1:length(vecTime)

% hertzian would be  t_contact propto (m^2/E^2 R v_0)^{1/5} or something


