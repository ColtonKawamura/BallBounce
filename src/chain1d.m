function varargout = chain1d(N, diameter, f_drive, t_end, mass, k_spring)
%CHAIN1D  Simulate a 1D particle chain with fixed last particle + animation.
%
%   chain1d(N, diameter, f_drive, t_end, mass, k_spring)
%
%   Defaults:
%     N=20, diameter=1.0, f_drive=5.0, t_end=2.0, mass=1.0, k_spring=1e3

    if nargin < 1, N        = 20;   end
    if nargin < 2, diameter = 1.0;  end
    if nargin < 3, f_drive  = 5.0;  end
    if nargin < 4, t_end    = 2.0;  end
    if nargin < 5, mass     = 1.0;  end
    if nargin < 6, k_spring = 1e3;  end

    %% Derived / fixed parameters
    gamma   = 0.5;
    A_drive = 0.3 * diameter;
    omega   = 2*pi * f_drive;
    dt      = 1e-4;
    t       = 0:dt:t_end;
    Nt      = length(t);

    c     = diameter * sqrt(k_spring / mass);
    f_cut = sqrt(k_spring / mass) / pi;
    fprintf('Wave speed ~ %.2f m/s   |   Cutoff freq ~ %.2f Hz\n', c, f_cut);
    if f_drive > f_cut
        warning('f_drive (%.2f Hz) exceeds cutoff (%.2f Hz) — wave will not propagate!', ...
                f_drive, f_cut);
    end

    %% Initial conditions
    x0 = (0:N-1) * diameter;
    x  = x0;
    v  = zeros(1, N);
    a  = zeros(1, N);
    X  = zeros(N, Nt);
    V  = zeros(N, Nt);
    X(:,1) = x;

    %% Time integration (Velocity Verlet)
    for step = 2:Nt
        tNow = t(step);

        F = zeros(1, N);
        for i = 1:N
            if i > 1
                delta_L = x(i) - x(i-1) - diameter;
                F(i) = F(i) - k_spring * delta_L;
            end
            if i < N
                delta_R = x(i+1) - x(i) - diameter;
                F(i) = F(i) + k_spring * delta_R;
            end
            F(i) = F(i) - gamma * v(i);
        end

        a_new = F / mass;
        x_new = x + v*dt + 0.5*a.*dt^2;
        v_new = v + 0.5*(a + a_new)*dt;

        x_new(1) = x0(1) + A_drive * sin(omega * tNow);
        v_new(1) = A_drive * omega * cos(omega * tNow);
        x_new(N) = x0(N);
        v_new(N) = 0;

        x = x_new;  v = v_new;  a = a_new;
        X(:, step) = x;
        V(:, step) = v;
    end

    U = X - x0';

    %% Detect wavefronts
    v_thresh = 0.01 * max(abs(V(:)));

    fwd_time = NaN(1, N);
    ref_time = NaN(1, N);

    for i = 2:N
        fwd_idx = find(V(i,:) >  v_thresh, 1, 'first');
        ref_idx = find(V(i,:) < -v_thresh, 1, 'first');
        if ~isempty(fwd_idx), fwd_time(i) = t(fwd_idx); end
        if ~isempty(ref_idx), ref_time(i) = t(ref_idx); end
    end

    %% Static plots
    figure(1);
    imagesc(t, 1:N, U); axis xy;
    colormap(redblue()); colorbar;
    xlabel('Time [s]'); ylabel('Particle index');
    title(sprintf('Displacement  (N=%d, f=%.2f Hz, m=%.2g, k=%.2g)', ...
                  N, f_drive, mass, k_spring));
    hold on;
    plot(fwd_time, 1:N, 'r.', 'MarkerSize', 6);
    plot(ref_time, 1:N, 'b.', 'MarkerSize', 6);
    legend('Fwd front','Ref front','Location','northeast');
    hold off;

    figure(2);
    plot(t, U(round(N/2), :));
    xlabel('Time [s]'); ylabel('Displacement [m]');
    title(sprintf('Midchain particle %d', round(N/2)));

    %% Animation
    anim_skip = 20;
    r     = diameter/2;
    theta = linspace(0, 2*pi, 40);
    cx = r*cos(theta);  cy = r*sin(theta);

    fig = figure(3);
    set(fig, 'Color', 'k', 'Name', 'Chain Animation');
    ax = axes('Parent', fig, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    xlim(ax, [x0(1)-diameter, x0(end)+diameter]);
    ylim(ax, [-diameter*1.5, diameter*1.5]);
    axis(ax, 'equal'); hold(ax, 'on'); axis(ax, 'off');

    patches = gobjects(N,1);
    for i = 1:N
        patches(i) = fill(ax, x0(i)+cx, cy, [0.4 0.6 1.0], ...
                          'EdgeColor','w','LineWidth',0.8);
    end

    wall_x = x0(N) + r;
    fill(ax, wall_x + [0 0.3 0.3 0]*diameter, [-1 -1 1 1]*diameter, ...
         [0.6 0.6 0.6], 'EdgeColor','none');

    fwd_dot = plot(ax, -999, 0, 'o', 'MarkerSize', 12, ...
                   'MarkerFaceColor','r', 'MarkerEdgeColor','w', 'LineWidth', 1.5);
    ref_dot = plot(ax, -999, 0, 'o', 'MarkerSize', 12, ...
                   'MarkerFaceColor','b', 'MarkerEdgeColor','w', 'LineWidth', 1.5);

    time_txt = text(ax, x0(1), diameter*1.3, '', 'Color','w','FontSize',12);
    title(ax, sprintf('N=%d  f=%.2f Hz  m=%.2g  k=%.2g  ● fwd  ● ref', ...
                      N, f_drive, mass, k_spring), 'Color','w','FontSize',11);

    for step = 1:anim_skip:Nt
        tNow = t(step);

        for i = 1:N
            patches(i).XData = X(i, step) + cx;
        end

        fwd_reached = find(fwd_time <= tNow, 1, 'last');
        if ~isempty(fwd_reached)
            fwd_dot.XData = x0(fwd_reached);
        else
            fwd_dot.XData = -999;
        end

        ref_reached = find(ref_time <= tNow, 1, 'first');
        if ~isempty(ref_reached)
            ref_dot.XData = x0(ref_reached);
        else
            ref_dot.XData = -999;
        end

        time_txt.String = sprintf('t = %.3f s', tNow);
        drawnow;
    end

    if nargout > 0, varargout{1} = t;  end
    if nargout > 1, varargout{2} = U;  end
    if nargout > 2, varargout{3} = x0; end
end

function cmap = redblue()
    n = 128;
    b = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)];
    r = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];
    cmap = [b; r];
end


