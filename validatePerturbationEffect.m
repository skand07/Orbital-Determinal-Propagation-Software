function validatePerturbationEffect(t, state, mu, time_unit)
    % Inputs:
    %   t         - Time vector (s)
    %   state     - Nx6 matrix containing [x, y, z, vx, vy, vz] at each time step
    %   mu        - Gravitational parameter (km^3/s^2)
    %   time_unit - 'minutes' or 'hours' to control time axis scaling
    % Convert time from seconds to minutes for better visualization

    % Convert time based on user input
    if strcmp(time_unit, 'minutes')
        t_scaled = t / 60; % Convert seconds to minutes
        time_label = 'Time (min)';
    elseif strcmp(time_unit, 'hours')
        t_scaled = t / 3600; % Convert seconds to hours
        time_label = 'Time (hours)';
    else
        error('Invalid time unit. Choose "minutes" or "hours".');
    end
    
    % Number of time steps
    N = length(t);

    % Initialize storage for orbital elements
    a = zeros(N,1);
    e = zeros(N,1);
    i = zeros(N,1);
    RAAN = zeros(N,1);
    omega = zeros(N,1);
    f = zeros(N,1);

    % Extract orbital elements at each timestep
    for k = 1:N
        r = state(k,1:3)';
        v = state(k,4:6)';
        
        try
            koe = rv2koe(r, v, mu, 'deg');
        catch
            warning('rv2koe failed at timestep %d', k);
            continue;
        end
        
        a(k) = koe(1);
        e(k) = koe(2);
        i(k) = koe(3); % Degrees
        RAAN(k) = koe(4); % Degrees
        omega(k) = koe(5); % Degrees
        f(k) = koe(6); % Degrees
    end

  
    % Define how many plots per row (to keep spacing clean)
    
    figure;
    
    subplot(3,2,1);
    plot(t_scaled, a, 'LineWidth', 1.5);
    xlabel(time_label); ylabel('Semi-major axis (km)');
    title('Semi-major axis (a)'); grid on;
    xlim([min(t_scaled), max(t_scaled)]);
    ylim([min(a)-5, max(a)+5]);

    subplot(3,2,2);
    plot(t_scaled, e, 'LineWidth', 1.5);
    xlabel(time_label); ylabel('Eccentricity');
    title('Eccentricity (e)'); grid on;
    xlim([min(t_scaled), max(t_scaled)]);
    ylim([min(e)-0.0001, max(e)+0.0001]);

    subplot(3,2,3);
    plot(t_scaled, i, 'LineWidth', 1.5);
    xlabel(time_label); ylabel('Inclination (deg)');
    title('Inclination (i)'); grid on;
    xlim([min(t_scaled), max(t_scaled)]);
    ylim([min(i)-0.05, max(i)+0.05]);

    subplot(3,2,4);
    plot(t_scaled, RAAN, 'LineWidth', 1.5);
    xlabel(time_label); ylabel('RAAN (deg)');
    title('Right Ascension of Ascending Node (Ω)'); grid on;
    xlim([min(t_scaled), max(t_scaled)]);
    ylim([min(RAAN)-0.1, max(RAAN)+0.1]);

    subplot(3,2,5);
    plot(t_scaled, omega, 'LineWidth', 1.5);
    xlabel(time_label); ylabel('Argument of Perigee (deg)');
    title('Argument of Perigee (ω)'); grid on;
    xlim([min(t_scaled), max(t_scaled)]);
    ylim([min(omega)-0.1, max(omega)+0.1]);

    subplot(3,2,6);
    plot(t_scaled, f, 'LineWidth', 1.5);
    xlabel(time_label); ylabel('True Anomaly (deg)');
    title('True Anomaly (f)'); grid on;
    xlim([min(t_scaled), max(t_scaled)]);
    ylim([min(f)-5, max(f)+5]);

    % Enable zooming for detailed inspection
    zoom on;
    disp('Use mouse scroll or drag to zoom in on the plots.');
end