function plotOrbitalElementsOnlyDownsampled(t, Y, mu, downsampleFactor)
    % Plots individual orbital elements over time using downsampled data
    % Inputs:
    %   t - time vector (s)
    %   Y - state history from ODE45
    %   mu - gravitational parameter (km^3/s^2)
    %   downsampleFactor - e.g., 10 means every 10th point plotted

    % Downsample time and states
    t_down = t(1:downsampleFactor:end);
    Y_down = Y(1:downsampleFactor:end,:);

    N = length(t_down);
    a = zeros(N,1); e = zeros(N,1); i = zeros(N,1);
    RAAN = zeros(N,1); omega = zeros(N,1); f = zeros(N,1);

    % Extract orbital elements from downsampled data
    for k = 1:N
        r = Y_down(k,1:3)';
        v = Y_down(k,4:6)';
        [a(k), e(k), i(k), RAAN(k), omega(k), f(k)] = orbitalElements(r, v, mu);
    end

    % Time in hours
    t_hours = t_down / 3600;

    % Plot settings
    font_title = 14; 
    font_labels = 12;
    line_w = 2;

    % --- Plot Each Element Individually ---

    % Semi-major axis
    figure('Name','Semi-major Axis','Color','w');
    plot(t_hours, a, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Semi-major axis (km)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Semi-major axis (a)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]); ylim padded;

    % Eccentricity
    figure('Name','Eccentricity','Color','w');
    plot(t_hours, e, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Eccentricity', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Eccentricity (e)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]); ylim padded;

    % Inclination
    figure('Name','Inclination','Color','w');
    plot(t_hours, i, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Inclination (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Inclination (i)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]); ylim padded;

    % RAAN
    figure('Name','RAAN','Color','w');
    plot(t_hours, RAAN, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('RAAN (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Right Ascension of Ascending Node (\Omega)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]); ylim padded;

    % Argument of Perigee
    figure('Name','Argument of Perigee','Color','w');
    plot(t_hours, omega, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Arg. of Perigee (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Argument of Perigee (\omega)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]); ylim padded;

    % True Anomaly
    figure('Name','True Anomaly','Color','w');
    plot(t_hours, f, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('True Anomaly (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('True Anomaly (f)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]); ylim padded;
end
