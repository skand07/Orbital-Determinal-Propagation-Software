function plotOrbitalElementsOnly(t, Y, mu)
    % Clean, polished plot of orbital elements over time in hours

    % Preallocate
    N = length(t);
    a = zeros(N,1); e = zeros(N,1); i = zeros(N,1);
    RAAN = zeros(N,1); omega = zeros(N,1); f = zeros(N,1);

    % Extract orbital elements
    for k = 1:N
        r = Y(k,1:3)';
        v = Y(k,4:6)';
        [a(k), e(k), i(k), RAAN(k), omega(k), f(k)] = orbitalElements(r, v, mu);
    end

    % Time in hours
    t_hours = t / 3600;

    % --- Plot Settings ---
    font_title = 15; 
    font_labels = 13;
    line_w = 2;

    figure('Name','Orbital Elements vs Time', 'Color','w', 'Position', [100 100 1200 800]);

    % --- Plot 1: Semi-major Axis ---
    subplot(3,2,1);
    plot(t_hours, a, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Semi-major axis (km)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Semi-major axis (a)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % --- Plot 2: Eccentricity ---
    subplot(3,2,2);
    plot(t_hours, e, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Eccentricity', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Eccentricity (e)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % --- Plot 3: Inclination ---
    subplot(3,2,3);
    plot(t_hours, i, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Inclination (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Inclination (i)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % --- Plot 4: RAAN ---
    subplot(3,2,4);
    plot(t_hours, RAAN, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('RAAN (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Right Ascension of Ascending Node (\Omega)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % --- Plot 5: Argument of Perigee ---
    subplot(3,2,5);
    plot(t_hours, omega, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('Arg. of Perigee (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('Argument of Perigee (\omega)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % --- Plot 6: True Anomaly ---
    subplot(3,2,6);
    plot(t_hours, f, 'b', 'LineWidth', line_w);
    xlabel('Time (hours)', 'FontSize', font_labels, 'FontWeight', 'bold');
    ylabel('True Anomaly (deg)', 'FontSize', font_labels, 'FontWeight', 'bold');
    title('True Anomaly (f)', 'FontSize', font_title, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % Overall Title
    sgtitle('Orbital Elements vs Time', 'FontSize', font_title+2, 'FontWeight', 'bold');
end
