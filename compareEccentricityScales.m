function compareEccentricityScales(t, Y, mu)
    % Extract eccentricity over time
    N = length(t);
    e = zeros(N,1);

    for k = 1:N
        r = Y(k,1:3)';
        v = Y(k,4:6)';
        [~, e_k, ~, ~, ~, ~] = orbitalElements(r, v, mu);
        e(k) = e_k;
    end

    t_hours = t / 3600;

    % Plot side-by-side comparison
    figure('Name','Eccentricity Linear vs Log Scale', 'Color','w', 'Position',[100 100 1000 400]);

    % --- Linear Scale ---
    subplot(1,2,1);
    plot(t_hours, e, 'b', 'LineWidth', 2);
    xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Eccentricity', 'FontSize', 12, 'FontWeight', 'bold');
    title('Eccentricity (Linear Scale)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    % --- Log Scale ---
    subplot(1,2,2);
    semilogy(t_hours, e, 'r', 'LineWidth', 2);
    xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Eccentricity (log)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Eccentricity (Log Scale)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    xlim([0 max(t_hours)]);
    ylim padded;

    sgtitle('Comparison of Eccentricity on Linear vs Log Scale', 'FontSize', 15, 'FontWeight', 'bold');
end
