function plotOrbitalElementsAndRates(t, Y, mu, time_unit)
    % Compute orbital elements and their rates over time
    N = length(t);
    a = zeros(N,1); e = zeros(N,1); i = zeros(N,1);
    RAAN = zeros(N,1); omega = zeros(N,1); f = zeros(N,1);

    % Extract orbital elements
    for k = 1:N
        r = Y(k,1:3)';
        v = Y(k,4:6)';
        [a(k), e(k), i(k), RAAN(k), omega(k), f(k)] = orbitalElements(r, v, mu);
    end

    % Time conversion
    switch time_unit
        case 'minutes'
            t_scaled = t / 60;
            t_label = 'Time (min)';
        case 'hours'
            t_scaled = t / 3600;
            t_label = 'Time (hours)';
        otherwise
            error('Invalid time unit. Choose "minutes" or "hours".');
    end

    % Numerical Derivatives (finite differences)
    dt = diff(t_scaled);
    da_dt = diff(a)./dt; de_dt = diff(e)./dt;
    di_dt = diff(i)./dt; dRAAN_dt = diff(RAAN)./dt;
    domega_dt = diff(omega)./dt; df_dt = diff(f)./dt;
    t_mid = t_scaled(1:end-1) + dt/2;  % Midpoints for derivative plots

    % Plot each element and its derivative
    figure('Name','Orbital Elements and Derivatives','Color','w');
    subplot(3,2,1);
    plot(t_scaled,a,'b','LineWidth',1.5); hold on;
    plot(t_mid,da_dt,'r','LineWidth',1.2);
    xlabel(t_label); ylabel('a [km] / da/dt'); title('Semi-major Axis'); grid on;

    subplot(3,2,2);
    plot(t_scaled,e,'b','LineWidth',1.5); hold on;
    plot(t_mid,de_dt,'r','LineWidth',1.2);
    xlabel(t_label); ylabel('e / de/dt'); title('Eccentricity'); grid on;

    subplot(3,2,3);
    plot(t_scaled,i,'b','LineWidth',1.5); hold on;
    plot(t_mid,di_dt,'r','LineWidth',1.2);
    xlabel(t_label); ylabel('i [deg] / di/dt'); title('Inclination'); grid on;

    subplot(3,2,4);
    plot(t_scaled,RAAN,'b','LineWidth',1.5); hold on;
    plot(t_mid,dRAAN_dt,'r','LineWidth',1.2);
    xlabel(t_label); ylabel('\Omega [deg] / d\Omega/dt'); title('RAAN'); grid on;

    subplot(3,2,5);
    plot(t_scaled,omega,'b','LineWidth',1.5); hold on;
    plot(t_mid,domega_dt,'r','LineWidth',1.2);
    xlabel(t_label); ylabel('\omega [deg] / d\omega/dt'); title('Arg. of Perigee'); grid on;

    subplot(3,2,6);
    plot(t_scaled,f,'b','LineWidth',1.5); hold on;
    plot(t_mid,df_dt,'r','LineWidth',1.2);
    xlabel(t_label); ylabel('f [deg] / df/dt'); title('True Anomaly'); grid on;

    legend('Value','Derivative');
end
