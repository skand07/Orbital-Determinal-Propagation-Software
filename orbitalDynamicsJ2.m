function dydt = orbitalDynamicsJ2(~, y, mu, Drag_switch)
    % Computes the dynamics including J2 perturbation
    % y = [r; v] where r = position vector (km), v = velocity vector (km/s)
    % Drag_swith - either 'yes' or 'no'
    
    % Constants
    J2 = 1.08263e-3;   % J2 perturbation coefficient
    R_earth = 6378.137; % Earth's radius (km)
    CD = 1.28;
    Am_ratio = 0.0123e3; 

    % Extract position and velocity
    r = y(1:3);
    v = y(4:6);
    
    r_mag = norm(r);
    v_mag = norm(v);
    
    % Keplerian acceleration
    a_kepler = -mu / r_mag^3 * r;

    % J2 perturbation acceleration
    z = r(3);
    factor = -1.5 * J2 * mu * R_earth^2 / r_mag^5;
    a_J2 = factor * [(1 - 5 * (z^2 / r_mag^2)) * r(1);
                    (1 - 5 * (z^2 / r_mag^2)) * r(2);
                    (3 - 5 * (z^2 / r_mag^2)) * z];

    % Drag pertubation acceleration
    h = r_mag - R_earth; % satellite's altitude
    [rho_atm, H] = atmosphere(h);

    factor_drag = -0.5*Am_ratio*CD*rho_atm*v_mag;

    a_drag = factor_drag * [v(1); v(2); v(3)];


    if strcmp(Drag_switch, 'yes')
        % Total acceleration
        a_total = a_kepler + a_J2 + a_drag;
    else
        % Total acceleration
        a_total = a_kepler + a_J2;
    end

    %% using orbit frame

    % oe = rv2koe(r,v,mu, 'deg');
    % angle = oe(5) + oe(6);
    % 
    % h_vec = cross(r, v);
    % h_mag = norm(h_vec);
    % 
    % r_hat = r / r_mag; % Radial unit vector
    % h_hat = h_vec / h_mag; % Normal unit vector
    % theta_hat = cross(h_hat, r_hat); % Transverse unit vector
    % 
    % J2 perturbation acceleration components in RTN (orbit frame)
    % z = r(3);
    % factor = -1.5 * J2 * mu * R_earth^2 / r_mag^4;
    % 
    % Compute components in RTN frame
    % a_J2_r = factor * 0.5*(1 - 3*(sind(oe(3))^2)*(sind(angle)^2));
    % a_J2_theta = factor * (sind(oe(3))^2)*(sind(angle)^2)*(cosd(angle));
    % a_J2_h = factor * sind(oe(3))*sind(angle)*cosd(angle);
    % 
    % a_J2 = a_J2_r * r_hat + a_J2_theta * theta_hat + a_J2_h * h_hat; %both RTN and ECI
    % 
    % Total acceleration in ECI
    % a_total = a_kepler + a_J2;

    % Time derivatives
    dydt = [v; a_total]; 
end