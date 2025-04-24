function [r0, v0] = Gauss_Oblate(lat_s, LST_array, alt_s, RA_array, Dec_array, JD_array, JD_prop)
% Gauss_Oblate - Angles-only orbit determination 

    mu = 398600.4354; % km^3/s^2

    % Oblate Earth model
    a_earth = 6378.137;               % km
    f = 1/298.257223563;
    e2 = 2*f - f^2;

    R_site_matrix = zeros(3,3);
    for i = 1:3
        lon_rad = deg2rad(LST_array(i));  % Convert LST (deg) to rad
        lat_rad = deg2rad(lat_s);        % Geodetic lat (deg) to rad

        N = a_earth / sqrt(1 - e2 * sin(lat_rad)^2);
        x = (N + alt_s) * cos(lat_rad) * cos(lon_rad);
        y = (N + alt_s) * cos(lat_rad) * sin(lon_rad);
        z = ((1 - e2) * N + alt_s) * sin(lat_rad);
        R_site_matrix(:, i) = [x; y; z];
    end

    rho_hat_matrix = zeros(3,3);
    for i = 1:3
        rho_hat_matrix(:, i) = [cosd(Dec_array(i)) * cosd(RA_array(i));
                                cosd(Dec_array(i)) * sind(RA_array(i));
                                sind(Dec_array(i))];
    end

    Tau1 = (JD_array(1) - JD_array(2)) * 86400;
    Tau3 = (JD_array(3) - JD_array(2)) * 86400;
    Tau  = Tau3 - Tau1;

    p1 = cross(rho_hat_matrix(:,2), rho_hat_matrix(:,3));
    p2 = cross(rho_hat_matrix(:,1), rho_hat_matrix(:,3));
    p3 = cross(rho_hat_matrix(:,1), rho_hat_matrix(:,2));

    D0  = dot(rho_hat_matrix(:,1), p1);
    D11 = dot(R_site_matrix(:,1), p1);
    D21 = dot(R_site_matrix(:,2), p1);
    D31 = dot(R_site_matrix(:,3), p1);

    D12 = dot(R_site_matrix(:,1), p2);
    D22 = dot(R_site_matrix(:,2), p2);
    D32 = dot(R_site_matrix(:,3), p2);

    D13 = dot(R_site_matrix(:,1), p3);
    D23 = dot(R_site_matrix(:,2), p3);
    D33 = dot(R_site_matrix(:,3), p3);

    A = (-D12 * (Tau3 / Tau) + D22 + D32 * (Tau1 / Tau)) / D0;
    B = (1 / (6 * D0)) * (D12 * (Tau3^2 - Tau^2) * (Tau3 / Tau) + ...
                          D32 * (Tau^2 - Tau1^2) * (Tau1 / Tau));
    E = dot(R_site_matrix(:,2), rho_hat_matrix(:,2));
    R2_mag_sq = norm(R_site_matrix(:,2))^2;

    a_poly = -(A^2 + 2*A*E + R2_mag_sq);
    b_poly = -2 * mu * B * (A + E);
    c_poly = -mu^2 * B^2;

    r2_guess = 7000; tol = 1e-8;
    for iter = 1:100
        F  = r2_guess^8 + a_poly*r2_guess^6 + b_poly*r2_guess^3 + c_poly;
        Fp = 8*r2_guess^7 + 6*a_poly*r2_guess^5 + 3*b_poly*r2_guess^2;
        r2_new = r2_guess - F / Fp;
        if abs(r2_new - r2_guess) < tol
            break;
        end
        r2_guess = r2_new;
    end
    r2_mag = r2_new;

    %Slant ranges
    rho2 = A + mu*B / r2_mag^3;
    rho1 = (1/D0) * ( ...
        (6*((D31*(Tau1/Tau3) + D21*(Tau/Tau3))*r2_mag^3) + ...
        mu*D31*(Tau^2 - Tau1^2)*(Tau1/Tau3)) / ...
        (6*r2_mag^3 + mu*(Tau^2 - Tau3^2)) - D11);

    rho3 = (1/D0) * ( ...
        (6*((D13*(Tau3/Tau1) - D23*(Tau/Tau1))*r2_mag^3) + ...
        mu*D13*(Tau^2 - Tau3^2)*(Tau3/Tau1)) / ...
        (6*r2_mag^3 + mu*(Tau^2 - Tau1^2)) - D33);

    %Position
    r1 = R_site_matrix(:,1) + rho1 * rho_hat_matrix(:,1);
    r2 = R_site_matrix(:,2) + rho2 * rho_hat_matrix(:,2);
    r3 = R_site_matrix(:,3) + rho3 * rho_hat_matrix(:,3);

    %Velocity 
    f1 = 1 - (mu * Tau1^2) / (2 * r2_mag^3);
    g1 = Tau1 - (mu * Tau1^3) / (6 * r2_mag^3);
    f3 = 1 - (mu * Tau3^2) / (2 * r2_mag^3);
    g3 = Tau3 - (mu * Tau3^3) / (6 * r2_mag^3);

    denom = f1 * g3 - f3 * g1;
    v2 = ( -f3 * r1 + f1 * r3 ) / denom;

    r0 = r2;
    v0 = v2;
end
