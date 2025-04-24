function [r0, v0] = Gauss_Sohan(lat, lst, alt, ra, dec, JD, JD_prop)
%GAUSS_SOHAN  Gauss angles-only orbit determination, accepting 3x1 arrays
%             for LST, RA, DEC, JD.
%
% Inputs:
%   lat  - Station latitude (scalar, deg)
%   lst  - 3x1 array of local sidereal times (deg)
%   alt  - Station altitude above Earth's surface (scalar, km)
%   ra   - 3x1 array of right ascensions (deg)
%   dec  - 3x1 array of declinations (deg)
%   JD   - 3x1 array of Julian dates
%   JD_prop - (Not explicitly used but included for consistency)
%
% Outputs:
%   r0, v0, - Position, velocity, at the
%                 second (middle) observation

    % --- Constants ---
    R_earth = 6378.137;      % Earth's equatorial radius [km]
    mu      = 3.98600442e5;  % Earth GM [km^3/s^2]
    
    % --- Station parameters ---
    lat_s  = lat;       % scalar latitude [deg]
    alt_s  = alt;       % scalar altitude [km]
    
    % --- Extract the 3×1 inputs ---
    %     Each of these is 3×1, e.g. lst(1), lst(2), lst(3).
    LST_array = lst;   % 3x1, in degrees
    RA_array  = ra;    % 3x1, in degrees
    Dec_array = dec;   % 3x1, in degrees
    JD_array  = JD;    % 3x1

    % --- Compute station position vectors for each observation ---
    %     (R_earth + alt_s) in spherical coords with (lat_s, LST).
    R_site_matrix = zeros(3,3);
    for i = 1:3
        R_site_matrix(:, i) = (R_earth + alt_s) * [
            cosd(lat_s) * cosd(LST_array(i));
            cosd(lat_s) * sind(LST_array(i));
            sind(lat_s)
        ];
    end
    
    % --- Line-of-sight unit vectors from RA, Dec (3x1) ---
    rho_hats_matrix = zeros(3,3);
    for i = 1:3
        rho_hats_matrix(:, i) = [
            cosd(Dec_array(i)) * cosd(RA_array(i));
            cosd(Dec_array(i)) * sind(RA_array(i));
            sind(Dec_array(i))
        ];
    end
    
    % --- Time differences (seconds) ---
    Tau_1 = (JD_array(1) - JD_array(2)) * 86400;  
    Tau_3 = (JD_array(3) - JD_array(2)) * 86400;
    Tau   = Tau_3 - Tau_1;
    
    % --- Scalar triple products used by Gauss's method ---
    p1 = cross(rho_hats_matrix(:,2), rho_hats_matrix(:,3));
    p2 = cross(rho_hats_matrix(:,1), rho_hats_matrix(:,3));
    p3 = cross(rho_hats_matrix(:,1), rho_hats_matrix(:,2));
    
    D0 = dot(rho_hats_matrix(:,1), p1);

    D11 = dot(R_site_matrix(:,1), p1);
    D21 = dot(R_site_matrix(:,2), p1);
    D31 = dot(R_site_matrix(:,3), p1);

    D12 = dot(R_site_matrix(:,1), p2);
    D22 = dot(R_site_matrix(:,2), p2);
    D32 = dot(R_site_matrix(:,3), p2);

    D13 = dot(R_site_matrix(:,1), p3);
    D23 = dot(R_site_matrix(:,2), p3);
    D33 = dot(R_site_matrix(:,3), p3);

    % --- A, B, E, etc. for polynomial in r2 ---
    A = ( -D12*(Tau_3/Tau) + D22 + D32*(Tau_1/Tau) ) / D0;
    B = (1/(6*D0)) * ( ...
          D12 * (Tau_3^2 - Tau^2) * (Tau_3/Tau) + ...
          D32 * (Tau^2 - Tau_1^2) * (Tau_1/Tau) ...
        );
    E     = dot(R_site_matrix(:,2), rho_hats_matrix(:,2));
    R2_sq = norm(R_site_matrix(:,2))^2;

    a_poly = -(A^2 + 2*A*E + R2_sq);
    b_poly = -2*mu * B * (A + E);
    c_poly = -mu^2 * B^2;

    % --- Solve for r2 (range at 2nd observation) via Newton's method ---
    r2_guess       = 9000;         % initial guess [km]
    tolerance      = 1e-6;
    max_iterations = 100;
    
    for iter = 1:max_iterations
        F_r2       = r2_guess^8 + a_poly*r2_guess^6 + b_poly*r2_guess^3 + c_poly;
        F_prime_r2 = 8*r2_guess^7 + 6*a_poly*r2_guess^5 + 3*b_poly*r2_guess^2;
        r2_new     = r2_guess - (F_r2 / F_prime_r2);

        if abs(r2_new - r2_guess)/r2_guess < tolerance
            break
        end
        r2_guess = r2_new;
    end
    
    r2_scalar = double(r2_new);
    
    % --- Compute slant ranges (rho) for each observation ---
    rho2 = A + (mu*B)/(r2_scalar^3);

    %   The expressions for rho1, rho3 come from the same Gauss method steps
    rho1 = (1/D0) * ( ...
           ( 6*( (D31*(Tau_1/Tau_3) + D21*(Tau/Tau_3))*r2_scalar^3 ) ...
             + mu*D31*(Tau^2 - Tau_1^2)*(Tau_1/Tau_3) ) ...
           / ( 6*r2_scalar^3 + mu*(Tau^2 - Tau_3^2) ) ...
           - D11 );
    
    rho3 = (1/D0) * ( ...
           ( 6*( (D13*(Tau_3/Tau_1) - D23*(Tau/Tau_1))*r2_scalar^3 ) ...
             + mu*D13*(Tau^2 - Tau_3^2)*(Tau_3/Tau_1) ) ...
           / ( 6*r2_scalar^3 + mu*(Tau^2 - Tau_1^2) ) ...
           - D33 );

    % --- Construct the geocentric position vectors r1, r2, r3 ---
    r1 = R_site_matrix(:,1) + rho1 * rho_hats_matrix(:,1);
    r2 = R_site_matrix(:,2) + rho2 * rho_hats_matrix(:,2);
    r3 = R_site_matrix(:,3) + rho3 * rho_hats_matrix(:,3);

    % --- f and g functions at Tau_1 and Tau_3 (approx 2-body) ---
    r2_mag = norm(r2);
    f1     = 1 - (mu * Tau_1^2) / (2 * r2_mag^3);
    g1     = Tau_1 - (mu * Tau_1^3) / (6 * r2_mag^3);
    f3     = 1 - (mu * Tau_3^2) / (2 * r2_mag^3);
    g3     = Tau_3 - (mu * Tau_3^3) / (6 * r2_mag^3);

    denominator = f1*g3 - f3*g1;
    v2 = (1/denominator) * ( -f3*r1 + f1*r3 );

    % --- Return the Gauss "solution" at the second (middle) observation ---
    r0 = r2;  
    v0 = v2;

end
