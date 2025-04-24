function [a, e, i, Omega, omega, f] = orbitalElements(r_vec, v_vec, mu)
%ORBITALELEMENTS Compute classical orbital elements from r, v, mu, 
% and display them with appropriate significant figures.
%
% Inputs:
%   r_vec  : Position vector [3 x 1] (e.g., in km)
%   v_vec  : Velocity vector [3 x 1] (e.g., in km/s)
%   mu     : Gravitational parameter (e.g., km^3/s^2)
%
% Outputs:
%   a      : Semi-major axis [km]
%   e      : Eccentricity
%   i      : Inclination [deg]
%   Omega  : Right ascension of ascending node (RAAN) [deg]
%   omega  : Argument of perigee [deg]
%   f      : True anomaly [deg]

    %-- 1. Magnitudes of r and v
    r = norm(r_vec);
    v = norm(v_vec);

    %-- 2. Angular momentum vector and its magnitude
    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);

    %-- 3. Eccentricity vector and scalar eccentricity
    e_vec = (1/mu) * ( cross(v_vec, h_vec) - mu*(r_vec/r) );
    e = norm(e_vec);

    %-- 4. Node vector (points to ascending node)
    k_vec = [0; 0; 1];         % Reference z-axis
    N_vec = cross(k_vec, h_vec);
    N = norm(N_vec);

    %-- 5. Specific orbital energy and semi-major axis
    energy = 0.5*(v^2) - mu/r;   % Specific orbital energy
    a = -mu / (2*energy);        % Semi-major axis

    %-- 6. Inclination (degrees)
    i = acosd(h_vec(3) / h);

    %-- 7. RAAN (Omega), in degrees
    if N_vec(2) >= 0
        Omega = acosd(N_vec(1)/N);
    else
        Omega = 360 - acosd(N_vec(1)/N);
    end

    %-- 8. Argument of Perigee (omega), in degrees
    if e_vec(3) >= 0
        omega = acosd(dot(N_vec, e_vec)/(N*e));
    else
        omega = 360 - acosd(dot(N_vec, e_vec)/(N*e));
    end

    %-- 9. True Anomaly (f), in degrees
    if dot(r_vec, v_vec) >= 0
        f = acosd(dot(e_vec, r_vec)/(e*r));
    else
        f = 360 - acosd(dot(e_vec, r_vec)/(e*r));
    end

end

