% Team K6
clear all; clc; format long;

%% Site Parameters
site_data  = [32.8811913; 
             -117.2336137; 
              0.111];
lat_s = site_data(1);
lon_s = site_data(2);
alt_s = site_data(3);
JD_prop = 2460542.7792;

%% Given Observation Data:
JD_array = [2460259.470365230; 
            2460259.470712452; 
            2460259.471059674];
RA_array = [274.943444732592; 
            288.272238630768; 
            301.895673996418];
Dec_array = [-34.275778246781; 
            -32.355828886996; 
            -29.159102805200];
LST_array = [281.953273625282; 
             282.078615856645; 
             282.203958087591];
[r0, v0, oe0, r, v, oe] = OrbitComp(lat_s, LST_array, alt_s, RA_array, Dec_array, JD_array, JD_prop);
[r0_oblate, v0_oblate] = Gauss_Oblate(lat_s, LST_array, alt_s, RA_array, Dec_array, JD_array, JD_prop);
delta_r = norm(r0 - r0_oblate);


%% After OrbitComp is called:
fprintf('Results from OrbitComp:\n\n');

% Display initial position (r0) and velocity (v0):
fprintf('Initial Position Vector (r0) [km]:\n');
disp(r0);

fprintf('Initial Velocity Vector (v0) [km/s]:\n');
disp(v0);


% Display each orbital element separately:
fprintf('Initial Semi-major axis (a0): %.4f km\n', oe0(1));
fprintf('Initial Eccentricity (e0):    %.6f\n', oe0(2));
fprintf('Initial Inclination (i0):     %.4f deg\n', oe0(3));
fprintf('Initial RAAN (Omega0):        %.4f deg\n', oe0(4));
fprintf('Initial Arg. of Perigee (omega0): %.4f deg\n', oe0(5));
fprintf('Initial True Anomaly (f0):    %.4f deg\n', oe0(6));

%% Display final (post-propagation) results
fprintf('\n Propagated Orbit Parameters \n\n');

fprintf('Propagated Position Vector (r) [km]:\n');
disp(r);

fprintf('Propagated Velocity Vector (v) [km/s]:\n');
disp(v);

fprintf('Propagated Semi-major axis (a):        %.4f km\n', oe(1));
fprintf('Propagated Eccentricity (e):           %.6f\n', oe(2));
fprintf('Propagated Inclination (i):            %.4f deg\n', oe(3));
fprintf('Propagated RAAN (Omega):               %.4f deg\n', oe(4));
fprintf('Propagated Argument of Perigee (omega): %.4f deg\n', oe(5));
fprintf('Propagated True Anomaly (f):           %.4f deg\n', oe(6));


fprintf('\nDifference in position (oblate vs spherical): %.4f km\n', delta_r);

disp('Running K6-S6 approach analysis and avoidance:');
compare();
