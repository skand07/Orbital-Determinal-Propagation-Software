function [r0, v0, oe0, r, v, oe] = OrbitComp(lat, lst, alt, ra, dec, JD, JD_prop)

    %Given constants:
    mu = 398600.4354; %km^3/s^2
    J2 = 1.08263e-3;   % J2 perturbation coefficient
    R_earth = 6378.137; % Earth's radius (km)
    CD = 1.28;
    Am_ratio = 0.0123e3; 
    wE = 7.2921159E-5;% Earth's inertial rotation rate

    %r0, v0, oe0:
    [r0,v0] = Gauss_Sohan(lat,lst,alt,ra,dec,JD,JD_prop);
    [a, e, i, Omega, omega, f] = orbitalElements(r0,v0,mu);
    oe0 = [a;
           e;
           i;
           Omega;
           omega;
           f];

    %state vector:
    state0 = [r0',v0'];

    %propagation:
    t_total = (JD_prop - JD(2))*86400;
    tt_span = [0 t_total];
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    [t, Y] = ode45(@(t, y) orbitalDynamicsJ2(t, y, mu,'yes'), tt_span, state0, options);

    r = Y(end,1:3);
    v = Y(end,4:6);
    r = r';
    v = v';

    [a_prop, e_prop, i_prop, Omega_prop, omega_prop, f_prop] = orbitalElements(r,v,mu);
    oe = [a_prop;
          e_prop;
          i_prop;
          Omega_prop;
          omega_prop;
          f_prop];

        figure('Name','Raw Orbit','Color','w');
    hold on; grid on; axis equal; view(3);
    
    load('earth_topo.mat','topo');  % Must be accessible in your MATLAB path
    R_earth = 6378.137;             % Earth radius [km]
    
    %--- Plot Earth ---
    [sx, sy, sz] = sphere(50);
    sx = R_earth*sx; sy = R_earth*sy; sz = R_earth*sz;
    hEarth = surf(sx, sy, sz);
    hEarth.CData = topo; 
    hEarth.FaceColor = 'texturemap';
    hEarth.EdgeColor = 'none';
    hEarth.FaceLighting = 'gouraud';
    hEarth.SpecularStrength = 0.4;
    hEarth.FaceAlpha = 0.8;
    
    %--- Plot Full Orbit (all points) ---
    r_trajectory_raw = Y(:,1:3);  % Nx3
    plot3(r_trajectory_raw(:,1), ...
          r_trajectory_raw(:,2), ...
          r_trajectory_raw(:,3), ...
          'r-', 'LineWidth',1.5, ...
          'DisplayName','Orbit Path');
    
    %--- Mark initial & final states ---
    plot3(r0(1), r0(2), r0(3), ...
          'bo','MarkerFaceColor','b','MarkerSize',8, ...
          'DisplayName','Initial Pos');
    plot3(r(1), r(2), r(3), ...
          'go','MarkerFaceColor','g','MarkerSize',8, ...
          'DisplayName','Final Pos');
    
    %--- (Optional) Velocity Arrows ---
    vScale = 1000;  % scale factor
    quiver3(r0(1), r0(2), r0(3), ...
            vScale*v0(1), vScale*v0(2), vScale*v0(3), ...
            'b','LineWidth',1.5, 'MaxHeadSize',2, ...
            'DisplayName','Initial Vel');
    quiver3(r(1),  r(2),  r(3), ...
            vScale*v(1),  vScale*v(2),  vScale*v(3), ...
            'g','LineWidth',1.5, 'MaxHeadSize',2, ...
            'DisplayName','Final Vel');
    
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('Satellite Orbit around Earth (Raw Data, No Downsampling)');
    legend('Location','best');
    hold off;
    
    figure('Name','Downsampled Orbit','Color','w');
    hold on; grid on; axis equal; view(3);
    
    % Plot Earth
    [sx, sy, sz] = sphere(50);
    sx = R_earth*sx;  sy = R_earth*sy;  sz = R_earth*sz;
    hEarthDown = surf(sx, sy, sz);
    hEarthDown.CData = topo;
    hEarthDown.FaceColor = 'texturemap';
    hEarthDown.EdgeColor = 'none';
    hEarthDown.FaceLighting = 'gouraud';
    hEarthDown.SpecularStrength = 0.4;
    hEarth.FaceAlpha = 0.8;
    
    % Downsample the orbit
    stepSize = 1000;
    idx = 1:stepSize:size(Y,1);
    Y_down = Y(idx,:);
    
    % Plot downsampled orbit
    plot3(Y_down(:,1), Y_down(:,2), Y_down(:,3), ...
          'r-o','LineWidth',1.2, 'MarkerSize',3, ...
          'MarkerFaceColor','r', 'DisplayName','Orbit (Downsampled)');
    
    % Mark initial/final
    plot3(r0(1), r0(2), r0(3), ...
          'bo','MarkerFaceColor','b','MarkerSize',8, ...
          'DisplayName','Initial Pos');
    plot3(r(1),  r(2),  r(3), ...
          'go','MarkerFaceColor','g','MarkerSize',8, ...
          'DisplayName','Final Pos');
    
    % Velocity arrows
    quiver3(r0(1), r0(2), r0(3), ...
            vScale*v0(1), vScale*v0(2), vScale*v0(3), ...
            'b','LineWidth',1.5,'MaxHeadSize',2, ...
            'DisplayName','Initial Vel');
    quiver3(r(1), r(2), r(3), ...
            vScale*v(1), vScale*v(2), vScale*v(3), ...
            'g','LineWidth',1.5,'MaxHeadSize',2, ...
            'DisplayName','Final Vel');
    
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('Satellite Orbit (Downsampled by 1000)');
    legend('Location','best');
    hold off;
    
        % --- Plot Orbital Elements and Derivatives ---
    plotOrbitalElementsOnly(t, Y, mu);
    
    % --- Compute Satellite Revolutions ---
    r_vectors = Y(:,1:3);
    num_points = size(Y,1);
    revolution_indices = [];  % To store indices where true anomaly crosses 0
    
    mu = 398600.4354;
    
    % Compute true anomaly at each time
    f_deg = zeros(num_points,1);
    for k = 1:num_points
        r_k = Y(k,1:3)';
        v_k = Y(k,4:6)';
        [~, ~, ~, ~, ~, f] = orbitalElements(r_k, v_k, mu);
        f_deg(k) = f;
    end
    
end






