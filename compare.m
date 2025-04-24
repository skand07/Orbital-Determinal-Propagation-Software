function compare()
    %% Constants
    mu = 398600.4354;     % km^3/s^2
    tspan = [0 2.4478e+07];  
    options = odeset('RelTol',1e-8,'AbsTol',1e-9);
    lat_s = 32.8801;      % deg
    alt_s = 0.111;        % km

    %% K6 Data
    JD_K6 = [2460259.470365230, 2460259.470712452, 2460259.471059674];
    RA_K6 = [274.943444732592, 288.272238630768, 301.895673996418];
    Dec_K6 = [-34.275778246781, -32.355828886996, -29.159102805200];
    LST_K6 = [281.953273625282, 282.078615856645, 282.203958087591];

    %% S6 Data
    JD_S6 = [2460255.857027431, 2460255.857374653, 2460255.857721875];
    RA_S6 = [88.466025147620, 88.617612938734, 89.023442430480];
    Dec_S6 = [-31.494647078655, -24.858385749194, -15.686618283112];
    LST_S6 = [57.590189090598, 57.715531321961, 57.840873552490];

    %% Compute r0, v0 for K6 and S6
    [r0_K6, v0_K6] = Gauss_Oblate(lat_s, LST_K6, alt_s, RA_K6, Dec_K6, JD_K6, JD_K6(2));
    [r0_S6, v0_S6] = Gauss_Oblate(lat_s, LST_S6, alt_s, RA_S6, Dec_S6, JD_S6, JD_S6(2));

    %% Propagate both satellites
    y0 = [r0_K6; v0_K6; r0_S6; v0_S6];
    [t, Y] = ode45(@(t,y) two_sat_orbit(t,y,mu), tspan, y0, options);

    %% Compute closest approach
    r_K6 = Y(:,1:3);
    r_S6 = Y(:,7:9);
    v_K6 = Y(:,4:6);
    v_S6 = Y(:,10:12);

    separation = vecnorm(r_K6 - r_S6, 2, 2);
    [min_dist, idx] = min(separation);
    t_min = t(idx)/3600;
    fprintf('Closest Approach: %.4f km at %.2f hours\n', min_dist, t_min);

    %% Optimized Δv along relative velocity at closest approach
    rel_v = v_S6(idx,:) - v_K6(idx,:);
    unit_rel_v = rel_v' / norm(rel_v);
    delta_v_mag = 0.02;  % km/s
    delta_v = delta_v_mag * unit_rel_v;

    % Apply impulse to S6 at closest approach
    y_burn = Y(idx,:)';
    y_burn(10:12) = y_burn(10:12) + delta_v;
    [t2, Y2] = ode45(@(t,y) two_sat_orbit(t,y,mu), tspan, y_burn, options);

    %% Separation post-burn:
    r_K6_post = Y2(:,1:3);
    r_S6_post = Y2(:,7:9);
    sep_post = vecnorm(r_K6_post - r_S6_post,2,2);
    
    %% 48-hour Window Around Closest Approach
    start_hr = 1628.57 - 24;
    end_hr   = 1628.57 + 24;
    t_start_sec = start_hr * 3600;
    t_end_sec   = end_hr * 3600;
    
    % Indices for pre-burn and post-burn
    idx_win_pre  = find(t >= t_start_sec & t <= t_end_sec);
    idx_win_post = find(t2 >= t_start_sec & t2 <= t_end_sec);
    
    t_win_hr  = t(idx_win_pre) / 3600;
    t2_win_hr = t2(idx_win_post) / 3600;
    
    sep_win      = separation(idx_win_pre);
    sep_post_win = sep_post(idx_win_post);
    
    %% Plot Focused Separation
    figure; hold on; grid on;
    plot(t_win_hr, sep_win, 'b', 'LineWidth',1.5);
    plot(t2_win_hr, sep_post_win, 'r--', 'LineWidth',1.5);
    xlabel('Time (hours)'); ylabel('Separation Distance (km)');
    title('Focused Separation Distance (±24 hrs of Closest Approach)');
    legend('Before Burn', 'After Burn');
    
    %% Load Earth Topography for 3D Plots
    load('earth_topo.mat','topo');
    [sx, sy, sz] = sphere(100);
    R_earth = 6378.137;
    sx = R_earth*sx; sy = R_earth*sy; sz = R_earth*sz;
    
    %% 3D Orbits Before Burn (Focused)
    r_K6_win = r_K6(idx_win_pre,:);
    r_S6_win = r_S6(idx_win_pre,:);
    
    figure('Name','3D Orbits Before Burn (Focused)','Color','w');
    surf(sx, sy, sz, topo, 'FaceColor','texturemap','EdgeColor','none');
    hold on; axis equal; view(3); grid on;
    plot3(r_K6_win(:,1), r_K6_win(:,2), r_K6_win(:,3), 'r', 'LineWidth',1.5);
    plot3(r_S6_win(:,1), r_S6_win(:,2), r_S6_win(:,3), 'b', 'LineWidth',1.5);
    
    % Collision point
    collision_pt_K6 = r_K6(idx,:);
    plot3(collision_pt_K6(1), collision_pt_K6(2), collision_pt_K6(3),'ko','MarkerFaceColor','k', 'MarkerSize',12, 'DisplayName','Collision Pt');
    
    xlabel('ECI X (km)'); ylabel('ECI Y (km)'); zlabel('ECI Z (km)');
    title('3D Orbits Before Burn (±24 hrs)');
    legend('Earth','K6','S6','Collision Pt');
    
    %% 3D Orbits After Burn (Focused)
    r_K6_post_win = r_K6_post(idx_win_post,:);
    r_S6_post_win = r_S6_post(idx_win_post,:);
    
    figure('Name','3D Orbits After Burn (Focused)','Color','w');
    surf(sx, sy, sz, topo, 'FaceColor','texturemap','EdgeColor','none');
    hold on; axis equal; view(3); grid on;
    plot3(r_K6_post_win(:,1), r_K6_post_win(:,2), r_K6_post_win(:,3), 'r', 'LineWidth',1.5);
    plot3(r_S6_post_win(:,1), r_S6_post_win(:,2), r_S6_post_win(:,3), 'b', 'LineWidth',1.5);
    
    % % Burn point (at t_min)
    % burn_pt_S6 = Y(idx,7:9);
    % plot3(burn_pt_S6(1), burn_pt_S6(2), burn_pt_S6(3), 'g*', 'MarkerSize',10, 'DisplayName','Burn Pt');
    % 
    % % Avoided closest approach (after burn)
    % [~, idx_post_min] = min(sep_post);
    % avoid_pt_S6 = r_S6_post(idx_post_min,:);
    % plot3(avoid_pt_S6(1), avoid_pt_S6(2), avoid_pt_S6(3), 'k*', 'MarkerSize',10, 'DisplayName','Avoided Closest Pt');
    
    xlabel('ECI X (km)'); ylabel('ECI Y (km)'); zlabel('ECI Z (km)');
    title('3D Orbits After Burn (±24 hrs)');
    legend('Earth','K6','S6 Post-Burn','Burn Pt','Avoided Closest Pt');

end

%% Two satellite dynamics with J2
function dydt = two_sat_orbit(~, y, mu)
    dydt = zeros(12,1);
    dydt(1:6) = orbitalDynamicsJ2([], y(1:6), mu, 'yes');
    dydt(7:12) = orbitalDynamicsJ2([], y(7:12), mu, 'yes');
end
