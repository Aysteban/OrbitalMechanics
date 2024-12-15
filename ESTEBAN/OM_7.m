clc
clear variables

% Physical parameters

mu_E = astroConstants(13);
R_E = 6378;
mu_S = astroConstants(4);
AU = astroConstants(2);

% Initial condition

vinfmin_vec = [15.1; 0; 0];
IP = 9200;
rE = [1; 0; 0]*AU;
rSOI = 145.3*R_E;

% Hyperbola solving

vinf = norm(vinfmin_vec);

a = -mu_E/vinf^2;
delta = 2*atan(-a/IP);
e = 1/sin(delta/2);
rp = a*(1-e);

dv = 2*vinf*sin(delta/2);

% Vinf+ computing

IP_infront = IP*[0; 1; 0];
IP_behind = IP*[0; -1; 0];  
IP_under = IP*[0; 0; -1];

u_infront = cross(IP_infront,vinfmin_vec)/norm(cross(IP_infront,vinfmin_vec));
u_behind = cross(IP_behind,vinfmin_vec)/norm(cross(IP_behind,vinfmin_vec));
u_under = cross(IP_under,vinfmin_vec)/norm(cross(IP_under,vinfmin_vec));

vinfplus_infront_vec = Rotate(vinfmin_vec, u_infront, delta);
vinfplus_behind_vec = Rotate(vinfmin_vec, u_behind, delta);
vinfplus_under_vec = Rotate(vinfmin_vec, u_under, delta);

% V+ and V- computing

Vp = sqrt(mu_S/norm(rE))*[0; 1; 0];

Vmin = Vp + vinfmin_vec;

Vplus_infront = Vp + vinfplus_infront_vec;
Vplus_behind = Vp + vinfplus_behind_vec;
Vplus_under = Vp + vinfplus_under_vec;

% Initial parameters

r0 = rE;

a_before = 1 / (2 / norm(r0) - dot(Vmin, Vmin) / mu_S);
T_before = 2 * pi * sqrt(a_before^3 / mu_S);

a_after_infront = 1 / (2 / norm(r0) - dot(Vplus_infront, Vplus_infront) / mu_S);
T_after_infront = 2 * pi * sqrt(a_after_infront^3 / mu_S);

a_after_behind = 1 / (2 / norm(r0) - dot(Vplus_behind, Vplus_behind) / mu_S);
T_after_behind = 2 * pi * sqrt(a_after_behind^3 / mu_S);

a_after_under = 1 / (2 / norm(r0) - dot(Vplus_under, Vplus_under) / mu_S);
T_after_under = 2 * pi * sqrt(a_after_under^3 / mu_S);

% Time span 

tspan_before = linspace(0, -T_before/3, 10000); 
tspan_after_infront = linspace(0, 2*T_after_infront/3, 10000); 
tspan_after_behind = linspace(0, T_after_behind/8, 10000); 
tspan_after_under = linspace(0, T_after_under/2, 10000); 

% Set options for ODE solver

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Integration of heliocentric trajectory

y0_before = [r0; Vmin];
[t1, Y1] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_before, y0_before, options);

y0_after_infront = [r0; Vplus_infront];
[t2a, Y2a] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_infront, y0_after_infront, options);

y0_after_behind = [r0; Vplus_behind];
[t2b, Y2b] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_behind, y0_after_behind, options);

y0_after_under = [r0; Vplus_under];
[t2c, Y2c] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_under, y0_after_under, options);

% Plot trajectories

% IN-FRONT HELIOCENTRIC

figure();
hold on;

plot(Y1(:, 1) / AU, Y1(:, 2) / AU, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Before flyby');
plot(Y2a(:, 1) / AU, Y2a(:, 2) / AU, 'r-', 'LineWidth', 1.5, 'DisplayName', 'After flyby (in front)');
plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');

xlabel('x [AU]');
ylabel('y [AU]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-1.5, 1.5]);
ylim([-1, 1.5]);
 
% BEHIND HELIOCENTRIC

figure();
hold on;

plot(Y1(:, 1) / AU, Y1(:, 2) / AU, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Before flyby');
plot(Y2b(:, 1) / AU, Y2b(:, 2) / AU, 'r-', 'LineWidth', 1.5, 'DisplayName', 'After flyby (behind)');
plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');

xlabel('x [AU]');
ylabel('y [AU]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-2.5, 2.5]);
ylim([-1, 3]);

% UNDER HELIOCENTRIC

figure();
hold on;

plot3(Y1(:, 1) / AU, Y1(:, 2) / AU, Y1(:, 3) / AU, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Before flyby');
plot3(Y2c(:, 1) / AU, Y2c(:, 2) / AU, Y2c(:, 3) / AU, 'r-', 'LineWidth', 1.5, 'DisplayName', 'After flyby (under)');
plot3(0, 0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');

xlabel('x [AU]');
ylabel('y [AU]');
zlabel('z [AU]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
zlim([-1.5, 1.5]);

view(3);

% Initial conditions planetocentric

dir_infront = Rotate([1; 0; 0], u_infront, delta/2);
dir_behind = Rotate([1; 0; 0], u_behind, delta/2);
dir_under = Rotate([1; 0; 0], u_under, delta/2);

beta = pi/2 - delta/2;

d_infront = Rotate([1; 0; 0], u_infront, -beta);
d_behind = Rotate([1; 0; 0], u_behind, -beta);
d_under = Rotate([1; 0; 0], u_under, -beta);

r0p_infront = rp * d_infront;
r0p_behind = rp * d_behind;
r0p_under = rp * d_under;

v_infront = sqrt(2*mu_E*(1./rp + 1/2/abs(a)))*dir_infront;
v_behind = sqrt(2*mu_E*(1./rp + 1/2/abs(a)))*dir_behind;
v_under = sqrt(2*mu_E*(1./rp + 1/2/abs(a)))*dir_under;

% Time span planetocentric

tspan_before_p = linspace(0, -5000, 100000);
tspan_after_p = linspace(0, 5000, 100000);

% Integration of heliocentric trajectory

y0_infront_p = [r0p_infront; v_infront];
y0_behind_p = [r0p_behind; v_behind];
y0_under_p = [r0p_under; v_under];

[t1ap, Y1ap] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_infront_p, options);
[t2ap, Y2ap] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_infront_p, options);

[t1bp, Y1bp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_behind_p, options);
[t2bp, Y2bp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_behind_p, options);

[t1cp, Y1cp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_under_p, options);
[t2cp, Y2cp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_under_p, options);

% IN-FRONT PLANETOCENTRIC

figure();
hold on;

plot(Y1ap(:, 1) / R_E, Y1ap(:, 2) / R_E, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Flyby hyperbola (infront)');
plot(Y2ap(:, 1) / R_E, Y2ap(:, 2) / R_E, 'b-', 'LineWidth', 1.5);
plot(0, 0, 'yo', 'MarkerSize', 60, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');

xlabel('x [Re]');
ylabel('y [Re]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-4.5, 4.5]);
ylim([-2.5, 5]);
 
% BEHIND PLANETOCENTRIC

figure();
hold on;

plot(Y1bp(:, 1) / R_E, Y1bp(:, 2) / R_E, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Flyby hyperbola (behind)');
plot(Y2bp(:, 1) / R_E, Y2bp(:, 2) / R_E, 'b-', 'LineWidth', 1.5);
plot(0, 0, 'yo', 'MarkerSize', 60, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');

xlabel('x [Re]');
ylabel('y [Re]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-4.5, 4.5]);
ylim([-3, 5]);

% UNDER PLANETOCENTRIC

figure();
hold on;

plot3(Y1cp(:, 1) / R_E, Y1cp(:, 2) / R_E, Y1cp(:, 3) / R_E, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Flyby hyperbola (under)');
plot3(Y2cp(:, 1) / R_E, Y2cp(:, 2) / R_E, Y2cp(:, 3) / R_E, 'b-', 'LineWidth', 1.5);
plot3(0, 0, 0, 'yo', 'MarkerSize', 60, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');

xlabel('x [Re]');
ylabel('y [Re]');
zlabel('z [Re]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-4.5, 4.5]);
ylim([-3, 5]);
zlim([-4, 4]);

view(3);
