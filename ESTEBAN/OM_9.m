clc
clear variables

% Physical parameters

mu_E = astroConstants(13);
R_E = 6378;
mu_S = astroConstants(4);
AU = astroConstants(2);

% Initial condition

rE = [0; -1; 0]*AU;
rSOI = 145.3*R_E;

% Powered gravity assist solver

Vpl = sqrt(mu_S/norm(rE))*[1; 0; 0];

Vm = [31.5; 5.2; 0];
Vp = [36.0; 0; 0];

vinfm_vec = Vm - Vpl;
vinfp_vec = Vp - Vpl;

[vinfm, vinfp, delta, rp, am, ap, em, ep, vpm, vpp, deltam, deltap, dv, dvp] = hyperbolic_powered(vinfm_vec, vinfp_vec, mu_E);

% Initial conditions planetocentric

u = cross(vinfm_vec,vinfp_vec)/norm(cross(vinfm_vec,vinfp_vec));

betam = pi/2 - deltam/2;

dir_vm = vinfm_vec/norm(vinfm_vec); % Vinf- velocity direction
dir_vp = vinfp_vec/norm(vinfp_vec); % Vinf+ velocity direction

dirm = Rotate(dir_vm, u, deltam/2); 
dirp = Rotate(dir_vp, u, -deltap/2);

r0 = rp * Rotate(dir_vm, u, -betam);

vm = vpm*dirm;
vp = vpp*dirp;

% Time span planetocentric

tspan_m = linspace(0, -50000, 100000);

tspan_p = linspace(0, 50000, 100000);

% Set options for ODE solver

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Integration of planetocentric trajectory

y0m = [r0; vm];
y0p = [r0; vp];

[t1, Y1] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_m, y0m, options);

[t2, Y2] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_p, y0p, options);

% PLANETOCENTRIC PLOT

figure();
hold on;

plot(Y1(:, 1) / R_E, Y1(:, 2) / R_E, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Flyby hyperbola (infront)');
plot(Y2(:, 1) / R_E, Y2(:, 2) / R_E, 'b-', 'LineWidth', 1.5);
plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');

xlabel('x [Re]');
ylabel('y [Re]');
title('Trajectory in Earth-centred frame parallel to (HECI)');
axis equal;
grid on;

xlim([-6, 9]);
ylim([-9, 3]);





