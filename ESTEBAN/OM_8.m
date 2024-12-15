clc
clear variables

% Physical parameters

mu_E = astroConstants(13);
R_E = 6378;
mu_S = astroConstants(4);
AU = astroConstants(2);

% Initial condition

vinfmin_vec = [15.1; 0; 0];
rE = [1; 0; 0]*AU;
rSOI = 145.3*R_E;

IP1 = 9200;
IP2 = 10200;
IP3 = 11200;
IP4 = 12200;
IP5 = 13200;

% Hyperbola solving

[vinf1, a1, delta1, e1, rp1, dv1] = hyperbolic(vinfmin_vec, IP1, mu_E);
[vinf2, a2, delta2, e2, rp2, dv2] = hyperbolic(vinfmin_vec, IP2, mu_E);
[vinf3, a3, delta3, e3, rp3, dv3] = hyperbolic(vinfmin_vec, IP3, mu_E);
[vinf4, a4, delta4, e4, rp4, dv4] = hyperbolic(vinfmin_vec, IP4, mu_E);
[vinf5, a5, delta5, e5, rp5, dv5] = hyperbolic(vinfmin_vec, IP5, mu_E);

% Vinf+ computing

IP_infront1 = IP1*[0; 1; 0];
IP_infront2 = IP2*[0; 1; 0];
IP_infront3 = IP3*[0; 1; 0];
IP_infront4 = IP4*[0; 1; 0];
IP_infront5 = IP5*[0; 1; 0];

u_infront1 = cross(IP_infront1,vinfmin_vec)/norm(cross(IP_infront1,vinfmin_vec));
u_infront2 = cross(IP_infront2,vinfmin_vec)/norm(cross(IP_infront2,vinfmin_vec));
u_infront3 = cross(IP_infront3,vinfmin_vec)/norm(cross(IP_infront3,vinfmin_vec));
u_infront4 = cross(IP_infront4,vinfmin_vec)/norm(cross(IP_infront4,vinfmin_vec));
u_infront5 = cross(IP_infront5,vinfmin_vec)/norm(cross(IP_infront5,vinfmin_vec));

vinfplus_infront_vec1 = Rotate(vinfmin_vec, u_infront1, delta1);
vinfplus_infront_vec2 = Rotate(vinfmin_vec, u_infront2, delta2);
vinfplus_infront_vec3 = Rotate(vinfmin_vec, u_infront3, delta3);
vinfplus_infront_vec4 = Rotate(vinfmin_vec, u_infront4, delta4);
vinfplus_infront_vec5 = Rotate(vinfmin_vec, u_infront5, delta5);

% V+ and V- computing

Vp = sqrt(mu_S/norm(rE))*[0; 1; 0];

Vmin = Vp + vinfmin_vec;

Vplus_infront1 = Vp + vinfplus_infront_vec1;
Vplus_infront2 = Vp + vinfplus_infront_vec2;
Vplus_infront3 = Vp + vinfplus_infront_vec3;
Vplus_infront4 = Vp + vinfplus_infront_vec4;
Vplus_infront5 = Vp + vinfplus_infront_vec5;

% Initial parameters

r0 = rE;

a_before = 1 / (2 / norm(r0) - dot(Vmin, Vmin) / mu_S);
T_before = 2 * pi * sqrt(a_before^3 / mu_S);

a_after_infront1 = 1 / (2 / norm(r0) - dot(Vplus_infront1, Vplus_infront1) / mu_S);
T_after_infront1 = 2 * pi * sqrt(a_after_infront1^3 / mu_S);

a_after_infront2 = 1 / (2 / norm(r0) - dot(Vplus_infront2, Vplus_infront2) / mu_S);
T_after_infront2 = 2 * pi * sqrt(a_after_infront2^3 / mu_S);

a_after_infront3 = 1 / (2 / norm(r0) - dot(Vplus_infront3, Vplus_infront3) / mu_S);
T_after_infront3 = 2 * pi * sqrt(a_after_infront3^3 / mu_S);

a_after_infront4 = 1 / (2 / norm(r0) - dot(Vplus_infront4, Vplus_infront4) / mu_S);
T_after_infront4 = 2 * pi * sqrt(a_after_infront4^3 / mu_S);

a_after_infront5 = 1 / (2 / norm(r0) - dot(Vplus_infront5, Vplus_infront5) / mu_S);
T_after_infront5 = 2 * pi * sqrt(a_after_infront5^3 / mu_S);


% Time span 

tspan_before = linspace(0, -T_before/3, 100000); 

tspan_after_infront = linspace(0, 2*T_after_infront1/3, 100000);

% Set options for ODE solver

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

% Integration of heliocentric trajectory

y0_before = [r0; Vmin];
[t1, Y1] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_before, y0_before, options);

y0_after_infront1 = [r0; Vplus_infront1];
[t2a, Y2a] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_infront, y0_after_infront1, options);

y0_after_infront2 = [r0; Vplus_infront2];
[t2b, Y2b] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_infront, y0_after_infront2, options);

y0_after_infront3 = [r0; Vplus_infront3];
[t2c, Y2c] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_infront, y0_after_infront3, options);

y0_after_infront4 = [r0; Vplus_infront4];
[t2d, Y2d] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_infront, y0_after_infront4, options);

y0_after_infront5 = [r0; Vplus_infront5];
[t2e, Y2e] = ode113(@(t, y) ode_2bp(t, y, mu_S), tspan_after_infront, y0_after_infront5, options);

% Plot trajectories

% IN-FRONT HELIOCENTRIC

figure();
hold on;

plot(Y1(:, 1) / AU, Y1(:, 2) / AU, 'c-', 'LineWidth', 2.5, 'DisplayName', 'Before flyby');

plot(Y2a(:, 1) / AU, Y2a(:, 2) / AU, 'b-', 'LineWidth', 2.5, 'DisplayName', 'IP = 9200');
plot(Y2b(:, 1) / AU, Y2b(:, 2) / AU, 'r-', 'LineWidth', 2.5, 'DisplayName', 'IP = 10200');
plot(Y2c(:, 1) / AU, Y2c(:, 2) / AU, 'LineWidth', 2, 'DisplayName', 'IP = 11200');
plot(Y2d(:, 1) / AU, Y2d(:, 2) / AU, 'm-', 'LineWidth', 2.5, 'DisplayName', 'IP = 12200');
plot(Y2e(:, 1) / AU, Y2e(:, 2) / AU, 'g-', 'LineWidth', 2.5, 'DisplayName', 'IP = 13200');

plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'yellow', 'DisplayName', 'Sun');

xlabel('x [AU]');
ylabel('y [AU]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
legend show;
axis equal;
grid on;

xlim([-1.5, 1.5]);
ylim([-1, 1.5]);

% Initial conditions planetocentric

dir_infront1 = Rotate([1; 0; 0], u_infront1, delta1/2);
dir_infront2 = Rotate([1; 0; 0], u_infront2, delta2/2);
dir_infront3 = Rotate([1; 0; 0], u_infront3, delta3/2);
dir_infront4 = Rotate([1; 0; 0], u_infront4, delta4/2);
dir_infront5 = Rotate([1; 0; 0], u_infront5, delta5/2);

beta1 = pi/2 - delta1/2;
beta2 = pi/2 - delta2/2;
beta3 = pi/2 - delta3/2;
beta4 = pi/2 - delta4/2;
beta5 = pi/2 - delta5/2;

d_infront1 = Rotate([1; 0; 0], u_infront1, -beta1);
d_infront2 = Rotate([1; 0; 0], u_infront2, -beta2);
d_infront3 = Rotate([1; 0; 0], u_infront3, -beta3);
d_infront4 = Rotate([1; 0; 0], u_infront4, -beta4);
d_infront5 = Rotate([1; 0; 0], u_infront5, -beta5);

r0p_infront1 = rp1 * d_infront1;
r0p_infront2 = rp2 * d_infront2;
r0p_infront3 = rp3 * d_infront3;
r0p_infront4 = rp4 * d_infront4;
r0p_infront5 = rp5 * d_infront5;

v_infront1 = sqrt(2*mu_E*(1./rp1 + 1/2/abs(a1)))*dir_infront1;
v_infront2 = sqrt(2*mu_E*(1./rp2 + 1/2/abs(a2)))*dir_infront2;
v_infront3 = sqrt(2*mu_E*(1./rp3 + 1/2/abs(a3)))*dir_infront3;
v_infront4 = sqrt(2*mu_E*(1./rp4 + 1/2/abs(a4)))*dir_infront4;
v_infront5 = sqrt(2*mu_E*(1./rp5 + 1/2/abs(a5)))*dir_infront5;

% Time span planetocentric

tspan_before_p = linspace(0, -5000, 100000);
tspan_after_p = linspace(0, 5000, 100000);

% Integration of planetocentric trajectory

y0_infront_p1 = [r0p_infront1; v_infront1];
y0_infront_p2 = [r0p_infront2; v_infront2];
y0_infront_p3 = [r0p_infront3; v_infront3];
y0_infront_p4 = [r0p_infront4; v_infront4];
y0_infront_p5 = [r0p_infront5; v_infront5];

[t1ap, Y1ap] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_infront_p1, options);
[t1bp, Y1bp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_infront_p2, options);
[t1cp, Y1cp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_infront_p3, options);
[t1dp, Y1dp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_infront_p4, options);
[t1ep, Y1ep] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_before_p, y0_infront_p5, options);

[t2ap, Y2ap] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_infront_p1, options);
[t2bp, Y2bp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_infront_p2, options);
[t2cp, Y2cp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_infront_p3, options);
[t2dp, Y2dp] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_infront_p4, options);
[t2ep, Y2ep] = ode113(@(t, y) ode_2bp(t, y, mu_E), tspan_after_p, y0_infront_p5, options);

% IN-FRONT PLANETOCENTRIC

figure();
hold on;

plot(Y1ap(:, 1) / R_E, Y1ap(:, 2) / R_E, 'b-', 'LineWidth', 2);
plot(Y1bp(:, 1) / R_E, Y1bp(:, 2) / R_E, 'r-', 'LineWidth', 2);
plot(Y1cp(:, 1) / R_E, Y1cp(:, 2) / R_E, 'y-', 'LineWidth', 2);
plot(Y1dp(:, 1) / R_E, Y1dp(:, 2) / R_E, 'm-', 'LineWidth', 2);
plot(Y1ep(:, 1) / R_E, Y1ep(:, 2) / R_E, 'g-', 'LineWidth', 2);

plot(Y2ap(:, 1) / R_E, Y2ap(:, 2) / R_E, 'b-', 'LineWidth', 2, 'DisplayName', 'IP = 9200');
plot(Y2bp(:, 1) / R_E, Y2bp(:, 2) / R_E, 'r-', 'LineWidth', 2, 'DisplayName', 'IP = 10200');
plot(Y2cp(:, 1) / R_E, Y2cp(:, 2) / R_E, 'y-', 'LineWidth', 2, 'DisplayName', 'IP = 11200');
plot(Y2dp(:, 1) / R_E, Y2dp(:, 2) / R_E, 'm-', 'LineWidth', 2, 'DisplayName', 'IP = 12200');
plot(Y2ep(:, 1) / R_E, Y2ep(:, 2) / R_E, 'g-', 'LineWidth', 2, 'DisplayName', 'IP = 13200');

plot(0, 0, 'yo', 'MarkerSize', 15, 'MarkerFaceColor', 'blue', 'DisplayName', 'Earth');

xlabel('x [Re]');
ylabel('y [Re]');
title('Trajectory in heliocentric ecliptic (HECI) frame');
axis equal;
grid on;

xlim([-6, 6]);
ylim([-3, 6]);