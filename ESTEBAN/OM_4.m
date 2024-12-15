clc
clear variables

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(13);

% Initial condition
r1 = [ -21800; 37900; 0 ]; % [km]
r2 = [ 27300; 27700; 0 ]; % [km]

% Set time quantities
t1 = 0;
t2= 54400;
dt = t2 - t1;

N = 50000;
tspan = linspace(t1, 4*t2, N);

% Solver
[a,P,E,ERROR,v1,v2,TPAR,THETA] = lambertMR(r1,r2,dt,mu_E,0,0,0,0);

disp(a);
disp(v1);
disp(v2);

y1 = [ r1; v1' ];

a = 1/( 2/norm(r1) - dot(v1,v1)/mu_E ); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]

% Orbit propagation
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

[ t1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y1, options );

% Plot the results

figure();
plot3(Y1(:,1), Y1(:,2), Y1(:,3), 'b-', 'LineWidth', 1); % Courbe bleue
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Two-body problem orbit');
grid on;