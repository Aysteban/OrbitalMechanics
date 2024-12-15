clear variables

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]

% Initial condition
r0 = [ 26578; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];

Re = 6378.1;
J2 = 0.00108263;

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E ); % Semi-major axis [km]
T1 = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
T2 = T1;
tspan1 = linspace( 0, 2000*T1, 1000 );
tspan2 = linspace( 0, 2000*T2, 1000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ t1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan1, y0, options );

[ t2, Y2 ] = ode113( @(t,y) ode_2bpp(t,y,mu_E, J2, Re), tspan2, y0, options );

% Plot the results
figure();
plot3(Y1(:,1), Y1(:,2), Y1(:,3), 'b-', 'LineWidth', 1); % Courbe bleue
hold on;
plot3(Y2(:,1), Y2(:,2), Y2(:,3), 'r-', 'LineWidth', 1); % Courbe rouge
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;
legend('Orbit 1', 'Orbit 2'); % Ajout d'une légende pour différencier les courbes


e=(1/mu_E)*cross(v0,cross(r0,v0))-r0/norm(r0); % Orbit's eccentricity