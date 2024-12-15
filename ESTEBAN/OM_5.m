clc
clear variables

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(13);

% Initial condition
kep1 = [12500; 0; 0; 0; 0; 120];
kep2 = [9500; 0.3; 0; 0; 0; 250];

[r1, v1i] = kep2cart(12500, 0, 0, 0, 0, deg2rad(120), mu_E);
[r2, v2f] = kep2cart(9500, 0.3, 0, 0, 0, deg2rad(250), mu_E);

% Set time quantities
t1 = 0;
t2= 3300;
dt = t2 - t1;

N = 50000;
tspan = linspace(t1, 50000, N);

% Solver
[a,P,E,ERROR,v1t,v2t,TPAR,THETA] = c

disp(a);
disp(v1t);
disp(v2t);

dv1 = v1t' - v1i;
dv2 = v2f - v2t';
dvtot = norm(dv1) + norm(dv2);

disp(norm(dv1));
disp(norm(dv2));
disp(dvtot);

% Orbit propagation
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

y1 = [ r1; v1t' ];
y2 = [ r1; v1i ];
y3 = [ r2; v2f ];

[ t1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y1, options );
[ t2, Y2 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y2, options );
[ t3, Y3 ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y3, options );

% Plot the results

figure();
plot(Y1(:,1), Y1(:,2), 'y-', 'LineWidth', 1); % Courbe jaune
hold on;
plot(Y2(:,1), Y2(:,2), 'b-', 'LineWidth', 1); % Courbe bleue
plot(Y3(:,1), Y3(:,2), 'g-', 'LineWidth', 1); % Courbe verte
xlabel('X [km]');
ylabel('Y [km]');
axis([-2E4 2E4 -2E4 2E4]);
title('Two-body problem orbit');
axis equal;
grid on;

legend('Transfer orbit', 'Initial orbit', 'Final orbit'); % Ajout d'une légende pour différencier les courbes
hold off;