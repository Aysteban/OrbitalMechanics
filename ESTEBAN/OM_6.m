clc
clear variables

% Initial condition
Ddate1 = [2003, 4, 1, 0, 0, 0];
Ddate2 = [2003, 8, 1, 0, 0, 0];
Adate1 = [2003, 9, 1, 0, 0, 0];
Adate2 = [2004, 3, 1, 0, 0, 0];

tD1 = date2mjd2000(Ddate1);
tA1 = date2mjd2000(Adate1);
tD2 = date2mjd2000(Ddate2);
tA2 = date2mjd2000(Adate2);

dt1 = tD2 - tD1 + 1;
dt2 = tA2 - tA1 + 1;

% Compute dv
dv= zeros(dt1, dt2);

for i = 1:dt1
    for j = 1:dt2
        t_depart = tD1 + (i - 1);  
        t_arrivee = tA1 + (j - 1);
        dv(i, j) = vtot(t_depart, t_arrivee); 
    end
end

% Create matrixes to plot
[X, Y] = meshgrid(tD1:tD2, tA1:tA2);

p = 31;

X_unique = [X(1, 1), X(1, 1+p:p:end), X(1, end)];
Y_unique = [Y(1, 1); Y(1+p:p:end, 1); Y(end, 1)];

X1 = arrayfun(@(x) mjd20002date(x), X_unique, 'UniformOutput', false);
Y1 = arrayfun(@(y) mjd20002date(y), Y_unique, 'UniformOutput', false);

X2 = cellfun(@(d) sprintf('%04d-%02d-%02d', d(1), d(2), d(3)), X1, 'UniformOutput', false);
Y2 = cellfun(@(d) sprintf('%04d-%02d-%02d', d(1), d(2), d(3)), Y1, 'UniformOutput', false);

% Porkchop plot
figure;
levels = 5:1:10;
[contourMatrix, h] = contour(X, Y, dv', levels, 'LineWidth', 0.5);

clabel(contourMatrix, h, 'FontSize', 8, 'Color', 'k');

colorbar; 
ylabel(colorbar, '\Delta v [km/s]');

% Add ticks and labels
xlabel('Date de départ');
ylabel('Date arrivée');
title('Porkchop plot');
grid on;

set(gca, 'XTick', X_unique, 'XTickLabel', X2);
set(gca, 'YTick', Y_unique, 'YTickLabel', Y2);

% Find minimum of Δv
dvmin = min(min(dv));
disp(dvmin);

for i = 1:dt1
    for j = 1:dt2
        if dv(i,j) == dvmin
            j0 = i + tD1 -1;
            jf = j + tA1 -1;
        end
    end
end

disp(mjd20002date(j0));
disp(mjd20002date(jf));

% Plot Δvmin
hold on;
plot(j0, jf, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b'); 
hold off;

% Orbit propagation
N = 50000;

t0 = j0*86400;
tf = jf*86400;

tspan = linspace(t0, tf, N);

[kep1,~] = uplanet(j0, 3);
[kep2,~] = uplanet(jf, 4);

mu = astroConstants(4);

[r1, v1i] = kep2cart(kep1(1), kep1(2), kep1(3), kep1(4), kep1(5), kep1(6), mu);
[r2, v2f] = kep2cart(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6),mu);

dt = (tf - t0);
[~, ~, ~, ~,v1t, ~, ~, ~] = lambertMR(r1,r2,dt,mu,0,0,0,0);

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

y1 = [ r1; v1t' ];
y2 = [ r1; v1i ];
y3 = [ r2; v2f ];

[ t1, Y1 ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y1, options );
[ t2, Y2 ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y2, options );
[ t3, Y3 ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y3, options );

n = 1.496e+8;

% Plot the transfer trajectory
figure();
plot(Y1(:,1)/n, Y1(:,2)/n, 'g-', 'LineWidth', 1); % Courbe verte
hold on;
plot(Y2(:,1)/n, Y2(:,2)/n, 'b-', 'LineWidth', 1); % Courbe bleue
plot(Y3(:,1)/n, Y3(:,2)/n, 'r-', 'LineWidth', 1); % Courbe rouge
xlabel('X [AU]');
ylabel('Y [AU]');
title('Two-body problem orbit');
axis([-1.75 1.75 -1.75 1.75]);
axis equal;
grid on;

legend('Transfer orbit', 'Initial orbit', 'Final orbit'); % Ajout d'une légende pour différencier les courbes
hold off;



