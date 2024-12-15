clear variables

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_E = astroConstants(23);
wE = 7.2921e-5;

earthImage = imread('EarthTexture.jpg');
earthImage = flipud(earthImage);

% Initial condition
r0 = [ 5493.312; -641.510; 4564.578 ]; % [km]
v0 = [ -4.792; -0.795; 5.656 ]; % [km/s]
t0 = 0;

e= norm((1/mu_E)*cross(v0,cross(r0,v0))-r0/norm(r0)); % Orbit's eccentricity
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E );% Semi-major axis [km]
Agw = 0;

% Set time quantities
n = sqrt(mu_E/a^3);
T = 2*pi/n; % Orbital period [s]
tf= 5*T;
N = 50000;
tspan = linspace(t0, tf, N);

% LAT and LON calculus
GT = groundtrack(r0, v0, Agw, tspan, wE, mu_E);
lon = GT(:, 3);
lat = GT(:, 4);

disp('Premières valeurs de Longitude et Latitude (en degrés) :');
disp(['Longitude (deg): ', num2str(lon(1:5)')]);
disp(['Latitude (deg): ', num2str(lat(1:5)')]);

figure(); 
hold on;

image([-180 180], [-90 90], earthImage); % Afficher l'image en prenant en compte les limites
axis([-180 180 -90 90]); % Limites de longitude et latitude
axis image; % Maintenir les proportions de l'image
set(gca, 'YDir', 'normal'); % Définit la direction de l'axe Y à normale

plot(lon, lat, 'Marker', '.', 'LineStyle', 'none'); % Utilise des marqueurs sans lignes
xlabel('Longitude (deg)'); 
ylabel('Latitude (deg)');
title('Ground track');
grid on; 
legend('');

xlim([-180 180]); % Limites de longitude
ylim([-90 90]);   % Limites de latitude

startPoint = 1; % Premier point
endPoint = length(lon); % Dernier point

% Points de départ et d'arrivée
hStart = plot(lon(startPoint), lat(startPoint), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); 
hEnd = plot(lon(endPoint), lat(endPoint), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 

legend([hStart, hEnd], {'Point de Départ', 'Point Arrivée'}, 'Location', 'best'); %

hold off;