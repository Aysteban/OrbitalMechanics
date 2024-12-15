clear variables

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]

% Initial condition
r0 = [ 26578; 0; 0 ]; % [km]
v0 = [ 0; 2.221; 3.173 ]; % [km/s]
y0 = [ r0; v0 ];
t0 = 0;
e= norm((1/mu_E)*cross(v0,cross(r0,v0))-r0/norm(r0)); % Orbit's eccentricity
a = 1/( 2/norm(r0) - dot(v0,v0)/mu_E );% Semi-major axis [km]

% Exercice 3b
a = 7000;
n = sqrt(mu_E/a^3); 

% Set time quantities
T = 2*pi/n; % Orbital period [s]
tf= 2*T;
t = tf - t0;
ka = floor(t/T);
t_a = t - ka*T;
N = 1000;
tspan = linspace(t0, tf, N);

% Solution with fzero (3a)
E0a = n*t + (e*sin(n*t))/(1-sin(n*t+e) + sin(n*t));
E_bar = fzero(@(E) fun(E, a, mu_E, e, t_a), E0a);
Ea = E_bar + 2*ka*pi;

% Solutions with fzero (3b)
Eb = zeros(size(tspan));
Eb2 = zeros(size(tspan));
Eb3 = zeros(size(tspan));
Eb4 = zeros(size(tspan));
Eb5 = zeros(size(tspan));
Eb6 = zeros(size(tspan));

for i = 1:N
    t_b = tspan(i);
    kb = floor(t_b/T);
    t_b = t_b - kb * T;
    e=0;
    E0b = n*t_b + (e*sin(n*t_b))/(1-sin(n*t_b+e) + sin(n*t_b));
    Eb(i) = (fzero(@(E) fun(E, a, mu_E, e, t_b), E0b) + 2*kb*pi) * 180/pi;
    e=0.2;
    E0b = n*t_b + (e*sin(n*t_b))/(1-sin(n*t_b+e) + sin(n*t_b));
    Eb2(i) = (fzero(@(E) fun(E, a, mu_E, e, t_b), E0b) + 2*kb*pi) * 180/pi;
    e=0.4;
    E0b = n*t_b + (e*sin(n*t_b))/(1-sin(n*t_b+e) + sin(n*t_b));
    Eb3(i) = (fzero(@(E) fun(E, a, mu_E, e, t_b), E0b) + 2*kb*pi) * 180/pi;
    e=0.6;
    E0b = n*t_b + (e*sin(n*t_b))/(1-sin(n*t_b+e) + sin(n*t_b));
    Eb4(i) = (fzero(@(E) fun(E, a, mu_E, e, t_b), E0b) + 2*kb*pi) * 180/pi;
    e=0.8;
    E0b = n*t_b + (e*sin(n*t_b))/(1-sin(n*t_b+e) + sin(n*t_b));
    Eb5(i) = (fzero(@(E) fun(E, a, mu_E, e, t_b), E0b) + 2*kb*pi) * 180/pi;
    e=0.95;
    E0b = n*t_b + (e*sin(n*t_b))/(1-sin(n*t_b+e) + sin(n*t_b));
    Eb6(i) = (fzero(@(E) fun(E, a, mu_E, e, t_b), E0b) + 2*kb*pi) * 180/pi;
end

figure (); 
plot(tspan/T, Eb, 'b-');
hold on;
plot(tspan/T, Eb2, 'r-');
plot(tspan/T, Eb3, 'g-');
plot(tspan/T, Eb4, 'b-');
plot(tspan/T, Eb5, 'r-');
plot(tspan/T, Eb6, 'g-');
xlabel('t[T]'); 
ylabel('E (deg)');
title('Anomalie Excentrique E en Fonction du Temps');
grid on; 
legend('e = 0', 'e = 0.2', 'e = 0.4', 'e = 0.6', 'e = 0.8', 'e = 0.95');