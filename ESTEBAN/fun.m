function F = fun(E, a, mu, e, t)

n = sqrt(mu/a^3); 

F = n*t + e*sin(E) - E;

end