function [vinf, a, delta, e, rp, dv] = hyperbolic(vinfmin_vec, IP, mu)

vinf = norm(vinfmin_vec);

a = -mu/vinf^2;

delta = 2*atan(-a/IP);

e = 1/sin(delta/2);

rp = a*(1-e);

dv = 2*vinf*sin(delta/2);

end