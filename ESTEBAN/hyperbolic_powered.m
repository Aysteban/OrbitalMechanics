function [vinfm, vinfp, delta, rp, am, ap, em, ep, vpm, vpp, deltam, deltap, dv, dvp] = hyperbolic_powered(vinfmin_vec, vinfplus_vec, mu)

vinfm = norm(vinfmin_vec);
vinfp = norm(vinfplus_vec);

delta = acos((dot(vinfmin_vec, vinfplus_vec))/(vinfm*vinfp));

rp_test = 6378.136;
options = optimset('Display','off');
rp = fzero( @(rp) -delta + asin(1/(1 + rp/mu*vinfm^2)) + asin(1/(1 + rp/mu*vinfp^2)), rp_test, options);

hatm = 200;
rp_critical = 6378.136 + hatm;

if rp > rp_critical
    am = -mu/vinfm^2;
    ap = -mu/vinfp^2;
    em = 1+rp*vinfm^2/mu;
    ep = 1+rp*vinfp^2/mu;
    vpm = sqrt(2*mu*(1./rp + 1/2/abs(am)));
    vpp = sqrt(2*mu*(1./rp + 1/2/abs(ap)));
    deltam = 2*asin(1/em);
    deltap = 2*asin(1/ep);
    dv = norm(vinfplus_vec - vinfmin_vec);
    dvp = abs(vpp - vpm);
else
    rp = NaN;
end

end