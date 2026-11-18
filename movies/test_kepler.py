#!/usr/bin/env python3
"""Validate fallback_stage._kepler_t_to_peri against direct integration.

The fallback stage integrates the debris only through the brief time-varying encounter
and then continues any not-yet-returned parcel analytically along its (now static) Kepler
ellipse.  This checks that the analytic time-to-next-pericentre matches a direct ODE
integration of the same two-body state, across a range of eccentricities and phases.

  python3 test_kepler.py        # prints a table; rel. errors should be ~1e-10
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fallback_stage as fs
from scipy.integrate import solve_ivp

G = 4*np.pi**2
GM = G*1.0                                    # 1 Msun central point mass (AU, yr units)


def numeric_t_to_peri(dr, dv, Pmax):
    def rhs(t, y):
        r = y[:3]; rr = (r @ r)**1.5
        return np.concatenate([y[3:], -GM*r/rr])
    tg = np.linspace(1e-6, 3*Pmax, 400000)
    sol = solve_ivp(rhs, (0, 3*Pmax), np.concatenate([dr, dv]), t_eval=tg,
                    method="DOP853", rtol=1e-11, atol=1e-13)
    r = np.linalg.norm(sol.y[:3], axis=0)
    d = np.diff(r); locmin = np.where((d[:-1] < 0) & (d[1:] >= 0))[0] + 1
    runmax = np.maximum.accumulate(r)
    for i in locmin:
        if runmax[i-1] > r[i]*1.0001:
            y0, y1, y2 = r[i-1], r[i], r[i+1]; den = y0 - 2*y1 + y2
            return tg[i] + (0.5*(y0 - y2)/den)*(tg[i+1]-tg[i]) if den != 0 else tg[i]
    return np.nan


def main():
    rng = np.random.default_rng(3)
    print(f"{'a':>6} {'e':>5} {'M0':>7} | {'analytic':>10} {'numeric':>10} {'rel.err':>9}")
    worst = 0.0
    for _ in range(8):
        a = rng.uniform(0.5, 3.0); e = rng.uniform(0.0, 0.85)
        P = 2*np.pi*np.sqrt(a**3/GM)
        M0 = rng.uniform(0.1, 2*np.pi-0.1)
        E = M0
        for _ in range(80):
            E = M0 + e*np.sin(E)
        r = a*(1 - e*np.cos(E))
        Edot = np.sqrt(GM/a**3)/(1 - e*np.cos(E))
        dr = np.array([a*(np.cos(E)-e), a*np.sqrt(1-e*e)*np.sin(E), 0.0])
        dv = np.array([-a*np.sin(E)*Edot, a*np.sqrt(1-e*e)*np.cos(E)*Edot, 0.0])
        ana = fs._kepler_t_to_peri(dr, dv, GM)
        num = numeric_t_to_peri(dr, dv, P)
        rel = abs(ana-num)/num if num and np.isfinite(num) else np.nan
        worst = max(worst, rel)
        print(f"{a:6.3f} {e:5.3f} {M0:7.3f} | {ana:10.4f} {num:10.4f} {rel:9.2e}")
    print(f"worst relative error: {worst:.2e}  ->  {'PASS' if worst < 1e-6 else 'FAIL'}")


if __name__ == "__main__":
    main()
