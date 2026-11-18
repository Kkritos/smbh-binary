#!/usr/bin/env python3
"""Peters (1964) GW inspiral (merger) timescale for a bound two-body orbit.

   python3 peters.py A_AU E [--m1 Msun --m2 Msun]

Prints the merger time in YEARS: the time for gravitational-wave emission to shrink
an orbit of semi-major axis a0 and eccentricity e0 to coalescence, from the exact
Peters (1964) result:

  T(a0,e0) = (12/19) c0^4 / beta
             * INT_0^e0  e^(29/19) (1 + 121/304 e^2)^(1181/2299) / (1-e^2)^(3/2) de
  c0   = a0 (1-e0^2) e0^(-12/19) (1 + 121/304 e0^2)^(-870/2299)
  beta = (64/5) G^3 m1 m2 (m1+m2) / c^5

Units: AU, yr, Msun  (G = 4 pi^2 AU^3 Msun^-1 yr^-2,  c = 63241.077 AU/yr), so T is in yr.
The circular limit (e0 -> 0) is T = a0^4 / (4 beta).  Returns nan for a<=0 or e>=1.
Used for the EMRI scene: the captured BH's (a_new, e_new) about the SMBH -> t_merge.
"""
import argparse, numpy as np

G = 4.0 * np.pi**2        # AU^3 / (Msun yr^2)
C = 63241.077             # AU / yr  (speed of light)


def t_merge(a_au, e, m1, m2):
    try:
        a_au = float(a_au); e = float(e)
    except (TypeError, ValueError):
        return float("nan")
    if not (np.isfinite(a_au) and np.isfinite(e)) or a_au <= 0 or e < 0 or e >= 1.0:
        return float("nan")
    beta = (64.0/5.0) * G**3 * m1*m2*(m1+m2) / C**5       # AU^4 / yr
    if e < 1e-6:
        return a_au**4 / (4.0*beta)                        # circular limit
    c0 = (a_au*(1.0-e*e) * e**(-12.0/19.0)
          * (1.0 + 121.0/304.0*e*e)**(-870.0/2299.0))
    ee = np.linspace(0.0, e, 40001)
    integ = (ee**(29.0/19.0) * (1.0 + 121.0/304.0*ee*ee)**(1181.0/2299.0)
             / (1.0 - ee*ee)**1.5)
    integ[0] = 0.0
    I = float(np.sum(0.5*(integ[1:] + integ[:-1]) * np.diff(ee)))   # trapezoid (np.trapz-free)
    return (12.0/19.0) * c0**4 / beta * I


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a", help="semi-major axis [AU]")
    ap.add_argument("e", help="eccentricity")
    ap.add_argument("--m1", type=float, default=1e6, help="mass 1 [Msun] (SMBH)")
    ap.add_argument("--m2", type=float, default=10.0, help="mass 2 [Msun] (BH)")
    a = ap.parse_args()
    print(f"{t_merge(a.a, a.e, a.m1, a.m2):.6g}")


if __name__ == "__main__":
    main()
