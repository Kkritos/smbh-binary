"""Few-body Newtonian integrator for SMBH + BH + star encounters.
Units: AU, yr, Msun -> G = 4 pi^2.  2D (coplanar).
Bodies order: 0 = SMBH, 1 = BH, 2 = star.
State vector: [x0,y0, x1,y1, x2,y2, vx0,vy0, vx1,vy1, vx2,vy2].
"""
import numpy as np
from scipy.integrate import solve_ivp

G = 4.0 * np.pi**2
RSUN_AU = 695700.0 / 1.495978707e8          # 0.004649 AU

def accel(pos, m, soft=0.0):
    n = len(m); a = np.zeros((n, 2))
    for i in range(n):
        for j in range(n):
            if i == j: continue
            d = pos[j] - pos[i]
            r2 = d @ d + soft*soft
            a[i] += G * m[j] * d / r2**1.5
    return a

def rhs(t, y, m, soft):
    pos = y[:6].reshape(3, 2)
    vel = y[6:].reshape(3, 2)
    a = accel(pos, m, soft)
    return np.concatenate([vel.ravel(), a.ravel()])

def kepler_state(Mcen, rp, e, nu, mu_total=None):
    """Position+velocity of a test body on a Kepler orbit about a mass at origin.
    rp = pericenter (AU), e = eccentricity, nu = true anomaly (rad).
    Returns (pos[2], vel[2]) in the orbital plane (periapsis along +x)."""
    mu = G * Mcen if mu_total is None else G * mu_total
    p = rp * (1.0 + e)                      # semi-latus rectum
    r = p / (1.0 + e*np.cos(nu))
    pos = np.array([r*np.cos(nu), r*np.sin(nu)])
    h = np.sqrt(mu * p)
    vr = (mu / h) * e * np.sin(nu)
    vt = (mu / h) * (1.0 + e*np.cos(nu))
    # radial + transverse -> cartesian
    er = np.array([np.cos(nu), np.sin(nu)])
    et = np.array([-np.sin(nu), np.cos(nu)])
    vel = vr*er + vt*et
    return pos, vel

def make_ic(M_smbh, m_bh, m_star, a_bin, rp_out, e_out, nu0, phase,
            incl_bin=0.0, bin_sense=+1, e_bin=0.0):
    """Build the 12-D initial state.
    - binary CoM placed on a Kepler orbit (rp_out,e_out,nu0) about the SMBH
    - inner binary: eccentricity e_bin, started at inner APOCENTER, orientation `phase`
      (e_bin=0 -> circular, separation a_bin)
    - bin_sense = +/-1 sets prograde/retrograde inner orbit
    """
    m_b = m_bh + m_star
    com_pos, com_vel = kepler_state(M_smbh + m_b, rp_out, e_out, nu0)
    # inner binary (relative to CoM): at apocenter, separation a(1+e), transverse v
    ang = phase + incl_bin
    u = np.array([np.cos(ang), np.sin(ang)])          # apocenter separation direction
    w = bin_sense * np.array([-np.sin(ang), np.cos(ang)])
    sep = a_bin * (1.0 + e_bin)
    v_bin = np.sqrt(G * m_b / a_bin * (1.0 - e_bin) / (1.0 + e_bin))   # apocenter speed
    # star at +(m_bh/m_b)*sep, BH at -(m_star/m_b)*sep
    r_star =  (m_bh / m_b) * sep * u
    r_bh   = -(m_star / m_b) * sep * u
    v_star =  (m_bh / m_b) * v_bin * w
    v_bh   = -(m_star / m_b) * v_bin * w
    pos = np.array([[0.0, 0.0], com_pos + r_bh, com_pos + r_star])
    vel = np.array([[0.0, 0.0], com_vel + v_bh, com_vel + v_star])
    # shift to system CoM rest frame (so SMBH recoil is handled cleanly)
    m = np.array([M_smbh, m_bh, m_star])
    vcom = (m[:, None] * vel).sum(0) / m.sum()
    pcom = (m[:, None] * pos).sum(0) / m.sum()
    pos -= pcom; vel -= vcom
    return np.concatenate([pos.ravel(), vel.ravel()]), m

def integrate(y0, m, t_span, n_frames, soft=0.0, rtol=1e-11, atol=1e-12,
              events=None):
    t_eval = np.linspace(t_span[0], t_span[1], n_frames)
    sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, args=(m, soft),
                    method="DOP853", rtol=rtol, atol=atol,
                    dense_output=True, events=events, max_step=(t_span[1]-t_span[0])/2000)
    return sol

def classify(sol, m):
    """Energies at the final time, each body relative to the SMBH."""
    y = sol.y[:, -1]
    pos = y[:6].reshape(3, 2); vel = y[6:].reshape(3, 2)
    out = {}
    for name, i in [("bh", 1), ("star", 2)]:
        dr = pos[i] - pos[0]; dv = vel[i] - vel[0]
        r = np.linalg.norm(dr); v2 = dv @ dv
        E = 0.5*v2 - G*m[0]/r
        out[name] = dict(E=E, bound=(E < 0), r=r, v=np.sqrt(v2))
    # inner binary binding
    dr = pos[2]-pos[1]; dv = vel[6:].reshape(3,2)[2]-vel[1] if False else (vel[2]-vel[1])
    r = np.linalg.norm(dr); v2 = dv@dv
    mb = m[1]+m[2]
    out["binary"] = dict(E=0.5*v2 - G*mb/r, bound=(0.5*v2 - G*mb/r < 0), sep=r)
    return out

# ----------------------------------------------------------------------------
# osculating elements (2D), Peters GW inspiral, debris  -- for preview rendering
# ----------------------------------------------------------------------------
C_AUYR = 299792.458 / (1.495978707e8 / (365.25*86400) )   # speed of light in AU/yr

def osculating(rvec, vvec, Mtot):
    """2D osculating elements of a body about a mass Mtot at origin."""
    mu = G * Mtot
    r = np.linalg.norm(rvec); v2 = vvec @ vvec
    E = 0.5*v2 - mu/r
    a = -mu/(2*E)
    evec = ((v2 - mu/r)*rvec - (rvec @ vvec)*vvec)/mu
    e = np.linalg.norm(evec)
    omega = np.arctan2(evec[1], evec[0])
    # true anomaly
    nu = np.arctan2(rvec[1], rvec[0]) - omega
    nu = (nu + np.pi) % (2*np.pi) - np.pi
    return dict(a=a, e=e, omega=omega, nu=nu, rp=a*(1-e), E=E)

def peters_inspiral(a0, e0, m1, m2, a_stop, n_steps=4000):
    """Orbit-averaged Peters(1964) a(t), e(t) from (a0,e0) until a < a_stop.
    PREVIEW ONLY -- the production pipeline integrates full 2.5PN in SpaceHub."""
    M = m1 + m2; c = C_AUYR
    K = (G**3) * m1 * m2 * M / c**5
    def dadt(a, e):
        return -(64/5)*K/(a**3*(1-e*e)**3.5)*(1 + 73/24*e*e + 37/96*e**4)
    def dedt(a, e):
        return -(304/15)*e*K/(a**4*(1-e*e)**2.5)*(1 + 121/304*e*e)
    # integrate in a as independent-ish var with small RK4 steps in time
    a, e, t = a0, e0, 0.0
    A, Ecc, T = [a], [e], [t]
    # adaptive-ish dt: fraction of current inspiral time
    while a > a_stop and e >= 0 and len(T) < 200000:
        tinsp = a/abs(dadt(a, e))
        dt = tinsp/300.0
        k1a, k1e = dadt(a, e), dedt(a, e)
        k2a, k2e = dadt(a+0.5*dt*k1a, e+0.5*dt*k1e), dedt(a+0.5*dt*k1a, e+0.5*dt*k1e)
        k3a, k3e = dadt(a+0.5*dt*k2a, e+0.5*dt*k2e), dedt(a+0.5*dt*k2a, e+0.5*dt*k2e)
        k4a, k4e = dadt(a+dt*k3a, e+dt*k3e), dedt(a+dt*k3a, e+dt*k3e)
        a += dt/6*(k1a+2*k2a+2*k3a+k4a)
        e += dt/6*(k1e+2*k2e+2*k3e+k4e); e = max(e, 0.0)
        t += dt
        A.append(a); Ecc.append(e); T.append(t)
    return np.array(T), np.array(A), np.array(Ecc)

def kepler_solve(M_anom, e):
    E = M_anom if e < 0.8 else np.pi
    for _ in range(80):
        dE = (E - e*np.sin(E) - M_anom)/(1 - e*np.cos(E))
        E -= dE
        if abs(dE) < 1e-13: break
    return E

def ellipse_point(a, e, omega, nu, Mtot):
    """Position (2D) on the osculating ellipse about origin."""
    r = a*(1-e*e)/(1+e*np.cos(nu))
    xp = r*np.cos(nu); yp = r*np.sin(nu)
    c_, s_ = np.cos(omega), np.sin(omega)
    return np.array([c_*xp - s_*yp, s_*xp + c_*yp])

def spawn_debris(r_star, v_star, r_bh, v_bh, m_bh, n=160, spread=0.35):
    """Tidal-disruption debris: frozen-in spread in specific orbital ENERGY about the
    BH -> kick along the orbital-velocity direction so the cloud stretches into a
    stream (bound debris fall back, unbound debris escape). 2D, in the input frame."""
    rng = np.random.default_rng(7)
    dr = r_star - r_bh
    r = np.linalg.norm(dr)
    vrel = v_star - v_bh
    vhat = vrel/np.linalg.norm(vrel)                 # tangential (orbital) direction
    v_orb = np.sqrt(G*m_bh/r)
    dv_scale = spread*v_orb
    kick = np.linspace(-1.6, 1.6, n) + rng.normal(0, 0.12, n)   # ordered energy ladder
    pos = np.tile(r_star, (n, 1)) + rng.normal(0, 0.06*r, (n, 2))
    vel = np.tile(v_star, (n, 1)) + kick[:, None]*dv_scale*vhat
    return pos, vel
