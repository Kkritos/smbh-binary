// ============================================================================
//  smbh_binary.cpp  --  binary + SMBH three-body scattering with SpaceHub
// ----------------------------------------------------------------------------
//  A single hierarchical three-body experiment: an inner binary (body1 + body2)
//  encounters a central massive object (body0, the "SMBH").  One run = one set
//  of initial conditions, integrated until a terminal event, an escape, or a
//  time backstop; the outcome is then classified and written as one CSV line.
//
//  Two ways to START a run:
//    (A) ELEMENT mode (default): the inner binary + outer orbit are built from
//        orbital elements (a_bin, e_bin, inclination, beta, v_inf, ...).  This
//        is what the sweep generator drives.
//    (B) RESTART mode (--restart FILE): the three bodies are seeded DIRECTLY
//        from Cartesian phase space (mass, radius, type, x,y,z, vx,vy,vz), e.g.
//        the final-state columns this code exports, so a follow-up run can begin
//        exactly where an earlier one ended.
//
//  UNIT SYSTEM (SpaceHub "unit" namespace, G = 1):
//    length  : AU         (so `x_AU * AU` -> internal length)
//    mass    : Msun       (`m_Msun * Ms`)
//    speed   : the unit in which sqrt(G Msun / AU) = V1 = 29.78 km/s; `v_kms*kms`
//    time    : chosen so the period of a 1 AU / 1 Msun orbit is 2*pi; hence
//              years = internal_time / (2*pi), and t_dyn = 2*pi*sqrt(r^3/M).
//    Radii carry Rs (Rsun).  C is the speed of light in internal units, so a
//    gravitational radius is r_g = M / C^2 (G = 1).
// ============================================================================

#include "spaceHub.hpp"   // path supplied by -I$(SPACEHUB), see Makefile
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>      // restart-mode IC file
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>    // std::max({...})
using namespace hub;
using namespace unit;
using namespace force;
using namespace callback;

// ----------------------------------------------------------------------------
//  Force law: full post-Newtonian (1 + 2 + 2.5) on top of Newtonian gravity.
//  PN2.5 is the leading radiation-reaction term, so GW inspiral/merger is
//  captured in-run.  Compile with -DNEWTONIAN_ONLY to drop all PN terms.
// ----------------------------------------------------------------------------
#ifdef NEWTONIAN_ONLY
using Forces = Interactions<NewtonianGrav>;
#else
using Forces = Interactions<NewtonianGrav, PN1, PN2, PN2p5>;
#endif
using Solver   = methods::DefaultMethod<Forces, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar   = Solver::Scalar;

// ----------------------------------------------------------------------------
//  pair_eps: specific relative orbital energy (per reduced mass, G = 1) of the
//  pair (i,j).  eps < 0 => the two are mutually BOUND; eps > 0 => UNBOUND.
//  Used by the classifier to decide which pairs survived as bound systems.
// ----------------------------------------------------------------------------
template<class P>
double pair_eps(P& p, size_t i, size_t j) {
    auto dr = p.pos(i) - p.pos(j);
    auto dv = p.vel(i) - p.vel(j);
    return 0.5*norm2(dv) - (p.mass(i)+p.mass(j))/norm(dr);   // v^2/2 - G M / r
}

// ----------------------------------------------------------------------------
//  eps_to_binary_com: specific relative energy (per reduced mass, G = 1) of body
//  k about the CENTRE OF MASS of the pair {a,b}.  eps < 0 => k is bound to that
//  pair.  This is the right test for "is the third body bound to the binary?"
//  (used for BOUND_TRIPLE) and for the escape criterion (k truly unbound).
// ----------------------------------------------------------------------------
template<class P>
double eps_to_binary_com(P& p, size_t k, size_t a, size_t b) {
    double ma = p.mass(a), mb = p.mass(b), mab = ma + mb;
    auto rcom = (p.pos(a)*ma + p.pos(b)*mb) / mab;           // COM position of {a,b}
    auto vcom = (p.vel(a)*ma + p.vel(b)*mb) / mab;           // COM velocity of {a,b}
    auto dr = p.pos(k) - rcom;
    auto dv = p.vel(k) - vcom;
    double M = p.mass(k) + mab;                              // total mass seen by k
    return 0.5*norm2(dv) - M/norm(dr);
}

// ----------------------------------------------------------------------------
//  pair_plunges: direct-capture (relativistic plunge) test for pair (i,j).
//  Reconstructs the osculating two-body orbit from (r, v) and asks whether its
//  pericenter rp is inside K_PLUNGE * r_g, i.e. the two would plunge together
//  rather than orbit.  GUARD avoids the expensive test when far apart.
//  ASSUMPTION: a geometric 4 r_g capture threshold stands in for the true
//  relativistic separatrix -- adequate for classification, not for waveform.
// ----------------------------------------------------------------------------
template<class P>
bool pair_plunges(P& p, size_t i, size_t j) {
    using namespace hub::consts;
    const double K_PLUNGE = 4.0;      // plunge if pericenter < 4 r_g
    const double GUARD    = 1000.0;   // skip test beyond 1000 r_g (cheap early-out)
    auto dr = p.pos(i) - p.pos(j);
    auto dv = p.vel(i) - p.vel(j);
    double r  = norm(dr);
    double v2 = norm2(dv);
    double M  = p.mass(i) + p.mass(j);
    double Rg = M / (C*C);            // gravitational radius (G = 1)
    if (r > GUARD * Rg) return false;
    double L2 = norm2(cross(dr, dv)); // specific angular momentum squared
    double E  = 0.5*v2 - M/r;         // specific orbital energy
    double rp;
    // pericenter of the conic: rp = (-M + sqrt(M^2 + 2 E L^2)) / (2 E), with the
    // parabolic (E ~ 0) limit rp = L^2 / (2 M) handled separately.
    if (std::abs(E) < 1e-30) {
        rp = L2 / (2.0*M);
    } else {
        rp = (-M + std::sqrt(std::max(0.0, M*M + 2.0*E*L2))) / (2.0*E);
        if (rp < 0) rp = (-M - std::sqrt(std::max(0.0, M*M + 2.0*E*L2))) / (2.0*E);
    }
    return rp < K_PLUNGE * Rg;
}

// ----------------------------------------------------------------------------
//  pair_event: unified terminal-event test for pair (i,j).  Returns
//     0 = nothing, 1 = TDE (sets victim), 2 = PLUNGE (sets victim), 3 = COLLISION
//  `ext[i]` flags an EXTENDED (tidally disruptable) body, i.e. not a BH.
//
//  Decision logic per pair (closest approach is what matters):
//    * r_surf  = R_i + R_j                       (physical contact)
//    * r_tide  = R_v (M_pair / m_v)^(1/3)         (tidal radius of the extended
//                                                  victim v; the easier-to-disrupt
//                                                  of the two if both extended)
//    * r_plunge= 4 r_g
//    - If the tidal radius sticks OUT beyond both the surface and the plunge
//      radius (r_tide > r_surf and > r_plunge) the star is disrupted first: TDE
//      when sep < r_tide.  This automatically switches OFF above the Hills mass:
//      for a heavy enough companion r_tide shrinks inside r_plunge, so the star
//      is swallowed whole (plunge) instead of disrupted.
//    - Otherwise, if the surface is reached before plunge, a COLLISION fires at
//      contact; if the plunge radius is the outer scale, a relativistic PLUNGE
//      fires (victim = the lighter body).
//  ASSUMPTION: disruption is treated as instantaneous at r_tide; no mass loss,
//  partial disruption, or GW-inspiral-to-contact chain is modelled here.
// ----------------------------------------------------------------------------
template<class P>
int pair_event(P& p, size_t i, size_t j, const bool* ext, size_t& victim) {
    using namespace hub::consts;
    double sep = distance(p.pos(i), p.pos(j));
    double mi = p.mass(i), mj = p.mass(j);
    double Ri = p.radius(i), Rj = p.radius(j);
    bool ei = ext[i], ej = ext[j];

    double r_surf = Ri + Rj;

    // tidal radius of each extended body; -1 means "not disruptable" (a BH)
    double rt_i = ei ? Ri * std::pow((mi + mj)/mi, 1.0/3.0) : -1.0;
    double rt_j = ej ? Rj * std::pow((mi + mj)/mj, 1.0/3.0) : -1.0;
    // if both are extended, the one with the LARGER tidal radius disrupts first
    size_t vic; double r_tide;
    if (rt_i >= rt_j) { vic = i; r_tide = rt_i; } else { vic = j; r_tide = rt_j; }

    double Rg = (mi + mj) / (C*C);
    double r_plunge = 4.0 * Rg;
    bool tde_possible = (r_tide > r_surf) && (r_tide > r_plunge);

    if (tde_possible) {                       // disruption is the outermost scale
        if (sep < r_tide) { victim = vic; return 1; }
        return 0;
    }
    if (r_surf >= r_plunge) {                 // two extended/large bodies: contact first
        if (sep < r_surf) return 3;
        return 0;
    }
    // compact pair: test for a relativistic plunge inside 4 r_g
    const double PLUNGE_CHECK = 50.0;         // only bother within 50 r_g
    if (sep < PLUNGE_CHECK * Rg && pair_plunges(p, i, j)) {
        victim = (mi <= mj ? i : j); return 2;   // lighter body is the "victim"
    }
    return 0;
}

int main(int argc, char** argv) {
    // ======================= defaults (ELEMENT mode) =======================
    double M_smbh_in = 1e6;                 // central mass [Msun]
    double m1_in = 10, m2_in = 1;           // inner binary component masses [Msun]
    std::string type1 = "bh", type2 = "star";
    double R1_in = -1, R2_in = -1;          // radii [Rsun]; -1 => use a mass-radius relation
    double a_bin_in = 1, e_bin_in = 0.7;    // inner semimajor axis [AU] / eccentricity
    double incl_in = 0, Omega_in = 0, omega_in = 0, nu_in = 0;   // inner orientation/phase [rad]
    std::string mode = "full";              // full = write trajectories; classify = CSV only
    double rtol_in = 1e-13;                 // integrator relative tolerance
    bool report = true;                     // human-readable summary to stdout

    // ---- outer-orbit controls (ELEMENT mode) ----
    std::string outer = "hyperbolic";       // hyperbolic | elliptic
    double v_inf_in = 100, beta_in = 1;     // hyperbolic: v at infinity [km/s]; beta = r_t_bin/rp
    double n_dyn_in = 1e4;                   // backstop = n_dyn * t_dyn(r_start)
    double a_out_in = 100.0;                 // elliptic: outer semimajor axis [AU]
    double e_out_in = 0.5;                   // elliptic: outer eccentricity
    double nu_out_in = M_PI;                 // elliptic: outer true anomaly at start (pi = apo)
    double n_orb_in = 1000.0;                // elliptic: backstop = n_orb * P_out

    // ---- RESTART-mode controls ----
    bool        restart     = false;         // set by --restart FILE
    std::string restart_file;
    double      t_end_yr_in = -1.0;          // optional integration-time cap [yr]; -1 => auto
    double      r_esc_AU_in = -1.0;          // optional escape radius [AU];      -1 => auto
    const double R_ESC_FAC  = 10.0;          // auto escape radius = 10 x initial system size

    // ============================ parse CLI ============================
    auto usage = [&](){
        std::cerr << "usage: " << argv[0] << " [--flag value ...]\n"
          << " ELEMENT mode (default):\n"
          << "  --M_smbh [Msun]   --m1 [Msun]   --m2 [Msun]\n"
          << "  --type1 [star|wd|ns|bh]  --type2 [...]   --R1 [Rsun,-1=auto]  --R2 [Rsun,-1=auto]\n"
          << "  --a_bin [AU]  --e_bin   --incl [rad]  --Omega [rad]  --omega [rad]  --nu [rad]\n"
          << "  --outer [hyperbolic|elliptic]\n"
          << "    hyperbolic: --beta [rt/rp]  --v_inf [km/s]  --n_dyn [t_dyn units, def 1e4]\n"
          << "    elliptic:   --a_out [AU]    --e_out         --nu_out [rad]  --n_orb [def 1000]\n"
          << " RESTART mode:\n"
          << "  --restart FILE     3 body lines: 'mass_Msun radius_Rsun type x_AU y_AU z_AU vx_kms vy_kms vz_kms'\n"
          << "                     order = body0(SMBH) body1 body2; '#' comments allowed\n"
          << "  --t_end_yr [yr]    integration-time cap (default: n_dyn * t_dyn of initial size)\n"
          << "  --r_esc_AU [AU]    escape radius      (default: 10 x initial system size)\n"
          << " common:\n"
          << "  --mode [full|classify]   --rtol   --report [0|1]\n";
    };
    // --header prints the CSV column names and exits (used to seed the output file)
    for (int a = 1; a < argc; ++a) if (std::string(argv[a]) == "--header") {
        std::cout << "CSV,M_smbh,m1,type1,R1,m2,type2,R2,a_bin,e_bin,outer,beta,v_inf,a_out,e_out,nu_out,"
                     "incl,Omega,omega,nu,rtol,"
                     "outcome,exit_mode,a_new_AU,e_new,v_esc_kms,beta_actual,"
                     "rp_body1_Rg,rp_body2_Rg,rp_inner_AU,R1_used_Rsun,R2_used_Rsun,Rsmbh_Rsun,t_sim_yr,wall_s,"
                     // final phase space of all 3 bodies at termination (system COM frame; AU, km/s)
                     "x_smbh_AU,y_smbh_AU,z_smbh_AU,vx_smbh_kms,vy_smbh_kms,vz_smbh_kms,"
                     "x_b1_AU,y_b1_AU,z_b1_AU,vx_b1_kms,vy_b1_kms,vz_b1_kms,"
                     "x_b2_AU,y_b2_AU,z_b2_AU,vx_b2_kms,vy_b2_kms,vz_b2_kms\n";
        return 0;
    }
    for (int a = 1; a < argc; ++a) {
        std::string key = argv[a];
        if (key == "--help" || key == "-h") { usage(); return 0; }
        if (a + 1 >= argc) { std::cerr << "missing value for " << key << "\n"; usage(); return 1; }
        std::string val = argv[++a];
        if      (key == "--M_smbh")  M_smbh_in = std::stod(val);
        else if (key == "--m1")      m1_in     = std::stod(val);
        else if (key == "--m2")      m2_in     = std::stod(val);
        else if (key == "--type1")   type1     = val;
        else if (key == "--type2")   type2     = val;
        else if (key == "--R1")      R1_in     = std::stod(val);
        else if (key == "--R2")      R2_in     = std::stod(val);
        else if (key == "--a_bin")   a_bin_in  = std::stod(val);
        else if (key == "--e_bin")   e_bin_in  = std::stod(val);
        else if (key == "--incl")    incl_in   = std::stod(val);
        else if (key == "--Omega")   Omega_in  = std::stod(val);
        else if (key == "--omega")   omega_in  = std::stod(val);
        else if (key == "--nu")      nu_in     = std::stod(val);
        else if (key == "--outer")   outer     = val;
        else if (key == "--v_inf")   v_inf_in  = std::stod(val);
        else if (key == "--beta")    beta_in   = std::stod(val);
        else if (key == "--n_dyn")   n_dyn_in  = std::stod(val);
        else if (key == "--a_out")   a_out_in  = std::stod(val);
        else if (key == "--e_out")   e_out_in  = std::stod(val);
        else if (key == "--nu_out")  nu_out_in = std::stod(val);
        else if (key == "--n_orb")   n_orb_in  = std::stod(val);
        else if (key == "--restart") { restart = true; restart_file = val; }
        else if (key == "--t_end_yr")t_end_yr_in = std::stod(val);
        else if (key == "--r_esc_AU")r_esc_AU_in = std::stod(val);
        else if (key == "--mode")    mode      = val;
        else if (key == "--rtol")    rtol_in   = std::stod(val);
        else if (key == "--report")  report    = (val=="1"||val=="true"||val=="on");
        else { std::cerr << "unknown flag: " << key << "\n"; usage(); return 1; }
    }
    bool classify = (mode == "classify");

    auto valid = [](const std::string& t){ return t=="star"||t=="wd"||t=="ns"||t=="bh"; };
    if (mode != "full" && mode != "classify") { std::cerr << "mode must be full|classify\n"; return 1; }
    if (outer != "hyperbolic" && outer != "elliptic") { std::cerr << "outer must be hyperbolic|elliptic\n"; return 1; }
    bool is_ell = (outer == "elliptic");     // only meaningful in ELEMENT mode

    using namespace hub::consts;             // G, C, pi (internal units)

    // These are filled by whichever setup branch runs, then consumed by the
    // shared solver / classifier / output code below.
    Scalar M_smbh = 0, m1 = 0, m2 = 0, R1 = 0, R2 = 0;
    Scalar m_bin = 0, M_tot = 0, r_t_bin = 0, Rg_smbh = 0;
    Scalar r_start = 0, t_backstop = 0, r_far = 0;
    Scalar rp_out = 0, b_out = 0, P_out = 0;      // reporting only
    bool   ext[3] = { false, false, false };      // extended (disruptable) flags
    Particle smbh{0.0, 0.0}, b1{0.0, 0.0}, b2{0.0, 0.0};
    const double NA = std::nan("");

    // small helpers shared by both branches
    auto star_type = [](const std::string& t) {
        if (t == "star") return stellar::StarType::STAR;
        if (t == "wd")   return stellar::StarType::WD;
        if (t == "ns")   return stellar::StarType::NS;
        return stellar::StarType::BH;
    };
    auto radius = [&](const std::string& type, double R_in, Scalar m) -> Scalar {
        if (R_in < 0) return stellar::stellar_radius(star_type(type), m);  // mass-radius relation
        return R_in * Rs;
    };

    if (restart) {
        // ============================ RESTART setup ============================
        // Seed the three bodies straight from Cartesian phase space.  Same physics
        // and output as ELEMENT mode; only the way the ICs are produced differs.
        if (!valid(type1)) type1 = "bh";     // overwritten per-line below
        std::ifstream fin(restart_file);
        if (!fin) { std::cerr << "cannot open restart file: " << restart_file << "\n"; return 1; }

        double       mm[3], RR[3], xx[3], yy[3], zz[3], ux[3], uy[3], uz[3];
        std::string  tt[3];
        int nb = 0; std::string line;
        // file format per body: mass[Msun] radius[Rsun] type x y z[AU] vx vy vz[km/s]
        // order MUST be body0 = SMBH (central), body1, body2 (same index convention
        // the classifier assumes: index 0 is the central mass).
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (nb >= 3) break;
            std::istringstream iss(line);
            if (!(iss >> mm[nb] >> RR[nb] >> tt[nb] >> xx[nb] >> yy[nb] >> zz[nb]
                      >> ux[nb] >> uy[nb] >> uz[nb])) {
                std::cerr << "bad restart line: " << line << "\n"; return 1;
            }
            ++nb;
        }
        if (nb != 3) { std::cerr << "restart file needs exactly 3 body lines (SMBH, b1, b2)\n"; return 1; }

        // echo variables for the CSV (element-only quantities are N/A here)
        M_smbh_in = mm[0]; m1_in = mm[1]; m2_in = mm[2];
        type1 = tt[1]; type2 = tt[2];
        R1_in = RR[1]; R2_in = RR[2];
        a_bin_in = e_bin_in = beta_in = v_inf_in = NA;      // undefined on a restart
        a_out_in = e_out_in = nu_out_in = NA;
        incl_in = Omega_in = omega_in = nu_in = NA;
        outer = "restart";                                  // marks the row's provenance

        // internal-unit masses / radii, and the tidal (extended) flags from type
        M_smbh = mm[0]*Ms; m1 = mm[1]*Ms; m2 = mm[2]*Ms;
        R1 = RR[1]*Rs;     R2 = RR[2]*Rs;
        Scalar R0 = RR[0]*Rs;                               // central-object radius from file
        ext[0] = (tt[0] != "bh"); ext[1] = (tt[1] != "bh"); ext[2] = (tt[2] != "bh");

        // build the particles directly: set mass, radius, then position & velocity
        // (converting AU -> internal length and km/s -> internal speed).
        smbh = Particle{M_smbh, R0}; b1 = Particle{m1, R1}; b2 = Particle{m2, R2};
        smbh.pos.x = xx[0]*AU; smbh.pos.y = yy[0]*AU; smbh.pos.z = zz[0]*AU;
        smbh.vel.x = ux[0]*kms; smbh.vel.y = uy[0]*kms; smbh.vel.z = uz[0]*kms;
        b1.pos.x = xx[1]*AU; b1.pos.y = yy[1]*AU; b1.pos.z = zz[1]*AU;
        b1.vel.x = ux[1]*kms; b1.vel.y = uy[1]*kms; b1.vel.z = uz[1]*kms;
        b2.pos.x = xx[2]*AU; b2.pos.y = yy[2]*AU; b2.pos.z = zz[2]*AU;
        b2.vel.x = ux[2]*kms; b2.vel.y = uy[2]*kms; b2.vel.z = uz[2]*kms;
        // NOTE: we deliberately do NOT re-centre to the COM frame -- the restart
        // reproduces the stored state verbatim (the export was already COM-frame).

        m_bin   = m1 + m2;
        M_tot   = M_smbh + m_bin;
        Rg_smbh = M_tot / (C*C);
        r_t_bin = NA;                        // no outer "binary tidal radius" concept here

        // characteristic system size now = max distance of any body from the COM,
        // used to set the escape radius and (if not overridden) the time backstop.
        Scalar mt = M_smbh + m1 + m2;
        auto Rc = (smbh.pos*M_smbh + b1.pos*m1 + b2.pos*m2) / mt;
        Scalar r_char = std::max({ norm(smbh.pos - Rc), norm(b1.pos - Rc), norm(b2.pos - Rc) });
        if (r_char <= 0) r_char = AU;        // degenerate guard

        r_start = r_char;                                   // reporting only
        r_far   = (r_esc_AU_in > 0) ? r_esc_AU_in*AU : R_ESC_FAC * r_char;
        Scalar t_dyn = 2.0 * pi * std::sqrt(r_char*r_char*r_char / M_tot);
        // backstop: explicit cap if given (yr -> internal via *2pi), else n_dyn * t_dyn
        t_backstop = (t_end_yr_in > 0) ? (t_end_yr_in * 2.0 * pi) : (n_dyn_in * t_dyn);

        if (report) std::cout << "RESTART from " << restart_file
                  << "  r_char=" << r_char/AU << " AU  r_esc=" << r_far/AU
                  << " AU  backstop=" << t_backstop/(2*pi) << " yr\n";

    } else {
        // ============================ ELEMENT setup ============================
        if (!valid(type1) || !valid(type2)) { std::cerr << "type must be star|wd|ns|bh\n"; return 1; }
        if (is_ell && !(e_out_in >= 0.0 && e_out_in < 1.0)) { std::cerr << "elliptic requires 0 <= e_out < 1\n"; return 1; }

        // convert CLI values to internal units
        M_smbh = M_smbh_in * Ms;  m1 = m1_in * Ms;  m2 = m2_in * Ms;
        Scalar a_bin  = a_bin_in  * AU;
        Scalar e_bin  = e_bin_in;
        Scalar incl   = incl_in, Omega = Omega_in, omega = omega_in, nu = nu_in;
        Scalar v_inf  = v_inf_in  * kms;
        Scalar beta   = beta_in;
        Scalar a_out  = a_out_in  * AU;
        Scalar e_out  = e_out_in;
        Scalar nu_out = nu_out_in;

        R1 = radius(type1, R1_in, m1);
        R2 = radius(type2, R2_in, m2);
        // only non-BH bodies can be tidally disrupted; the SMBH (index 0) never is
        ext[0] = false; ext[1] = (type1 != "bh"); ext[2] = (type2 != "bh");

        // derived scales
        m_bin   = m1 + m2;
        M_tot   = M_smbh + m_bin;
        r_t_bin = a_bin * std::pow(M_smbh / m_bin, 1.0/3.0);  // binary's tidal radius about SMBH
        Rg_smbh = M_tot / (C*C);

        if (!is_ell) {
            // HYPERBOLIC outer orbit: pericenter from penetration beta = r_t_bin/rp,
            // impact parameter b from (rp, v_inf) via energy+angular-momentum.
            rp_out  = r_t_bin / beta;
            b_out   = std::sqrt(rp_out*rp_out + 2.0*M_tot*rp_out/(v_inf*v_inf));
            // Start a few tens of binary-tidal-radii out: far enough that the SMBH
            // tide on the inner binary is negligible (~(r_t_bin/r_start)^3 ~ 1e-5 at
            // 50 r_t_bin), but NOT ~5 b_out -- for slow, strongly-focused encounters
            // b_out >> r_t_bin, and starting there forces thousands of wasted inner
            // orbits during a long drift where nothing happens.  The encounter is
            // fully set by (rp_out, b_out), which are unchanged; the inner phase is
            // randomized, so starting closer introduces no bias, only huge speedup.
            r_start = std::max<Scalar>(50.0*r_t_bin, 10.0*rp_out);
            Scalar t_dyn = 2.0 * pi * std::sqrt(r_start*r_start*r_start / M_tot);
            t_backstop = n_dyn_in * t_dyn;    // time cap for captured/undecided cases
            r_far      = 2.0 * r_start;       // escape distance
        } else {
            // ELLIPTIC outer orbit: start at outer true anomaly nu_out (pi = apocenter)
            rp_out  = a_out * (1.0 - e_out);
            r_start = a_out * (1.0 + e_out);
            P_out   = 2.0 * pi * std::sqrt(a_out*a_out*a_out / M_tot);
            t_backstop = n_orb_in * P_out;
            r_far      = 2.0 * r_start;       // (escape stop is disabled for elliptic below)
        }

        if (report) {
            std::cout << "outer=" << outer
                      << "  r_t_bin=" << r_t_bin/AU << " AU  rp_out=" << rp_out/AU << " AU";
            if (is_ell) std::cout << "  a_out=" << a_out/AU << " AU  e_out=" << e_out
                                  << "  P_out=" << P_out/(2*pi) << " yr  backstop=" << t_backstop/(2*pi) << " yr\n";
            else        std::cout << "  b_out=" << b_out/AU << " AU  r_start=" << r_start/AU
                                  << " AU  backstop=" << t_backstop/(2*pi) << " yr\n";
        }

        // ---------- hierarchical initial conditions ----------
        smbh = Particle{M_smbh, 2.0*M_smbh/(C*C)};   // SMBH radius = 2 r_g (Schwarzschild)
        b1   = Particle{m1, R1};
        b2   = Particle{m2, R2};

        // inner binary oriented by the CLI angles (incl = mutual inclination vs outer plane)
        auto inner = orbit::Elliptic(m1, m2, a_bin, e_bin, incl, Omega, omega, nu);
        orbit::move_particles(inner, b2);
        orbit::move_to_COM_frame(b1, b2);

        // outer orbit placed in the reference plane (its own angles fixed to 0)
        if (!is_ell) {
            auto outer_orb = orbit::Hyperbolic(M_smbh, m_bin, v_inf, b_out,
                                               0.0, 0.0, 0.0, r_start, orbit::Hyper::in);
            orbit::move_particles(outer_orb, b1, b2);
        } else {
            auto outer_orb = orbit::Elliptic(M_smbh, m_bin, a_out, e_out,
                                             0.0, 0.0, 0.0, nu_out);
            orbit::move_particles(outer_orb, b1, b2);
        }
        orbit::move_to_COM_frame(smbh, b1, b2);   // integrate in the system COM frame
    }

    // ======================== solver + stop conditions ========================
    Solver solver{0, smbh, b1, b2};               // three-body solver, t0 = 0
    Solver::RunArgs args;
    args.rtol = rtol_in;

    std::string term_event = "none";              // set by the terminal-event stop
    std::string exit_mode  = "backstop";          // how the run ended; see below
    size_t term_victim = 99, term_i = 99, term_j = 99;

    // Track the closest approach of each pair over the whole run (for reporting the
    // achieved pericenters / penetration, independent of where the run stopped).
    double min_sep01 = 1e300, min_sep02 = 1e300, min_sep12 = 1e300;
    args.add_operation([&min_sep01,&min_sep02,&min_sep12](auto& p, auto h) {
        double s01 = distance(p.pos(0), p.pos(1));
        double s02 = distance(p.pos(0), p.pos(2));
        double s12 = distance(p.pos(1), p.pos(2));
        if (s01 < min_sep01) min_sep01 = s01;
        if (s02 < min_sep02) min_sep02 = s02;
        if (s12 < min_sep12) min_sep12 = s12;
    });

    // (1) TIME BACKSTOP: hard cap on integrated time.  A run that neither resolves
    //     nor escapes within this budget stops here with exit_mode = "backstop"
    //     (these become BOUND_TRIPLE after classification).
    args.add_stop_condition((Scalar)t_backstop);

    // (2) TERMINAL EVENT: scan all pairs every step for TDE / PLUNGE / COLLISION.
    args.add_stop_condition([&ext,&term_event,&term_victim,&term_i,&term_j,&exit_mode](auto& p, auto h) -> bool {
        for (size_t i = 0; i < p.number(); ++i)
          for (size_t j = i+1; j < p.number(); ++j) {
            size_t v = 99;
            int e = pair_event(p, i, j, ext, v);
            if (e) {
                term_i = i; term_j = j; term_victim = v;
                term_event = (e==1 ? "TDE" : (e==2 ? "PLUNGE" : "COLLISION"));
                exit_mode  = "terminal";
                return true;
            }
          }
        return false;
    });

    // (3) ESCAPE (enabled except for a bound ELLIPTIC outer orbit, which never
    //     escapes).  A body is a TRUE escaper only if it is (i) beyond r_far from
    //     the COM of the other two, (ii) receding, AND (iii) energetically UNBOUND
    //     from them (eps_to_binary_com > 0).  Clause (iii) stops a body caught at
    //     the apocenter of a wide but BOUND orbit from ending the run prematurely:
    //     such a body instead swings back for more pericenter passages (possibly
    //     resolving to plunge/exchange/collision) or runs on to the backstop.
    if (!is_ell) {
        Scalar r_far_local = r_far;
        args.add_stop_condition([r_far_local,&exit_mode](auto& p, auto h) -> bool {
            size_t n = p.number();                 // == 3
            for (size_t k = 0; k < n; ++k) {
                size_t a = (k == 0 ? 1 : 0);        // the other two bodies
                size_t b = (k == 2 ? 1 : 2);
                double mo = 0;
                typename std::remove_reference_t<decltype(p)>::Vector rc{0,0,0}, vc{0,0,0};
                for (size_t j = 0; j < n; ++j) if (j != k) {
                    mo += p.mass(j); rc += p.pos(j)*p.mass(j); vc += p.vel(j)*p.mass(j);
                }
                rc /= mo; vc /= mo;
                auto dr = p.pos(k) - rc;
                auto dv = p.vel(k) - vc;
                if (norm(dr) > r_far_local && dot(dr, dv) > 0
                    && eps_to_binary_com(p, k, a, b) > 0.0) {
                    exit_mode = "escape"; return true;
                }
            }
            return false;
        });
    }

    // ---------- optional trajectory output (full mode only) ----------
    std::ostringstream tag;
    tag << "bss_" << outer << "_" << type1 << m1_in << "_" << type2 << m2_in
        << "_i" << incl_in << "_" << mode;
    if (!classify) {
        args.add_operation(DefaultWriter(tag.str() + ".txt"));
        args.add_operation(EnergyErrWriter(tag.str() + ".err"));
    }

    // ---------- run (wall-clock timed) ----------
    auto t0 = std::chrono::steady_clock::now();
    solver.run(args);
    auto t1 = std::chrono::steady_clock::now();
    double wall = std::chrono::duration<double>(t1 - t0).count();

    // ============================ classify outcome ============================
    auto& p = solver.particles();
    auto name = [&](size_t id) -> std::string {
        if (id == 0) return "SMBH";
        if (id == 1) return type1 + "1";
        return type2 + "2";
    };

    std::string outcome;
    if (term_event == "TDE")
        outcome = "TDE:" + name(term_victim) + "_by_" + name(term_i==term_victim ? term_j : term_i);
    else if (term_event == "PLUNGE")
        outcome = "PLUNGE:" + name(term_victim) + "_into_" + name(term_i==term_victim ? term_j : term_i);
    else if (term_event == "COLLISION")
        outcome = "COLLISION:" + name(term_i) + "+" + name(term_j);
    else {
        // No terminal event: decide from the final energies which bodies are bound.
        //   inner bound + bound to SMBH  -> BOUND_TRIPLE (hierarchical, unresolved)
        //   inner bound + free of SMBH   -> FLYBY (binary survives, flies off)
        //   inner unbound, one body kept -> EXCHANGE (SMBH captured that body)
        //   inner unbound, both kept     -> BOUND_TRIPLE (both separately bound)
        //   inner unbound, none kept     -> IONIZATION (all three unbound)
        bool binary_bound = pair_eps(p,1,2) < 0;
        if (binary_bound) {
            bool smbh_binary = eps_to_binary_com(p, 0, 1, 2) < 0;
            outcome = smbh_binary ? "BOUND_TRIPLE" : "FLYBY";
        } else {
            bool b01 = pair_eps(p,0,1) < 0;
            bool b02 = pair_eps(p,0,2) < 0;
            if      (b01 && b02) outcome = "BOUND_TRIPLE";
            else if (b01)        outcome = "EXCHANGE_" + type1 + "1";
            else if (b02)        outcome = "EXCHANGE_" + type2 + "2";
            else                 outcome = "IONIZATION";
        }
    }

    // ============================ derived quantities ============================
    auto ae = [&](size_t i, size_t j) {   // osculating (a,e) of a pair, phase-independent
        double u = p.mass(i) + p.mass(j);
        auto [a, e] = orbit::calc_a_e(u, p.pos(i)-p.pos(j), p.vel(i)-p.vel(j));
        return std::make_pair((double)a, (double)e);
    };
    auto vinf_kms = [&](size_t k) { return norm(p.vel(k)) / kms; };

    double a_new = NA, e_new = NA, v_esc = NA, beta_actual = NA;
    double rp_bh1 = min_sep01/Rg_smbh;            // body1 closest approach to SMBH, in r_g
    double rp_bh2 = min_sep02/Rg_smbh;            // body2 closest approach to SMBH, in r_g
    double rp_inner_AU = min_sep12/AU;            // inner-pair closest approach, in AU

    double rp_outer_min = std::min(min_sep01, min_sep02);
    beta_actual = r_t_bin / rp_outer_min;         // achieved penetration (NA in restart)

    // surviving-orbit elements + escaper speed, for a non-terminal outcome
    if (term_event == "none") {
        if (pair_eps(p,1,2) < 0) {                        // inner binary survived
            auto [a,e] = ae(1,2); a_new = a/AU; e_new = e;
        } else if (pair_eps(p,0,1) < 0) {                 // SMBH captured body1
            auto [a,e]=ae(0,1); a_new=a/AU; e_new=e; v_esc=vinf_kms(2);   // body2 ejected
        } else if (pair_eps(p,0,2) < 0) {                 // SMBH captured body2
            auto [a,e]=ae(0,2); a_new=a/AU; e_new=e; v_esc=vinf_kms(1);   // body1 ejected
        } else {                                          // full ionization
            v_esc = std::max(vinf_kms(1), vinf_kms(2));
        }
    }

    // ---------- human-readable summary (optional) ----------
    if (report) std::cout << "\n"
      << "========================================================\n"
      << "  inner binary [" << type1 << " + " << type2 << "]  --  SMBH encounter (" << outer << ")\n"
      << "========================================================\n"
      << "  OUTCOME    : " << outcome << "\n"
      << "  exit mode  : " << exit_mode << (exit_mode=="backstop" ? "  (WARNING: hit time limit, may be unsettled)" : "") << "\n"
      << "  sim time   : " << p.time()/(2*pi) << " yr\n"
      << "  wall time  : " << wall << " s\n"
      << "  new/surviving orbit  a = " << a_new << " AU   e = " << e_new << "\n"
      << "  achieved penetration beta_actual = " << beta_actual << "\n"
      << "========================================================\n\n";

    // ============================ CSV output line ============================
    // Columns mirror the --header list.  In RESTART mode the element-only inputs
    // (a_bin, beta, v_inf, angles, ...) print as nan; everything downstream of the
    // run (outcome, surviving orbit, final phase space) is computed identically.
    std::cout << "CSV," << M_smbh_in << "," << m1_in << "," << type1 << "," << R1_in << ","
              << m2_in << "," << type2 << "," << R2_in << "," << a_bin_in << "," << e_bin_in << ","
              << outer << "," << beta_in << "," << v_inf_in << ","
              << a_out_in << "," << e_out_in << "," << nu_out_in << ","
              << incl_in << "," << Omega_in << "," << omega_in << "," << nu_in << "," << rtol_in << ","
              << outcome << "," << exit_mode << ","
              << a_new << "," << e_new << "," << v_esc << "," << beta_actual << ","
              << rp_bh1 << "," << rp_bh2 << "," << rp_inner_AU << ","
              << R1/Rs << "," << R2/Rs << "," << (2.0*M_smbh/(C*C))/Rs << ","
              << p.time()/(2*pi) << "," << wall;

    // final phase space of each body at termination (system COM frame; AU, km/s).
    // For a (mu)TDE / COLLISION these are the positions/velocities at the moment the
    // terminal radius was crossed -> ready to seed a follow-up run via --restart.
    std::cout << std::setprecision(12);
    for (size_t i = 0; i < p.number(); ++i) {
        auto r = p.pos(i);
        auto v = p.vel(i);
        std::cout << "," << r.x/AU  << "," << r.y/AU  << "," << r.z/AU
                  << "," << v.x/kms << "," << v.y/kms << "," << v.z/kms;
    }
    std::cout << "\n";
    return 0;
}
