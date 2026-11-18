#include "spaceHub.hpp"
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace hub;
using namespace unit;
using namespace force;
using namespace callback;

#ifdef NEWTONIAN_ONLY
using Forces = Interactions<NewtonianGrav>;
#else
using Forces = Interactions<NewtonianGrav, PN1, PN2, PN2p5>;
#endif
using Solver   = methods::DefaultMethod<Forces, particles::SizeParticles>;
using Particle = Solver::Particle;
using Scalar   = Solver::Scalar;

// pairwise relative energy per reduced mass (G=1): eps<0 bound, eps>0 unbound
template<class P>
double pair_eps(P& p, size_t i, size_t j) {
    auto dr = p.pos(i) - p.pos(j);
    auto dv = p.vel(i) - p.vel(j);
    return 0.5*norm2(dv) - (p.mass(i)+p.mass(j))/norm(dr);
}

// Relative energy (per reduced mass, G=1) of body k about the COM of bodies in group {a,b}.
template<class P>
double eps_to_binary_com(P& p, size_t k, size_t a, size_t b) {
    double ma = p.mass(a), mb = p.mass(b), mab = ma + mb;
    auto rcom = (p.pos(a)*ma + p.pos(b)*mb) / mab;
    auto vcom = (p.vel(a)*ma + p.vel(b)*mb) / mab;
    auto dr = p.pos(k) - rcom;
    auto dv = p.vel(k) - vcom;
    double M = p.mass(k) + mab;
    return 0.5*norm2(dv) - M/norm(dr);
}

// Direct-capture (plunge) test.
template<class P>
bool pair_plunges(P& p, size_t i, size_t j) {
    using namespace hub::consts;
    const double K_PLUNGE = 4.0;
    const double GUARD    = 1000.0;
    auto dr = p.pos(i) - p.pos(j);
    auto dv = p.vel(i) - p.vel(j);
    double r  = norm(dr);
    double v2 = norm2(dv);
    double M  = p.mass(i) + p.mass(j);
    double Rg = M / (C*C);
    if (r > GUARD * Rg) return false;
    double L2 = norm2(cross(dr, dv));
    double E  = 0.5*v2 - M/r;
    double rp;
    if (std::abs(E) < 1e-30) {
        rp = L2 / (2.0*M);
    } else {
        rp = (-M + std::sqrt(std::max(0.0, M*M + 2.0*E*L2))) / (2.0*E);
        if (rp < 0) rp = (-M - std::sqrt(std::max(0.0, M*M + 2.0*E*L2))) / (2.0*E);
    }
    return rp < K_PLUNGE * Rg;
}

// Unified terminal-event test for a pair (i,j).
// returns: 0 none, 1 TDE (sets victim), 2 PLUNGE (sets victim), 3 COLLISION
template<class P>
int pair_event(P& p, size_t i, size_t j, const bool* ext, size_t& victim) {
    using namespace hub::consts;
    double sep = distance(p.pos(i), p.pos(j));
    double mi = p.mass(i), mj = p.mass(j);
    double Ri = p.radius(i), Rj = p.radius(j);
    bool ei = ext[i], ej = ext[j];

    double r_surf = Ri + Rj;

    double rt_i = ei ? Ri * std::pow((mi + mj)/mi, 1.0/3.0) : -1.0;
    double rt_j = ej ? Rj * std::pow((mi + mj)/mj, 1.0/3.0) : -1.0;
    size_t vic; double r_tide;
    if (rt_i >= rt_j) { vic = i; r_tide = rt_i; } else { vic = j; r_tide = rt_j; }

    double Rg = (mi + mj) / (C*C);
    double r_plunge = 4.0 * Rg;
    bool tde_possible = (r_tide > r_surf) && (r_tide > r_plunge);

    if (tde_possible) {
        if (sep < r_tide) { victim = vic; return 1; }
        return 0;
    }
    if (r_surf >= r_plunge) {
        if (sep < r_surf) return 3;
        return 0;
    }
    const double PLUNGE_CHECK = 50.0;
    if (sep < PLUNGE_CHECK * Rg && pair_plunges(p, i, j)) {
        victim = (mi <= mj ? i : j); return 2;
    }
    return 0;
}

int main(int argc, char** argv) {
    // ---------- defaults ----------
    double M_smbh_in = 1e6;
    double m1_in = 10, m2_in = 1;
    std::string type1 = "bh", type2 = "star";
    double R1_in = -1, R2_in = -1;
    double a_bin_in = 1, e_bin_in = 0.7;
    double incl_in = 0, Omega_in = 0, omega_in = 0, nu_in = 0;
    std::string mode = "full";
    double rtol_in = 1e-13;
    bool report = true;

    // ---- outer-orbit controls ----
    std::string outer = "hyperbolic";     // hyperbolic | elliptic
    // hyperbolic params:
    double v_inf_in = 100, beta_in = 1;
    double n_dyn_in = 1e4;                 // backstop = n_dyn * t_dyn
    // elliptic params:
    double a_out_in = 100.0;               // outer semi-major axis [AU]
    double e_out_in = 0.5;                 // outer eccentricity
    double nu_out_in = M_PI;               // outer true anomaly at start (pi = apocenter)
    double n_orb_in = 1000.0;              // backstop = n_orb * P_out

    // ---------- parse ----------
    auto usage = [&](){
        std::cerr << "usage: " << argv[0] << " [--flag value ...]\n"
          << "  --M_smbh [Msun]   --m1 [Msun]   --m2 [Msun]\n"
          << "  --type1 [star|wd|ns|bh]  --type2 [...]\n"
          << "  --R1 [Rsun,-1=auto]      --R2 [Rsun,-1=auto]\n"
          << "  --a_bin [AU]  --e_bin\n"
          << "  --incl [rad]  --Omega [rad]  --omega [rad]  --nu [rad]   (inner binary)\n"
          << "  --outer [hyperbolic|elliptic]\n"
          << "    hyperbolic: --beta [rt/rp]  --v_inf [km/s]  --n_dyn [t_dyn units, def 1e4]\n"
          << "    elliptic:   --a_out [AU]    --e_out         --nu_out [rad, def pi=apo]  --n_orb [def 1000]\n"
          << "  --mode [full|classify]   --rtol\n"
          << "  (params for the non-selected outer mode are ignored; --help prints this)\n";
    };
    for (int a = 1; a < argc; ++a) if (std::string(argv[a]) == "--header") {
        std::cout << "CSV,M_smbh,m1,type1,R1,m2,type2,R2,a_bin,e_bin,outer,beta,v_inf,a_out,e_out,nu_out,"
                     "incl,Omega,omega,nu,rtol,"
                     "outcome,exit_mode,a_new_AU,e_new,v_esc_kms,beta_actual,"
                     "rp_body1_Rg,rp_body2_Rg,rp_inner_AU,R1_used_Rsun,R2_used_Rsun,Rsmbh_Rsun,t_sim_yr,wall_s\n";
        return 0;
    }
    for (int a = 1; a < argc; ++a) {
        std::string key = argv[a];
        if (key == "--help" || key == "-h") { usage(); return 0; }
        if (a + 1 >= argc) { std::cerr << "missing value for " << key << "\n"; usage(); return 1; }
        std::string val = argv[++a];
        if      (key == "--M_smbh") M_smbh_in = std::stod(val);
        else if (key == "--m1")     m1_in     = std::stod(val);
        else if (key == "--m2")     m2_in     = std::stod(val);
        else if (key == "--type1")  type1     = val;
        else if (key == "--type2")  type2     = val;
        else if (key == "--R1")     R1_in     = std::stod(val);
        else if (key == "--R2")     R2_in     = std::stod(val);
        else if (key == "--a_bin")  a_bin_in  = std::stod(val);
        else if (key == "--e_bin")  e_bin_in  = std::stod(val);
        else if (key == "--incl")   incl_in   = std::stod(val);
        else if (key == "--Omega")  Omega_in  = std::stod(val);
        else if (key == "--omega")  omega_in  = std::stod(val);
        else if (key == "--nu")     nu_in     = std::stod(val);
        else if (key == "--outer")  outer     = val;
        else if (key == "--v_inf")  v_inf_in  = std::stod(val);
        else if (key == "--beta")   beta_in   = std::stod(val);
        else if (key == "--n_dyn")  n_dyn_in  = std::stod(val);
        else if (key == "--a_out")  a_out_in  = std::stod(val);
        else if (key == "--e_out")  e_out_in  = std::stod(val);
        else if (key == "--nu_out") nu_out_in = std::stod(val);
        else if (key == "--n_orb")  n_orb_in  = std::stod(val);
        else if (key == "--mode")   mode      = val;
        else if (key == "--rtol")   rtol_in   = std::stod(val);
        else if (key == "--report") report    = (val=="1"||val=="true"||val=="on");
        else { std::cerr << "unknown flag: " << key << "\n"; usage(); return 1; }
    }
    bool classify = (mode == "classify");

    auto valid = [](const std::string& t){ return t=="star"||t=="wd"||t=="ns"||t=="bh"; };
    if (mode != "full" && mode != "classify") { std::cerr << "mode must be full|classify\n"; return 1; }
    if (!valid(type1) || !valid(type2))       { std::cerr << "type must be star|wd|ns|bh\n"; return 1; }
    if (outer != "hyperbolic" && outer != "elliptic") { std::cerr << "outer must be hyperbolic|elliptic\n"; return 1; }
    bool is_ell = (outer == "elliptic");
    if (is_ell && !(e_out_in >= 0.0 && e_out_in < 1.0)) { std::cerr << "elliptic requires 0 <= e_out < 1\n"; return 1; }

    using namespace hub::consts;                       // G, C, pi

    // ---------- convert to internal units ----------
    Scalar M_smbh = M_smbh_in * Ms;
    Scalar m1     = m1_in     * Ms;
    Scalar m2     = m2_in     * Ms;
    Scalar a_bin  = a_bin_in  * AU;
    Scalar e_bin  = e_bin_in;
    Scalar incl   = incl_in;
    Scalar Omega  = Omega_in;
    Scalar omega  = omega_in;
    Scalar nu     = nu_in;
    Scalar v_inf  = v_inf_in  * kms;
    Scalar beta   = beta_in;
    Scalar a_out  = a_out_in  * AU;
    Scalar e_out  = e_out_in;
    Scalar nu_out = nu_out_in;

    auto star_type = [](const std::string& t) {
        if (t == "star") return stellar::StarType::STAR;
        if (t == "wd")   return stellar::StarType::WD;
        if (t == "ns")   return stellar::StarType::NS;
        return stellar::StarType::BH;
    };
    auto radius = [&](const std::string& type, double R_in, Scalar m) -> Scalar {
        if (R_in < 0) return stellar::stellar_radius(star_type(type), m);
        return R_in * Rs;
    };
    Scalar R1 = radius(type1, R1_in, m1);
    Scalar R2 = radius(type2, R2_in, m2);

    bool ext[3] = { false, (type1 != "bh"), (type2 != "bh") };

    // ---------- derived scales ----------
    Scalar m_bin   = m1 + m2;
    Scalar M_tot   = M_smbh + m_bin;
    Scalar r_t_bin = a_bin * std::pow(M_smbh / m_bin, 1.0/3.0);
    Scalar Rg_smbh = M_tot / (C*C);

    // outer-orbit-dependent scales
    Scalar rp_out, b_out, r_start, t_backstop, P_out = 0;
    if (!is_ell) {
        // hyperbolic: pericenter from beta, impact parameter from (rp,v_inf)
        rp_out  = r_t_bin / beta;
        b_out   = std::sqrt(rp_out*rp_out + 2.0*M_tot*rp_out/(v_inf*v_inf));
        r_start = std::max<Scalar>(100.0*r_t_bin,
                  std::max<Scalar>(100.0*rp_out, 5.0*b_out));
        Scalar t_dyn = 2.0 * pi * std::sqrt(r_start*r_start*r_start / M_tot);
        t_backstop = n_dyn_in * t_dyn;
    } else {
        // elliptic: pericenter and start radius from (a_out,e_out); start at nu_out
        rp_out  = a_out * (1.0 - e_out);
        b_out   = 0;                                   // not used
        // start separation for reporting: apocenter if nu_out=pi
        r_start = a_out * (1.0 + e_out);               // apocenter (nu_out=pi default)
        P_out   = 2.0 * pi * std::sqrt(a_out*a_out*a_out / M_tot);
        t_backstop = n_orb_in * P_out;
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
    Particle smbh{M_smbh, 2.0*M_smbh/(C*C)};
    Particle b1  {m1, R1};
    Particle b2  {m2, R2};

    // inner binary: oriented by the CLI angles (incl = mutual inclination vs outer plane)
    auto inner = orbit::Elliptic(m1, m2, a_bin, e_bin, incl, Omega, omega, nu);
    orbit::move_particles(inner, b2);
    orbit::move_to_COM_frame(b1, b2);

    // outer orbit: fixed in reference plane (angles 0,0,0); branch on outer mode
    if (!is_ell) {
        auto outer_orb = orbit::Hyperbolic(M_smbh, m_bin, v_inf, b_out,
                                           0.0, 0.0, 0.0, r_start, orbit::Hyper::in);
        orbit::move_particles(outer_orb, b1, b2);
    } else {
        // NOTE: verify this SpaceHub Elliptic signature matches the inner one:
        //   Elliptic(m_of_body_a, m_of_body_b, a, e, incl, Omega, omega, nu)
        auto outer_orb = orbit::Elliptic(M_smbh, m_bin, a_out, e_out,
                                         0.0, 0.0, 0.0, nu_out);
        orbit::move_particles(outer_orb, b1, b2);
    }
    orbit::move_to_COM_frame(smbh, b1, b2);

    Solver solver{0, smbh, b1, b2};

    // ---------- solver args + stop conditions ----------
    Solver::RunArgs args;
    args.rtol = rtol_in;

    std::string term_event = "none";
    std::string exit_mode  = "backstop";
    size_t term_victim = 99, term_i = 99, term_j = 99;

    double min_sep01 = 1e300, min_sep02 = 1e300, min_sep12 = 1e300;
    args.add_operation([&min_sep01,&min_sep02,&min_sep12](auto& p, auto h) {
        double s01 = distance(p.pos(0), p.pos(1));
        double s02 = distance(p.pos(0), p.pos(2));
        double s12 = distance(p.pos(1), p.pos(2));
        if (s01 < min_sep01) min_sep01 = s01;
        if (s02 < min_sep02) min_sep02 = s02;
        if (s12 < min_sep12) min_sep12 = s12;
    });

    // time backstop (both modes; value differs by mode)
    args.add_stop_condition((Scalar)t_backstop);

    // unified terminal-event stop (both modes)
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

    // escape stop: HYPERBOLIC ONLY (a bound outer orbit never escapes)
    if (!is_ell) {
        Scalar r_far = 2.0 * r_start;
        args.add_stop_condition([r_far,&exit_mode](auto& p, auto h) -> bool {
            size_t n = p.number();
            for (size_t k = 0; k < n; ++k) {
                double mo = 0;
                typename std::remove_reference_t<decltype(p)>::Vector rc{0,0,0}, vc{0,0,0};
                for (size_t j = 0; j < n; ++j) if (j != k) {
                    mo += p.mass(j); rc += p.pos(j)*p.mass(j); vc += p.vel(j)*p.mass(j);
                }
                rc /= mo; vc /= mo;
                auto dr = p.pos(k) - rc;
                auto dv = p.vel(k) - vc;
                if (norm(dr) > r_far && dot(dr, dv) > 0) { exit_mode = "escape"; return true; }
            }
            return false;
        });
    }

    // ---------- outputs ----------
    std::ostringstream tag;
    tag << "bss_" << outer << "_" << type1 << m1_in << "_" << type2 << m2_in
        << "_i" << incl_in << "_" << mode;
    if (!classify) {
        args.add_operation(DefaultWriter(tag.str() + ".txt"));
        args.add_operation(EnergyErrWriter(tag.str() + ".err"));
    }

    // ---------- run (timed) ----------
    auto t0 = std::chrono::steady_clock::now();
    solver.run(args);
    auto t1 = std::chrono::steady_clock::now();
    double wall = std::chrono::duration<double>(t1 - t0).count();

    // ---------- classify ----------
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

    // ---------- derived quantities ----------
    const double NA = std::nan("");
    auto ae = [&](size_t i, size_t j) {
        double u = p.mass(i) + p.mass(j);
        auto [a, e] = orbit::calc_a_e(u, p.pos(i)-p.pos(j), p.vel(i)-p.vel(j));
        return std::make_pair((double)a, (double)e);
    };
    auto vinf_kms = [&](size_t k) { return norm(p.vel(k)) / kms; };

    double a_new = NA, e_new = NA;
    double v_esc = NA;
    double beta_actual = NA;
    double rp_bh1 = min_sep01/Rg_smbh;
    double rp_bh2 = min_sep02/Rg_smbh;
    double rp_inner_AU = min_sep12/AU;

    double rp_outer_min = std::min(min_sep01, min_sep02);
    beta_actual = r_t_bin / rp_outer_min;

    if (term_event == "none") {
        if (pair_eps(p,1,2) < 0) {
            auto [a,e] = ae(1,2); a_new = a/AU; e_new = e;
        } else {
            if (pair_eps(p,0,1) < 0) { auto [a,e]=ae(0,1); a_new=a/AU; e_new=e; v_esc=vinf_kms(2); }
            else if (pair_eps(p,0,2) < 0) { auto [a,e]=ae(0,2); a_new=a/AU; e_new=e; v_esc=vinf_kms(1); }
            else { v_esc = std::max(vinf_kms(1), vinf_kms(2)); }
        }
    }

    // ---------- formatted report (optional) ----------
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

    // ---------- CSV line ----------
    std::cout << "CSV," << M_smbh_in << "," << m1_in << "," << type1 << "," << R1_in << ","
              << m2_in << "," << type2 << "," << R2_in << "," << a_bin_in << "," << e_bin_in << ","
              << outer << "," << beta_in << "," << v_inf_in << ","
              << a_out_in << "," << e_out_in << "," << nu_out_in << ","
              << incl_in << "," << Omega_in << "," << omega_in << "," << nu_in << "," << rtol_in << ","
              << outcome << "," << exit_mode << ","
              << a_new << "," << e_new << "," << v_esc << "," << beta_actual << ","
              << rp_bh1 << "," << rp_bh2 << "," << rp_inner_AU << ","
              << R1/Rs << "," << R2/Rs << "," << (2.0*M_smbh/(C*C))/Rs << ","
              << p.time()/(2*pi) << "," << wall << "\n";
    return 0;
}
