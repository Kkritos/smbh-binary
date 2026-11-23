# Theory Behind the smbh-binary Project

*A self-contained tutorial on the physics and numerical methods behind the simulations in this repository: encounters of a
stellar binary (a black hole plus a star, or two stars) with a massive black hole. It covers what happens to the binary,
when the star is torn apart by its companion black hole (a "micro-tidal disruption event", μTDE), how the debris falls
back to produce two flares, and how all of this is computed with SpaceHub and the analysis code.*

**Who this is for.** Anyone with a solid grasp of calculus (derivatives, integrals, vectors) and high-school physics
(Newton's laws, energy). Every new idea is introduced from scratch, with a worked numerical example and its units.

**How to read it.**
- Part I builds the toolkit: units, energy, orbits.
- Part II covers the astrophysical objects: black holes, stars, binaries, and the three-body encounter.
- Part III covers tidal disruption and the debris.
- Part IV covers the numerical methods: how SpaceHub integrates orbits, and how the analysis code turns thousands of
  simulations into cross sections and rates.
- Part V connects the simulations to the real Universe: event rates.

Boxes marked **Worked example** contain numbers; boxes marked **In plain words** restate the physics without equations.
Math is written in LaTeX and renders on GitHub.

---

## Contents

**Part 0: How to think like an astrophysicist**

**Part I: Foundations**
1. Units: CGS, astronomical units, and the code units used by SpaceHub
2. Energy, specific energy, and angular momentum
3. The two-body problem: Kepler orbits
4. Hyperbolic encounters, impact parameter, and gravitational focusing

**Part II: The objects and the encounter**

5. Black holes
6. Beyond Newton: post-Newtonian (PN) gravity
7. Stars
8. Binaries
9. The three-body encounter: binary meets massive black hole

**Part III: Tidal disruption and debris**

10. Tidal forces and the tidal radius
11. Tidal disruption: frozen-in energies and the fallback rate
12. The μTDE inside an SMBH encounter: two flares and the delay
13. Ballistic debris physics: how the debris code works

**Part IV: Numerical methods**

14. Integrating orbits: from Euler to symplectic methods
15. Why few-body problems are hard: close encounters and regularization
16. SpaceHub: the integrator used in this project
17. The driver `smbh_binary.cpp`: setup, stopping conditions, classification
18. Monte Carlo scattering experiments and cross sections

**Part V: From simulations to the Universe**

19. Event rates: loss-cone theory
20. The f_μ fit and the population average

Appendices: A. Symbols · B. Constants · C. References · D. Changelog

---

# Part 0: How to think like an astrophysicist

Astrophysicists study systems they can't touch, over times they can't wait for, with physics spanning many orders of
magnitude. The habit that makes this possible is **estimating before calculating**. Here is the toolkit, applied
throughout this document.

1. **Know the scales.** Before any calculation, ask: what is the typical size, mass, speed and time? For our encounters:
   - a binary of size $a\sim0.1$ AU;
   - a black hole of mass $M_\bullet\sim10^6\thinspace M_\odot$;
   - the binary's orbital speed $\sim300$ km/s;
   - the speed at closest approach to the black hole $\sim10^4$ km/s.
2. **Compare things that have the same units.** A number like "1000 km/s" means nothing alone. It means a lot once you
   compare it with the speed of light (300 000 km/s), or with the escape speed. Ratios of like quantities, such as $v/c$,
   $r_p/r_g$ or $v_\infty/v_{\rm orb}$, decide which physics matters. They are called **dimensionless numbers**.
3. **Compare timescales.** The fastest process usually controls what happens; much slower processes can be treated as
   frozen. *Example:* the binary orbits in ~3 days, but the encounter with the black hole takes a fraction of a day near closest
   approach. Comparable timescales mean the binary is strongly disturbed.
4. **Find the balance.** Most structures are a balance between two effects: gravity versus pressure (a star), the black
   hole's tidal pull versus the binary's own gravity (the tidal radius), kinetic versus potential energy (bound or not).
   Write both effects with scales plugged in, set them equal, and you have the key scale.
5. **Use conservation laws.** Energy, momentum and angular momentum often give the answer without solving any equations
   of motion. The impact-parameter formula (§4) and the fallback time (§11) come from conservation alone.
6. **Check limits.** Does your formula behave sensibly when a mass goes to zero, or a speed to infinity?
7. **Check the numbers.** Is the answer reasonable? A star moving faster than light, or a black hole smaller than an
   atom, means a mistake.

**Thinking about simulations.** A simulation is an experiment on a computer. Treat it like one:
- **Calibrate:** test the code on problems with known answers (a two-body orbit should stay a perfect ellipse).
- **Measure the error:** track a quantity that should be conserved, such as total energy.
- **Test convergence:** change a numerical choice (tolerance, sampling range) and check the answer doesn't move.
- **Count statistics:** a fraction measured from $N$ trials has an uncertainty of about $\sqrt{f(1-f)/N}$.
- **Never trust a single run:** chaotic three-body encounters must be studied statistically.

---

# Part I: Foundations

## 1. Units

### 1.1 Why units matter

A physical quantity is a number **times a unit**: "3 km" and "3 AU" differ by a factor of $1.5\times10^8$. Every equation
must have the same units on both sides. Checking that (dimensional analysis) catches most algebra mistakes.

### 1.2 CGS and astronomical units

Astrophysicists mostly use **CGS** units: centimetres, grams, seconds.

| Quantity | CGS unit | In base units | SI equivalent |
|---|---|---|---|
| length | cm | cm | 0.01 m |
| mass | g | g | 0.001 kg |
| time | s | s | s |
| speed | cm/s | cm s⁻¹ | 0.01 m/s |
| force | dyne | g cm s⁻² | 10⁻⁵ N |
| energy | erg | g cm² s⁻² | 10⁻⁷ J |
| power (luminosity) | erg/s | g cm² s⁻³ | 10⁻⁷ W |

On astronomical scales, natural units are built from the Sun and the Earth's orbit:

| Unit | Value (CGS) | Meaning |
|---|---|---|
| $M_\odot$ (solar mass) | 1.989×10³³ g | mass of the Sun |
| $R_\odot$ (solar radius) | 6.957×10¹⁰ cm | radius of the Sun |
| AU (astronomical unit) | 1.496×10¹³ cm = 215 $R_\odot$ | Earth–Sun distance |
| pc (parsec) | 3.086×10¹⁸ cm = 206 265 AU | typical distance between stars |
| yr (year) | 3.156×10⁷ s | |
| $G$ | 6.674×10⁻⁸ cm³ g⁻¹ s⁻² | gravitational constant |
| $c$ | 2.998×10¹⁰ cm/s = 299 792 km/s | speed of light |
| $GM_\odot$ | 1.327×10²⁶ cm³ s⁻² | appears in every orbit formula |

### 1.3 The code units of SpaceHub (and of this project)

Computers handle numbers near 1 best, so the code chooses units in which typical quantities are of order unity, and sets
**$G = 1$**:

| Code unit | Value | Why |
|---|---|---|
| mass | $1\thinspace M_\odot$ | |
| length | 1 AU | |
| speed | $V_1 = \sqrt{GM_\odot/\mathrm{AU}} = 29.78$ km/s | the Earth's orbital speed |
| time | $\sqrt{\mathrm{AU}^3/GM_\odot} = 5.02\times10^6$ s $= 1\thinspace\mathrm{yr}/2\pi = 0.159$ yr | |
| speed of light | $c = 299\thinspace 792/29.78 = 1.006\times10^4$ in code units | |

> **Why the time unit is 1 yr / 2π.** Kepler's third law with $G=1$ reads $P = 2\pi\sqrt{a^3/M}$. For the Earth
> ($a = 1$, $M = 1$), $P = 2\pi$ code time units. That must equal 1 year, so one code time unit is $1/(2\pi)$ yr.
> That's why the code converts times to years by dividing by $2\pi$ (`p.time()/(2*pi)` in `smbh_binary.cpp`).

> **Worked example: converting a speed.** The binary's orbital speed for $a = 0.1$ AU and $m_{\rm bin} = 11\thinspace M_\odot$ is,
> in code units, $v = \sqrt{m/a} = \sqrt{110} = 10.5$. In km/s: $10.5\times29.78 = 312$ km/s.
> **Units check:** $\sqrt{GM/a}$ has units $\sqrt{(\mathrm{cm^3\thinspace g^{-1}\thinspace s^{-2}})(\mathrm{g})/\mathrm{cm}} = \mathrm{cm/s}$ ✓.

> **Worked example: a gravitational radius.** With $G = 1$, the gravitational radius $r_g = GM/c^2$ becomes $M/C^2$ with
> $C = 1.006\times10^4$. For $M = 10^6$: $r_g = 10^6/(1.006\times10^4)^2 = 9.87\times10^{-3}$ AU $= 1.48\times10^{11}$ cm.

---

## 2. Energy, specific energy, and angular momentum

### 2.1 What is energy?

Energy is a number you can compute from the positions and velocities of a system that **stays constant as the system
evolves**, provided nothing outside acts on it. It's a bookkeeping device: it can change form (motion into height, motion
into heat), but the total doesn't change. Deep down, energy conservation follows from the fact that the laws of physics
don't change with time (Noether's theorem).

For a body of mass $m$ moving at speed $v$ at distance $r$ from a point mass $M$:

$$E = \underbrace{\tfrac12mv^2}_{\text{kinetic}} \ \underbrace{-\ \frac{GMm}{r}}_{\text{potential}} .$$

- **Kinetic energy** is energy of motion, always positive.
- **Gravitational potential energy** is negative. The zero is chosen at infinite separation, and bringing bodies closer
  *releases* energy (they speed up). So a body deep in a gravitational well has a large negative potential energy.

> **Why the potential energy is −GMm/r.** The gravitational force is $F = GMm/r^2$, directed inward. The work you must do
> to pull $m$ from $r$ out to infinity is $\int_r^\infty\frac{GMm}{r'^2}dr' = \frac{GMm}{r}$. Ending with zero potential
> energy at infinity means you started with $-GMm/r$.

**Units:** erg = g cm² s⁻². Check: $GMm/r$ has units (cm³ g⁻¹ s⁻²)(g)(g)/cm = g cm² s⁻² ✓.

### 2.2 Specific energy: energy per unit mass

In gravity, the motion of a small body doesn't depend on its own mass. A feather and a hammer fall the same way (the
equivalence principle). So it's natural to divide by $m$ and work with the **specific orbital energy**:

$$\varepsilon \equiv \frac{E}{m} = \frac12v^2 - \frac{GM}{r}\qquad[\mathrm{erg/g} = \mathrm{cm^2/s^2} = (\mathrm{km/s})^2\times10^{10}].$$

- **Specific energy has units of velocity squared.** You can quote it in (km/s)², which makes it easy to picture.
- **Its sign tells you whether the orbit is bound:**
  - $\varepsilon\lt 0$: **bound**. The body can never escape to infinity, because it would need positive kinetic energy there.
  - $\varepsilon\gt 0$: **unbound**. It escapes, with leftover speed $v_\infty = \sqrt{2\varepsilon}$ at infinity.
  - $\varepsilon = 0$: exactly **parabolic**, the boundary case.

**For two bodies of comparable mass**, use their *relative* position $\mathbf r$, relative velocity $\mathbf v$, and the
total mass $M = m_1+m_2$: $\varepsilon = \tfrac12v^2 - GM/r$. This is exactly what `pair_eps` computes in
`smbh_binary.cpp`, and the classifier uses its sign to decide which pairs remain bound.

> **Worked example: the Earth.** $v = 29.78$ km/s, $r = 1$ AU, $GM_\odot/r = V_1^2 = (29.78)^2 = 887$ (km/s)².
> So $\varepsilon = \tfrac12(887) - 887 = -443$ (km/s)² $= -4.4\times10^{12}$ erg/g. Negative: the Earth is bound.
> To escape it would need $\varepsilon\ge0$, i.e. $v\ge\sqrt{2}\times29.78 = 42.1$ km/s, the solar escape speed at 1 AU.

> **Worked example: is the star bound to the 10 M☉ black hole?** Take a separation of 0.1 AU and a relative speed of
> 300 km/s. Then $GM/r = 10\times887/0.1 = 88\thinspace 700$ (km/s)² and $\tfrac12v^2 = 45\thinspace 000$ (km/s)².
> $\varepsilon = -43\thinspace 700$ (km/s)² < 0: bound.

**Specific energy relative to a pair.** To ask whether body $k$ is bound to a *binary* $\lbrace a,b\rbrace$, use $k$'s position and
velocity relative to the binary's **centre of mass**, and the total mass. That's `eps_to_binary_com` in the code.

### 2.3 Angular momentum

The **specific angular momentum** $\mathbf h = \mathbf r\times\mathbf v$ (units cm²/s) measures how much "sideways" motion
a body has. Gravity between two bodies acts along the line joining them, so it exerts no torque, and $\mathbf h$ is
conserved. Two consequences:
1. **The orbit lies in a fixed plane**, the plane perpendicular to $\mathbf h$.
2. **The body sweeps equal areas in equal times** (Kepler's second law): $dA/dt = h/2$.

The magnitude $h = r v_\perp$ controls how close the body can get. At closest approach (pericentre) all the motion is
sideways, so $h = r_pv_p$.

---

## 3. The two-body problem: Kepler orbits

### 3.1 Reduction to one body

Two bodies $m_1$, $m_2$ attract each other. Their centre of mass moves at constant velocity, so we can sit in its frame
(SpaceHub does this: `move_to_COM_frame`). The relative vector $\mathbf r = \mathbf r_2-\mathbf r_1$ obeys

$$\ddot{\mathbf r} = -\frac{G(m_1+m_2)}{r^3}\mathbf r,$$

which is the motion of a single body around a fixed mass $M = m_1+m_2$. **Everything about a two-body orbit depends only on
the total mass.**

### 3.2 Conic sections and orbital elements

Solving this equation (Kepler, Newton) gives a **conic section**:

$$r(\nu) = \frac{p}{1+e\cos\nu},\qquad p = \frac{h^2}{GM}.$$

| $e$ | Shape | Energy |
|---|---|---|
| 0 | circle | $\varepsilon\lt 0$ |
| 0 < e < 1 | ellipse | $\varepsilon\lt 0$ |
| 1 | parabola | $\varepsilon = 0$ |
| > 1 | hyperbola | $\varepsilon\gt 0$ |

Six numbers, the **orbital elements**, specify an orbit completely:
- **semi-major axis** $a$: the size;
- **eccentricity** $e$: the shape;
- **inclination** $i$, **longitude of the ascending node** $\Omega$, **argument of pericentre** $\omega$: the orientation in space;
- **true anomaly** $\nu$: where the body is along the orbit right now (the angle from pericentre).

These are exactly the inputs `--a_bin --e_bin --incl --Omega --omega --nu` of the simulation driver.

**Key relations** (all derivable from energy and angular-momentum conservation):

$$\varepsilon = -\frac{GM}{2a},\qquad h^2 = GMa(1-e^2),\qquad r_p = a(1-e),\qquad r_a = a(1+e),\qquad P = 2\pi\sqrt{\frac{a^3}{GM}}.$$

The **vis-viva equation** gives the speed anywhere on the orbit: $v^2 = GM\left(\frac2r-\frac1a\right)$.

> **Worked example: the BH–star binary.** $m_{\rm bin} = 11\thinspace M_\odot$, $a = 0.1$ AU, $e = 0.7$.
> - Period: $P = \sqrt{a^3/m}$ yr $= \sqrt{10^{-3}/11} = 9.5\times10^{-3}$ yr $= 3.5$ days.
> - Pericentre: $r_p = 0.03$ AU $= 6.4\thinspace R_\odot$; apocentre $r_a = 0.17$ AU.
> - Orbital energy per unit (reduced) mass: $\varepsilon = -GM/2a = -(11/0.2)(887)$ (km/s)² $= -4.9\times10^4$ (km/s)².

### 3.3 Mean anomaly, eccentric anomaly, and Kepler's equation

A body on an eccentric orbit moves fast near pericentre and slowly near apocentre, so the true anomaly $\nu$ does *not*
grow uniformly in time. What does grow uniformly is the **mean anomaly** $\mathcal M = 2\pi t/P$. The link goes through the
**eccentric anomaly** $E$:

$$\mathcal M = E - e\sin E\qquad(\text{Kepler's equation}),\qquad \tan\frac\nu2 = \sqrt{\frac{1+e}{1-e}}\tan\frac E2 .$$

Kepler's equation has no closed-form solution for $E$, so it is solved numerically by **Newton's method**:
$E_{k+1} = E_k - (E_k - e\sin E_k - \mathcal M)/(1-e\cos E_k)$, which converges in a handful of iterations.

**Where this is used.** The sweep generator `smbh_bhstar_sweep.py` picks the binary's starting phase **uniformly in mean
anomaly** (uniform in time), then converts to $\nu$ with exactly this Newton iteration (`nu_from_uniform_time`). Why
uniform in time? If you photograph a binary at a random moment, you're more likely to catch it near apocentre, where it
spends more time. Choosing $\nu$ uniformly would over-represent pericentre and bias the results.

### 3.4 From position and velocity to orbital elements

Given $\mathbf r$ and $\mathbf v$ at one instant (the *osculating* orbit: the Kepler orbit the pair would follow if all
other forces vanished), compute $\varepsilon = \tfrac12v^2-GM/r$ and $h = |\mathbf r\times\mathbf v|$. Then
$a = -GM/2\varepsilon$ and $e = \sqrt{1+2\varepsilon h^2/(GM)^2}$.

The pericentre distance follows from the quadratic obtained by setting $v_r = 0$ in the energy equation,
$\varepsilon = h^2/(2r_p^2) - GM/r_p$:

$$r_p = \frac{-GM+\sqrt{(GM)^2+2\varepsilon h^2}}{2\varepsilon}\qquad(\to h^2/2GM\ \text{for}\ \varepsilon\to0).$$

This is the formula in `pair_plunges` in `smbh_binary.cpp`. It predicts whether two bodies are *about to* pass closer
than the capture radius, even before they get there.

---

## 4. Hyperbolic encounters, impact parameter, and gravitational focusing

### 4.1 The setup

The binary approaches the massive black hole from far away, with relative speed $v_\infty$ "at infinity". If gravity were
switched off, it would miss the black hole by the **impact parameter** $b$. Gravity bends the path toward the black hole,
so the actual closest approach $r_p$ is *smaller* than $b$.

### 4.2 The impact parameter–pericentre relation

Two conservation laws link $b$ and $r_p$:
- **angular momentum:** $h = bv_\infty = r_pv_p$;
- **energy:** $\tfrac12v_\infty^2 = \tfrac12v_p^2 - GM/r_p$.

Eliminate $v_p = bv_\infty/r_p$:

$$\boxed{b^2 = r_p^2 + \frac{2GM\thinspace r_p}{v_\infty^2} = r_p^2 + r_{\rm foc}\thinspace r_p,\qquad r_{\rm foc}\equiv\frac{2GM}{v_\infty^2}.}$$

This exact formula appears in `smbh_binary.cpp` (`b_out = sqrt(rp^2 + 2 M rp / v_inf^2)`) and in `analysis/bhstar_mutde/pkg/scales.py`.

### 4.3 Gravitational focusing

The cross section for passing within $r_p$ is the area of the disk of radius $b$:

$$\Sigma(\lt r_p) = \pi b^2 = \pi r_p^2\left(1+\frac{2GM}{r_pv_\infty^2}\right) = \pi r_p^2\left(1+\frac{v_{\rm esc}^2(r_p)}{v_\infty^2}\right).$$

The second term is **gravitational focusing**: gravity pulls in trajectories that would otherwise miss. When $v_\infty$ is
small compared with the speed at pericentre (always true here), focusing dominates by an enormous factor.

> **Worked example: focusing for an SMBH encounter.** $M = 10^6\thinspace M_\odot$, binary with $a = 0.1$ AU,
> $v_\infty = 0.01\thinspace v_{\rm orb} = 3.1$ km/s (a cold, slow approach).
> - $r_{\rm foc} = 2GM/v_\infty^2 = 2\times10^6\times887/3.1^2\approx1.8\times10^8$ AU. That's an enormous length, about 900 pc.
> - For $r_p = 45$ AU (10 binary tidal radii; §9), $b = \sqrt{45^2 + 1.8\times10^8\times45}\approx9\times10^4$ AU ≈ 0.4 pc.
>
> Gravity funnels binaries arriving within 0.4 pc down to within 45 AU of the black hole. The cross section is enhanced by
> $(b/r_p)^2\approx4\times10^6$.

**In plain words:** a slow-moving binary has plenty of time to be pulled toward the black hole, so even binaries aimed far
off end up passing close. That's why the rates of close encounters with massive black holes are much larger than the
geometric size of the target would suggest.

---

# Part II: The objects and the encounter

## 5. Black holes

### 5.1 Escape speed and the Schwarzschild radius

From §2, the escape speed from radius $r$ is $v_{\rm esc} = \sqrt{2GM/r}$. Setting $v_{\rm esc} = c$ gives a radius from
which not even light escapes. General relativity (Schwarzschild 1916) gives exactly the same answer:

$$r_s = \frac{2GM}{c^2} = 2r_g,\qquad r_g\equiv\frac{GM}{c^2}.$$

$r_g$, the **gravitational radius**, is the natural length scale of a black hole.

| Black hole | Mass | $r_g$ | $r_s$ |
|---|---|---|---|
| stellar-mass (BH1 here) | 10 $M_\odot$ | 14.8 km | 29.5 km |
| intermediate-mass (IMBH) | 10³ $M_\odot$ | 1480 km | 2.0×10⁻⁵ AU |
| Sgr A* (Milky Way centre) | 4×10⁶ $M_\odot$ | 5.9×10⁶ km | 0.08 AU (17 $R_\odot$) |
| SMBH | 10⁸ $M_\odot$ | 0.99 AU | 2 AU |

In the code, the SMBH is given radius $2r_g$ (`Particle{M_smbh, 2.0*M_smbh/(C*C)}`). Rule of thumb:
**$r_g\approx1.5$ km per solar mass**.

### 5.2 Special orbits near a black hole

General relativity changes orbits near a black hole (§6). Two radii matter:
- **ISCO (innermost stable circular orbit), $r_{\rm ISCO} = 6r_g$** for a non-spinning black hole. Inside it, circular orbits
  are unstable and gas plunges in.
- **Marginally bound orbit, $r_{\rm mb} = 4r_g$.** A body arriving on a nearly parabolic orbit ($\varepsilon\approx0$) whose
  pericentre would be inside $4r_g$ doesn't swing back out. It is captured directly (Bardeen, Press & Teukolsky 1972;
  Misner, Thorne & Wheeler 1973).

That second fact is why the code flags a **PLUNGE** when the osculating pericentre of a compact pair falls below $4r_g$
(`K_PLUNGE = 4.0` in `pair_plunges`). Far-from-parabolic orbits have slightly different capture radii; the code's
comments note that 4 $r_g$ is a geometric stand-in, adequate for classifying outcomes.

### 5.3 Swallowed whole or torn apart? The Hills mass

A star approaching a black hole is torn apart at the tidal radius $r_t = R_\star(M/m_\star)^{1/3}$ (§10). That grows only as
$M^{1/3}$, while the capture radius grows as $M$. So for heavy enough black holes, the star is **swallowed whole before it
can be disrupted**. The crossover is the **Hills mass** (Hills 1975), where $r_t = 4r_g$:

$$R_\star\left(\frac{M_H}{m_\star}\right)^{1/3} = \frac{4GM_H}{c^2}\ \Rightarrow\ M_H = \left(\frac{R_\star c^2}{4Gm_\star^{1/3}}\right)^{3/2}\approx4\times10^7\thinspace M_\odot\ \ (\text{Sun-like star}).$$

> **Worked example: tidal radius versus capture radius.**
>
> | $M$ ($M_\odot$) | $r_t$ (AU) | $r_g$ (AU) | $r_t/r_g$ | what happens |
> |---|---|---|---|---|
> | 10³ | 0.047 | 1.0×10⁻⁵ | 4700 | disruption, far from the hole |
> | 10⁶ | 0.47 | 9.9×10⁻³ | 47 | disruption |
> | 10⁸ | 2.2 | 0.99 | 2.2 | swallowed whole |
>
> The code handles this automatically. `pair_event` declares a TDE only if the tidal radius lies outside both the
> surfaces and $4r_g$ (`tde_possible`). Otherwise the star plunges.

---

## 6. Beyond Newton: post-Newtonian (PN) gravity

### 6.1 The small parameter

Newtonian gravity is the low-speed, weak-field limit of general relativity. How "relativistic" a two-body orbit is depends
on one number:

$$\frac{v^2}{c^2}\sim\frac{GM}{rc^2} = \frac{r_g}{r}.$$

The **post-Newtonian (PN) expansion** writes the acceleration as a power series in this small number:

$$\mathbf a = \underbrace{\mathbf a_{\rm N}}_{\text{Newton}}\left[1 + \underbrace{O\left(\tfrac{v^2}{c^2}\right)}_{\rm 1PN} + \underbrace{O\left(\tfrac{v^4}{c^4}\right)}_{\rm 2PN} + \underbrace{O\left(\tfrac{v^5}{c^5}\right)}_{\rm 2.5PN} + \cdots\right].$$

The "$n{\rm PN}$" term is smaller than Newton by $(v/c)^{2n}$.

| Term | Physics | Conserves energy? |
|---|---|---|
| 1PN | relativistic corrections to the orbit: pericentre **precession** (the Mercury effect) | yes |
| 2PN | the next correction, important close in | yes |
| 2.5PN | **gravitational-wave radiation reaction**: the orbit loses energy and angular momentum to gravitational waves | **no** (energy is radiated away) |

The 2.5PN term is the first odd power, which is what makes it break time-reversal symmetry: it describes energy leaving
the system. Standard references are Blanchet (2014, Living Reviews in Relativity) and, for the two-body equations of
motion used in N-body codes, Kidder (1995). This project runs with **1PN + 2PN + 2.5PN** (`Interactions<NewtonianGrav, PN1, PN2, PN2p5>`),
and also with Newtonian gravity alone (`-DNEWTONIAN_ONLY`) for comparison.

### 6.2 How big are PN effects in our encounters?

> **Worked example: the binary's pass by the SMBH.** $M = 10^6\thinspace M_\odot$, $a = 0.1$ AU, penetration $\beta = 1.5$ (§9).
> The pericentre is $r_p = r_{t,\rm bin}/\beta = 3.0$ AU $= 304\thinspace r_g$.
> - The speed there is $v_p = \sqrt{2GM/r_p}\approx24\thinspace 000$ km/s, i.e. $v_p/c = 0.08$.
> - 1PN corrections: $r_g/r_p\approx3\times10^{-3}$, a 0.3% effect.
> - Precession per passage: $\Delta\phi\approx3\pi r_g/r_p\approx0.03$ rad, about 2°.
>
> Small, but the encounters are **chaotic** (§9): small differences grow exponentially, so a 0.3% force change can flip the
> outcome of an individual encounter. Statistically the effect is modest except for the closest passages (tens of $r_g$),
> where PN terms decide whether a body plunges or survives. The star–star campaign (Newtonian vs PN) measures exactly this.

### 6.3 Gravitational-wave inspiral (2.5PN)

A binary emitting gravitational waves shrinks. For a circular orbit, the time to merge is (Peters 1964)

$$t_{\rm GW} = \frac{5}{256}\frac{c^5a^4}{G^3m_1m_2(m_1+m_2)}.$$

> **Worked example.** For the BH–star binary ($10+1\thinspace M_\odot$, $a = 0.1$ AU), $t_{\rm GW}\approx3\times10^{11}$ yr, longer than
> the age of the Universe. Gravitational waves don't matter for the binary itself. They matter for a body captured onto a
> tight orbit around the SMBH (an **extreme-mass-ratio inspiral**, EMRI): $t_{\rm GW}\propto a^4$ falls steeply as the orbit tightens.

---

## 7. Stars

### 7.1 What holds a star up

A star is a ball of hot gas in which **pressure balances gravity** (hydrostatic equilibrium): $dP/dr = -Gm(r)\rho/r^2$.
Combined with an equation of state that links pressure to density (and temperature), this sets the star's internal
structure. Two numbers characterize a star for our purposes: its **mass** and its **radius**.

- **Mean density:** $\bar\rho_\odot = M_\odot/(\tfrac43\pi R_\odot^3) = 1.4$ g/cm³.
- **Dynamical time:** $t_{\rm dyn} = \sqrt{R^3/GM}\approx1600$ s ≈ 27 min for the Sun. This is how fast a star can respond
  to any disturbance (free-fall and sound-crossing times are both about this long).

### 7.2 Polytropes

A useful simplified model assumes $P = K\rho^{1+1/n}$ everywhere: a **polytrope** of index $n$. Then hydrostatic
equilibrium plus Newton's law of gravity combine into the **Lane–Emden equation** for the dimensionless density profile
$\theta$ (with $\rho = \rho_c\theta^n$ and $r = \alpha\xi$):

$$\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2\frac{d\theta}{d\xi}\right) = -\theta^n,\qquad\theta(0)=1,\ \theta'(0)=0 .$$

The surface is the first zero, $\theta(\xi_1) = 0$.
- $n = 3$ ($\xi_1 = 6.90$) is a good model for a Sun-like star, with central density 54 times the mean.
- $n = 1.5$ ($\xi_1 = 3.65$) describes fully convective, low-mass stars.

**Where this is used.** The debris code (`debris_tracer.py`, `mutde_lab.py`) integrates the Lane–Emden equation for $n=3$
(with `scipy.integrate.solve_ivp`, starting from the series $\theta\approx1-\xi^2/6$ to avoid the $1/\xi$ singularity at the
centre). It then places debris particles inside the star so that each carries equal mass. Particles are densely packed in
the core and sparse in the envelope, just like the real mass.

> **Inverse-CDF sampling (how equal-mass particles are drawn).** The enclosed mass $M(\lt\xi)\propto\int_0^\xi\xi'^2\theta^n\thinspace d\xi'$
> rises from 0 to 1 (normalized). Draw a uniform random number $u\in[0,1]$ and find the radius where $M(\lt\xi) = u$: the
> radius *enclosing a fraction u of the mass*. Repeat $N$ times. Radii where lots of mass sits get proportionally more
> particles. Directions are drawn uniformly on the sphere ($\cos\theta$ uniform in [−1, 1], $\phi$ uniform in [0, 2π]).

### 7.3 Mass–radius relations

Main-sequence stars have $R\approx R_\odot(M/M_\odot)^{0.8}$ roughly. In this project the star is fixed at
$1\thinspace M_\odot$, $1\thinspace R_\odot$ for the BH–star campaign. Compact objects are far smaller:
- **white dwarfs:** about the size of the Earth, ~0.01 $R_\odot$;
- **neutron stars:** ~10 km;
- **black holes:** $2r_g$.

When no radius is given, the driver uses SpaceHub's built-in `stellar::stellar_radius` relations (§16.7), including
$R = R_\odot(m/M_\odot)^{0.75}$ for stars.

---

## 8. Binaries

### 8.1 Basic scales

For a binary of total mass $m_{\rm bin}$ and semi-major axis $a$:
- **orbital speed:** $v_{\rm orb} = \sqrt{Gm_{\rm bin}/a}$;
- **binding energy:** $E_{\rm bin} = -Gm_1m_2/2a$.

> **Worked example.** For $m_{\rm bin} = 11\thinspace M_\odot$:
>
> | $a$ (AU) | $v_{\rm orb}$ (km/s) | $P$ |
> |---|---|---|
> | 0.03 | 570 | 0.57 d |
> | 0.1 | 312 | 3.5 d |
> | 1 | 99 | 110 d |
> | 30 | 18 | 49 yr |
>
> The campaign grid spans exactly this range (13 values from 0.03 to 30 AU).

### 8.2 Hard and soft binaries

Compare the binary's orbital speed with the speed of the passer-by:
- **hard** binary: $v_{\rm orb}\gg v_\infty$. Too tightly bound for a typical encounter to break it.
- **soft** binary: $v_{\rm orb}\ll v_\infty$. Easily broken.

Heggie's law (Heggie 1975) says that on average hard binaries get harder and soft binaries get softer. The campaign
parameter $v_\infty/v_{\rm orb}$ (`frac`, from 0.001 to 10) scans from very hard to soft. (Against a massive black hole every
binary that gets close enough is eventually broken by the tide; §9.)

### 8.3 Eccentricity distribution

What eccentricities do real binaries have? If a population has been shuffled by many encounters until its energy states
are uniformly populated ("thermal equilibrium", Jeans 1919; Heggie 1975), the distribution is

$$p(e)\thinspace de = 2e\thinspace de\qquad(\text{the "thermal" distribution}).$$

High eccentricities are favoured: half the binaries have $e\gt 0.71$. The campaign doesn't draw $e$ randomly. It runs a
**grid** of 12 eccentricity values (0 to 0.99), then **weights** the results by $p(e)$ in post-processing. This is
*stratified sampling*: just as accurate for the same number of runs, and it can be re-weighted to any assumed $p(e)$
without re-running.

### 8.4 The survival cut

A BH–star binary with too small a pericentre would tear the star apart *on its own*, before any encounter. The star is
disrupted by its companion if the binary's pericentre $a(1-e)$ is inside the star's tidal radius about the BH (§10):

$$r_{t,\rm inner} = R_\star\left(\frac{m_{\rm bin}}{m_\star}\right)^{1/3} = 1\thinspace R_\odot\times11^{1/3} = 0.0103\ \text{AU}.$$

So only $e\lt e_{\max}(a) = 1 - r_{t,\rm inner}/a$ is physical. The sweep skips "born-disrupted" cells (`SURVIVAL_CUT`).

> **Worked example.** For $a = 0.03$ AU, $e_{\max} = 1 - 0.0103/0.03 = 0.66$: only 7 of the 12 eccentricity nodes survive.
> For $a = 1$ AU, $e_{\max} = 0.9897$, so the $e = 0.99$ node is skipped and 11 of the 12 survive; from $a\approx2$ AU on, all 12 do.

---

## 9. The three-body encounter: binary meets massive black hole

### 9.1 The binary's tidal radius

The black hole pulls harder on the nearer member of the binary than on the farther one. This *difference* (the tide)
competes with the binary's own gravity. Balance them (§10 shows the derivation):

$$r_{t,\rm bin} = a\left(\frac{M_\bullet}{m_{\rm bin}}\right)^{1/3}.$$

Inside about $r_{t,\rm bin}$, the black hole's tide is stronger than the binary's self-attraction, so the binary gets
pulled apart. The **penetration factor** measures how deep the encounter goes:

$$\beta = \frac{r_{t,\rm bin}}{r_p}.$$

$\beta\gt 1$ means a pass inside the tidal radius.

> **Worked example.** $M_\bullet = 10^6\thinspace M_\odot$, $a = 0.1$ AU: $r_{t,\rm bin} = 0.1\times(10^6/11)^{1/3} = 4.5$ AU. At $\beta = 1.5$,
> $r_p = 3.0$ AU.
>
> **Timescale check (Part 0, step 3).** Passing pericentre takes about $\sqrt{r_p^3/GM}\approx\sqrt{27/10^6}$ code units
> ≈ 0.005 × 0.159 yr ≈ 0.3 days. That's shorter than the binary's 3.5-day orbit, so the binary is hit *impulsively*: the
> black hole's tide acts before the binary can react.

### 9.2 Possible outcomes

After the encounter, the code classifies the result (`smbh_binary.cpp`, see §17 for the logic):

| Outcome | What happened | Physical significance |
|---|---|---|
| **FLYBY** | binary survives and leaves (possibly with changed $a$, $e$) | the common case for distant passes |
| **EXCHANGE** | binary broken; one member captured by the black hole, the other ejected | the **Hills mechanism** (§9.3) |
| **IONIZATION** | all three end up mutually unbound | rare for a massive black hole |
| **BOUND_TRIPLE** | binary still bound *and* bound to the black hole, or both members separately bound | unresolved within the time limit |
| **TDE:star2_by_bh1** | the star is torn apart by its companion black hole: a **μTDE** | the focus of the project |
| **TDE:star2_by_SMBH** | the star is torn apart by the massive black hole | ordinary TDE |
| **PLUNGE** | a body falls straight into a black hole | above the Hills mass, or a direct hit |
| **COLLISION** | two bodies touch | e.g. star–star collisions |

### 9.3 The Hills mechanism

When a binary is broken by a massive black hole, one member is captured onto a tight, highly eccentric orbit and the
other flies away at high speed (Hills 1988). Energy bookkeeping explains why.

At pericentre the binary's centre of mass moves at $v_p$. The member that escapes moves relative to the centre of mass at
$v_c = v_{\rm orb}\thinspace m_c/m_{\rm bin}$, where $m_c$ is the mass of its companion (the member that is captured). When
the binary is split, the escaping member keeps its total velocity, about $v_p + v_c$, and its specific energy relative to
the black hole changes by about

$$\Delta\varepsilon\sim v_pv_c.$$

That is enormous, because $v_p\gg v_c$. The captured member ends up bound; the other leaves with $v_{\rm ej}\approx\sqrt{2\Delta\varepsilon}$.
Writing $v_p = (2GM_\bullet\beta/r_{t,\rm bin})^{1/2}$ with $r_{t,\rm bin} = a(M_\bullet/m_{\rm bin})^{1/3}$ gives

$$v_{\rm ej}\approx2^{3/4}\beta^{1/4}\left(\frac{Gm_c}{a}\right)^{1/2}\left(\frac{M_\bullet}{m_{\rm bin}}\right)^{1/6}.$$

> **Worked example.** A star ejected from a BH–star binary ($m_c = 10\thinspace M_\odot$, $m_{\rm bin} = 11\thinspace M_\odot$,
> $a = 0.1$ AU) by a $10^6\thinspace M_\odot$ black hole at $\beta = 1$: $(Gm_c/a)^{1/2} = 298$ km/s and $(9\times10^4)^{1/6} = 6.7$, so
> $v_{\rm ej}\approx3400$ km/s. For a star–star binary ($m_c = 1$, $m_{\rm bin} = 2$) at the same $a$ and $\beta$ the speed is
> $\sqrt{10}\thinspace(2/11)^{1/6}\approx2.4$ times lower. That's faster than the Milky Way's escape speed (~550 km/s near the Sun). The Hills mechanism
> is the leading explanation for observed **hypervelocity stars** (Brown 2015).

**This same energy kick sets the second flare in a μTDE** (§12): debris inherit roughly the same spread of energies relative
to the massive black hole.

### 9.4 Chaos and why we need statistics

Three-body motion is generally **chaotic**: two encounters starting with almost identical conditions can end completely
differently, because small differences grow exponentially during close approaches. There is no formula giving the
outcome of one encounter. Instead we ask *statistical* questions:
- What is the probability of a μTDE for binaries of this type?
- What is the **cross section**, the effective target area for producing it?

That requires thousands of simulated encounters per parameter set (§18), and modern statistical theories of the
three-body problem (e.g. Stone & Leigh 2019).

---

# Part III: Tidal disruption and debris

## 10. Tidal forces and the tidal radius

### 10.1 What a tide is

A body of size $R$ at distance $r$ from a mass $M$. Its near side is pulled by $GM/(r-R)^2$, its centre by $GM/r^2$. The
**difference** stretches the body. Taylor-expand $g(r) = GM/r^2$:

$$\Delta g\approx-R\frac{dg}{dr} = \frac{2GMR}{r^3}.$$

The tide falls off as $1/r^3$, faster than gravity itself ($1/r^2$), so it switches on abruptly at close range. A uniform
pull would just move the whole body; only the *difference* deforms it. It's the same effect that raises ocean tides on Earth.

### 10.2 The tidal radius

The body is torn apart when the tidal stretching exceeds its own surface gravity, $Gm/R^2$:

$$\frac{2GMR}{r^3}\gtrsim\frac{Gm}{R^2}\ \Rightarrow\ r_t\approx R\left(\frac{M}{m}\right)^{1/3}.$$

The factor $2^{1/3}$ is conventionally dropped, since stellar structure changes it at order unity anyway. Equivalently,
$M/r_t^3 = m/R^3$: **at the tidal radius, the black hole's "mean density" inside the orbit equals the body's mean density.**
The same formula gives the binary's tidal radius (§9.1), with $R\to a$ and $m\to m_{\rm bin}$.

> **Worked example: the μTDE tidal radius.** A Sun-like star near its 10 $M_\odot$ companion:
> $r_t = R_\odot\times10^{1/3} = 2.15\thinspace R_\odot = 0.0100$ AU. That's `R_TDE` in `pkg/scales.py`. (The code's `pair_event`
> uses the pair's total mass, $R_\star[(m_{\rm BH}+m_\star)/m_\star]^{1/3} = 0.0103$ AU, a 3% difference.)

---

## 11. Tidal disruption: frozen-in energies and the fallback rate

### 11.1 The frozen-in approximation

When a star is disrupted, every bit of gas keeps the velocity it had at the moment of disruption, but different bits sit
at different depths in the black hole's potential (Rees 1988; Evans & Kochanek 1989). With the star's centre on an orbit of
energy $\varepsilon_\star$, a gas element offset by $\Delta r$ along the line to the black hole has energy

$$\varepsilon\approx\varepsilon_\star + \frac{GM}{r_d^2}\Delta r,\qquad\Delta r\in[-R_\star, R_\star],$$

where $r_d$ is the distance at disruption. So the debris spans an energy range

$$\Delta\varepsilon\approx\frac{GMR_\star}{r_d^2}.$$

For a parabolic star ($\varepsilon_\star = 0$), the near half is bound and the far half escapes.

> **In plain words:** at the moment of disruption all parts of the star move together, but the near side sits deeper in
> the black hole's gravity well. Same speed but deeper means less energy, so the near side is bound and falls back.

### 11.2 Fallback time and the t^(−5/3) law

A bound debris element with energy $\varepsilon\lt 0$ follows a Kepler ellipse with $a = GM/2|\varepsilon|$, and returns to
pericentre after one period:

$$t(\varepsilon) = 2\pi GM\thinspace (2|\varepsilon|)^{-3/2}.$$

The **most bound** debris ($|\varepsilon| = \Delta\varepsilon$) returns first, at $t_{\min}$. Debris with smaller $|\varepsilon|$
returns later. The mass return rate follows from the chain rule:

$$\dot M = \frac{dM}{dt} = \frac{dM}{d\varepsilon}\frac{d\varepsilon}{dt},\qquad\frac{d\varepsilon}{dt} = \frac13(2\pi GM)^{2/3}t^{-5/3}.$$

If the debris is spread roughly evenly in energy ($dM/d\varepsilon\approx$ const), then $\dot M\propto t^{-5/3}$ (Phinney 1989).
Debris returning to pericentre shocks, forms an accretion disk and shines, so $\dot M(t)$ roughly traces the light curve.

> **Worked example: the μTDE fallback.** For a star of $1\thinspace M_\odot$, $1\thinspace R_\odot$ disrupted by a 10 $M_\odot$ black hole at
> $r_d = r_t$: $\Delta\varepsilon = GM R_\star/r_t^2$, and the algebra gives
> $t_{\min} = \frac{\pi}{\sqrt2}\left(\frac{M}{m_\star}\right)^{1/2}t_{\rm dyn,\star} = 2.22\times3.16\times1593\ \mathrm{s}\approx3$ hours.
> Compare a $10^6\thinspace M_\odot$ black hole: ~40 days. μTDEs are fast.

---

## 12. The μTDE inside an SMBH encounter: two flares and the delay

### 12.1 The scenario

The binary (a 10 $M_\odot$ black hole, BH1, plus a Sun-like star) dives toward the massive black hole (BH2, mass $M_\bullet$).
The massive black hole's tide distorts the binary, often driving the star onto a closer approach to BH1, until BH1 tears it
apart: a μTDE. The code stops at the instant the star crosses its tidal radius around BH1 and records the positions and
velocities of all three bodies. That snapshot is the starting point for the debris calculation.

### 12.2 Two energies, two flares

Each debris element has two specific energies:
- **relative to BH1:** $\varepsilon_1 = \tfrac12|\mathbf v-\mathbf v_{\rm BH1}|^2 - Gm_{\rm BH1}/|\mathbf r-\mathbf r_{\rm BH1}|$;
- **relative to BH2:** $\varepsilon_2 = \tfrac12|\mathbf v-\mathbf v_{\rm BH2}|^2 - GM_\bullet/|\mathbf r-\mathbf r_{\rm BH2}|$.

The debris splits into three populations (`classify` in `debris_tracer.py`), assigned **hierarchically**:
1. $\varepsilon_1\lt 0$: **bound to BH1.** This debris falls back onto the 10 $M_\odot$ black hole within hours: the **first flare**.
2. otherwise, $\varepsilon_2\lt 0$: **escaped BH1 but bound to the massive black hole.** This debris travels out on wide orbits
   around BH2 and returns months later: a **second flare**.
3. otherwise: **unbound** from both.

> **Why hierarchical and not "most bound"?** BH2's potential well is vastly deeper than BH1's. Near BH1 essentially every
> particle has $\varepsilon_2\lt 0$, simply because BH1 itself orbits BH2. Asking "which is more negative?" would wrongly
> assign everything to BH2. The physical question is nested: first, is it still held by its local disruptor? Only if not,
> is it held by the big one?

### 12.3 The delay between the flares

How long after the μTDE does BH2-bound debris come back? Its energy spread relative to BH2 comes from the Hills-type kick
(§9.3): $\Delta\varepsilon_2\sim v_pv_{\rm orb}$, with $v_p$ the binary's speed at pericentre and $v_{\rm orb}$ the speed within
the binary. Plug this into the Kepler return time:

$$t_{\rm peak,2}\sim2\pi GM_\bullet(2\Delta\varepsilon_2)^{-3/2}.$$

> **Derivation of the scaling.** With $r_p = r_{t,\rm bin}/\beta = a(M_\bullet/m)^{1/3}/\beta$ (writing $m = m_{\rm bin}$):
> - $v_p = \sqrt{2GM_\bullet/r_p}$;
> - $v_{\rm orb} = \sqrt{Gm/a}$;
> - so $\Delta\varepsilon_2 = v_pv_{\rm orb} = \frac{G}{a}\sqrt{2\beta}\thinspace m^{2/3}M_\bullet^{1/3}$.
>
> Substituting:
>
> $$t_{\rm peak,2}\sim\frac{2\pi GM_\bullet}{\left[2\frac Ga\sqrt{2\beta}\thinspace m^{2/3}M_\bullet^{1/3}\right]^{3/2}} = \frac{P_{\rm bin}}{(8\beta)^{3/4}}\left(\frac{M_\bullet}{m}\right)^{1/2},\qquad P_{\rm bin} = 2\pi\sqrt{\frac{a^3}{Gm}}.$$
>
> This is exactly the estimate `Dt_pred` in `mutde_lab.py`.

> **Worked example.** $M_\bullet = 10^6\thinspace M_\odot$, $a = 0.1$ AU ($P_{\rm bin} = 3.5$ d), $\beta = 1.5$:
> $t_{\rm peak,2}\approx3.5\ \mathrm d\times\sqrt{10^6/11}/(12)^{3/4} = 3.5\times301/6.45\approx160$ days.
> So the signature is a **double flare**: a short μTDE flare lasting hours to days, followed ~5 months later by a second
> flare when debris falls onto the massive black hole.

**How the code measures it.** `reduce_delay.py` doesn't integrate orbits. For each μTDE it seeds 500 polytrope debris
particles at the disruption state, computes $\varepsilon_1$ and $\varepsilon_2$, and takes the **90th percentile** of
$|\varepsilon|$ among the bound particles as the "most-bound edge". It then converts each edge into a peak time with
$t = 2\pi M(2|\varepsilon|)^{-3/2}$, and defines

$$\Delta t_{\rm delay} = t_{\rm peak,BH2} - t_{\rm peak,BH1}.$$

The 90th percentile (rather than the single most bound particle) makes the estimate robust against random outliers.
`survey_delay.ipynb` then fits how $\Delta t_{\rm delay}$ scales with the parameters, using only events with pericentre
$r_p\geq100\thinspace r_g$ (closer passages need general relativity). The black-hole mass and the binary separation carry
almost all the information; with $x_M = \log_{10}(M_\bullet/M_\odot)$ and $x_a = \log_{10}(a/{\rm AU})$,

$$\log_{10}\frac{\Delta t_{\rm delay}}{\rm yr} = 0.2511-0.1181\thinspace x_M+2.0696\thinspace x_a+0.0353\thinspace x_M^2+0.1251\thinspace x_a^2-0.1473\thinspace x_Mx_a ,$$

with a scatter of 0.68 dex between individual events ($R^2 = 0.80$). Adding $e$ reduces the residual variance by only 3.7%,
and $v_\infty$ has no effect. For $M_\bullet = 10^6\thinspace M_\odot$ and $a = 0.1$ AU the fit gives about 200 days.

---

## 13. Ballistic debris physics: how the debris code works

### 13.1 The model

After disruption, the star is replaced by $N$ **massless test particles** ("debris"). Each feels the gravity of the two
black holes but not of each other, and doesn't pull on the black holes. This is the **ballistic approximation**: it ignores
gas pressure, shocks and the debris' self-gravity.

Why is this reasonable?
- Once the star is torn apart, the debris expands and its internal pressure becomes negligible compared with the black
  holes' gravity. Its orbits are then essentially free-fall orbits.
- The return times are set by the energies frozen in at disruption (§11), which ballistic motion captures exactly.

Neglecting the debris' pull on BH1 is about a 10% approximation, since the debris carries ~1 $M_\odot$ compared with BH1's
10 $M_\odot$. Hydrodynamics (shocks, disk formation) would be needed for the light curve itself.

### 13.2 Initial conditions: position offsets plus spin

Each debris particle $k$ gets (`make_debris_ic`):

$$\mathbf r_k = \mathbf r_\star + \delta\mathbf r_k,\qquad\mathbf v_k = \mathbf v_\star + \boldsymbol\omega\times\delta\mathbf r_k,$$

- $\delta\mathbf r_k$ is drawn from the $n=3$ polytrope (§7.2);
- $\boldsymbol\omega$ is the star's spin. By default ("tidal" mode) it equals the star's **orbital angular velocity about BH1 at
  pericentre**, $h/r_{\rm peri}^2$, along the orbital angular momentum. Tidal torques act most strongly at pericentre
  ($\propto1/r^3$), so the star is spun up to roughly co-rotation there.

**Both ingredients matter.** Spin with every particle at the star's centre would give no energy spread at all. Position
offsets alone give the "geometric" spread $GMR_\star/r^2$, and spin adds to it.

> **Worked example: is the spin physical?** Break-up spin is $\Omega_{\rm break} = \sqrt{Gm_\star/R_\star^3}$, about
> $6\times10^{-4}$ rad/s for the Sun (an equatorial speed of 440 km/s). The pericentre angular velocity about BH1 grows as
> $\beta^{3/2}$, and exceeds break-up for penetrations $\beta\gtrsim0.8$. Deeply plunging stars are spun past break-up as they
> are torn apart, which is consistent with their being disrupted.

### 13.3 Integration

The two black holes follow their own two-body orbit, unaffected by the massless debris, so all debris can be advanced
together in one big set of differential equations.

- **`debris_tracer.py`** packs every position and velocity into one vector and calls `scipy.integrate.solve_ivp` with the
  **DOP853** method (§14.4), tolerances `rtol = 1e-9`, `atol = 1e-12`. The force on the debris is **softened**:
  $1/r^2\to r/(r^2+s^2)^{3/2}$ with $s\approx R_\star$. That prevents infinite accelerations for particles passing very close to
  a black hole (such material would be swallowed anyway). The accuracy check is the black holes' two-body energy
  conservation (`dE_bh`).
- **`mutde_lab.py`** uses **velocity Verlet** (§14.3) in **two phases**, because the problem has two very different
  timescales:
  - **Phase A:** small time steps (1/240 of the mini-disk orbital period at $r_t$) to resolve the disruption and the
    formation of the "mini-disk" around BH1 (orbital period ~ hours).
  - **Freeze:** BH1-bound debris within $5r_t$ is then **frozen**: it rides rigidly with BH1, since nothing new happens to it.
  - **Phase B:** large time steps to follow the BH2-bound debris on its months-long return.

  This is a standard trick for **stiff** problems (timescales differing by ~1000 times). The measured first-return times
  give the numerical $\dot M(t)$, which is compared with the frozen-in estimate.

---

# Part IV: Numerical methods

## 14. Integrating orbits: from Euler to symplectic methods

### 14.1 The problem

Newton's law is a set of ordinary differential equations (ODEs):

$$\frac{d\mathbf r_i}{dt} = \mathbf v_i,\qquad\frac{d\mathbf v_i}{dt} = \mathbf a_i(\mathbf r_1,\dots,\mathbf r_N) = -\sum_{j\ne i}\frac{Gm_j(\mathbf r_i-\mathbf r_j)}{|\mathbf r_i-\mathbf r_j|^3}.$$

Given positions and velocities now, we want them later. Except for two bodies, there is no formula, so we step forward in
small time steps $\Delta t$.

### 14.2 Euler's method and why it fails

The simplest scheme: $\mathbf r_{n+1} = \mathbf r_n+\mathbf v_n\Delta t$, $\mathbf v_{n+1} = \mathbf v_n+\mathbf a_n\Delta t$. Its error per step is
$O(\Delta t^2)$, so the accumulated error is $O(\Delta t)$: "first order". Worse, for orbits it **steadily gains energy**: a
planet integrated with Euler's method spirals outward. You can see this geometrically. The step follows the tangent line,
which always points slightly outside a circle.

### 14.3 Leapfrog (velocity Verlet): a symplectic method

**Kick–drift–kick:**

$$\mathbf v_{n+1/2} = \mathbf v_n+\tfrac12\mathbf a(\mathbf r_n)\Delta t,\qquad\mathbf r_{n+1} = \mathbf r_n+\mathbf v_{n+1/2}\Delta t,\qquad\mathbf v_{n+1} = \mathbf v_{n+1/2}+\tfrac12\mathbf a(\mathbf r_{n+1})\Delta t .$$

- It is **second order** (error $O(\Delta t^2)$).
- It is **time-reversible**: run it backwards and you return exactly to the start.
- It is **symplectic**: it exactly conserves a slightly perturbed version of the true energy. So the energy error
  *oscillates* but doesn't drift, even over millions of orbits.

This is why leapfrog is the workhorse of gravitational dynamics, and why it forms the core of SpaceHub's integrator (§16).

> **In plain words:** Euler's method takes a step using only information at the start, so it always overshoots in the
> same direction. Leapfrog splits each step symmetrically (half a kick, a drift, half a kick), so errors made in the first
> half are undone in the second half.

### 14.4 Adaptive Runge–Kutta (DOP853)

Runge–Kutta methods evaluate the forces at several points within a step and combine them cleverly to cancel error terms.
**DOP853** (Dormand & Prince; Hairer, Nørsett & Wanner 1993) is an 8th-order method. It also computes an embedded
lower-order estimate, and the difference between the two estimates the error of each step. If the error exceeds the
tolerance (`rtol`, `atol`), the step is retried with a smaller $\Delta t$. Otherwise $\Delta t$ grows. The debris tracer uses it (the Lane–Emden solver uses SciPy's default
fifth-order RK45 with a tight tolerance). It's very accurate per step but not symplectic, so energy slowly drifts. That's acceptable
for the short debris integrations, not for long few-body runs.

### 14.5 Tolerances and error control

A **relative tolerance** `rtol` asks each step to be accurate to that fraction of the solution's size. The main
simulations use `rtol = 1e-13`, close to the limit of double-precision arithmetic (~$10^{-16}$). The energy error, the
fractional change of total energy $|E(t)-E(0)|/|E(0)|$, is the standard health check; `EnergyErrWriter` records it in
`--mode full`.

---

## 15. Why few-body problems are hard: close encounters and regularization

### 15.1 The singularity

The force $\propto1/r^2$ blows up as $r\to0$. During a close approach:
- velocities become huge and change direction within a short time, of order $\sqrt{r^3/GM}$;
- a fixed time step that's fine far away is catastrophically large near pericentre;
- even adaptive steps shrink so much that **round-off errors** accumulate.

In our problem, the separations range from $10^{-2}$ AU (star near BH1) to $10^5$–$10^6$ AU (start of the encounter), and the
timescales from hours to centuries. That's a spread of $10^7$ in both.

### 15.2 Time transformation and regularization

The key idea: **step in a different "time" variable** $s$, related to physical time by $dt = g(\mathbf r)\thinspace ds$, where $g$ is
small when bodies are close. Equal steps in $s$ then automatically become tiny steps in $t$ near pericentre.
- **Kustaanheimo–Stiefel (KS) regularization** (1965) goes further: a change of coordinates turns the singular Kepler problem
  into a harmonic oscillator, which has no singularity at all.
- **Algorithmic regularization (AR)** (Mikkola & Tanikawa 1999; Preto & Tremaine 1999) achieves a similar effect without
  changing coordinates. It uses a time transformation chosen so that **leapfrog with a constant step in $s$ integrates
  two-body collisions exactly** (apart from a phase error). A common choice is the *logarithmic Hamiltonian*: the step is
  controlled by the potential energy $U$, with $dt\propto ds/U$.
- **Chain coordinates** (Mikkola & Aarseth 1993) describe the system by the vectors between neighbouring bodies along a
  "chain" rather than by absolute positions. That avoids subtracting two large, nearly equal numbers (losing precision)
  when two bodies are close.

### 15.3 Extrapolation: Bulirsch–Stoer

Leapfrog's error is a series in **even** powers of the step, $\Delta t^2, \Delta t^4,\dots$. **Gragg–Bulirsch–Stoer (GBS)
extrapolation** exploits this:
1. Integrate the same interval with 2, 4, 6, … substeps.
2. Fit the results as a function of the step size squared.
3. Extrapolate to **zero step size**.

Each extra level removes another error term. This turns a second-order method into a very high-order one, with automatic
error control: stop when successive extrapolations agree to within `rtol` (Bulirsch & Stoer 1966; Press et al., *Numerical Recipes*).

**AR + chain + GBS** is the standard recipe for accurate few-body integrations with close encounters (Mikkola & Merritt
2008). It's what makes SpaceHub suitable for this project.

### 15.4 Velocity-dependent forces (PN terms)

Leapfrog's "drift then kick" structure assumes the force depends only on positions. PN accelerations depend on velocities
too, which breaks the split. Known fixes:
- an **auxiliary velocity** that is evolved alongside the real one (Hellström & Mikkola 2010);
- an **extended phase space** with duplicated coordinates (Pihajoki 2015).

---

## 16. SpaceHub: the integrator used in this project

SpaceHub (Wang, Leigh, Liu & Perna 2021, MNRAS 505, 1053) is a C++ library of few-body integrators. The project uses
commit `190777d6`. Everything below comes from that version of the source. File names refer to SpaceHub's `src/` directory.

### 16.1 What the driver asks for

```cpp
using Forces = Interactions<NewtonianGrav, PN1, PN2, PN2p5>;
using Solver = methods::DefaultMethod<Forces, particles::SizeParticles>;
args.rtol = 1e-13;
```

`DefaultMethod` is defined as **`AR_Chain_Plus`** (`spaceHub.hpp`), which is:
- **AR-chain:** chain coordinates plus algorithmic regularization (§15.2);
- base step: **drift–kick–drift leapfrog** (§14.3);
- **Gragg–Bulirsch–Stoer extrapolation** (§15.3);
- an error estimate based on the "worst offender" variable, and a step-size controller;
- numbers stored as **Kahan-compensated doubles**. Kahan summation keeps a running correction for the round-off lost in
  each addition, which matters over billions of small increments.

### 16.2 Chain coordinates

The bodies are ordered into a chain by sorting the pairwise distances, so that the closest pairs are neighbours. The state
is stored as the vectors between consecutive bodies, $\mathbf X_k = \mathbf r_{k+1}-\mathbf r_k$, plus the centre of mass. The
chain is rebuilt after every accepted step, so it follows the system as different pairs become close. For our three
bodies, when the star is near BH1, the BH1–star separation is stored *directly*. It isn't computed as the difference of
two large positions measured from the SMBH, which would lose digits.

> **Worked example: why that matters.** At the moment of a μTDE, BH1 and the star are ~0.01 AU apart, while both are ~3 AU
> from the SMBH. With positions stored relative to the centre of mass, their separation is a difference of two numbers
> near 3, accurate to about $3\times10^{-16}$ in absolute terms. That's a relative accuracy of $3\times10^{-14}$ on the
> separation, close to `rtol` itself. Storing the separation directly keeps the full $10^{-16}$ relative precision.

### 16.3 The logarithmic-Hamiltonian time transformation

SpaceHub's default regularization is the **logarithmic Hamiltonian (LogH)** scheme (`ReguType::LogH` in `archain.hpp`).
Let $T$ be the kinetic energy and $U\lt 0$ the Newtonian potential energy, and let $B$ be minus the total energy (stored at the
start and updated as external forces such as PN do work). The two leapfrog sub-steps advance physical time differently:

$$\text{drift:}\quad dt = \frac{ds}{T+B},\qquad\text{kick:}\quad dt = \frac{ds}{-U}.$$

Both denominators grow like $1/r$ during a close approach (since $T+B = T - E = -U$ along the exact solution). So a fixed step
$ds$ automatically becomes a physical step $dt\propto r$, tiny near pericentre and large far away. The remarkable property of
this transformation (Mikkola & Tanikawa 1999; Preto & Tremaine 1999) is that leapfrog then follows a two-body Kepler orbit
*exactly*, even through a head-on collision, apart from a timing (phase) error. The extrapolation then removes that phase error.

> **In plain words:** the integrator automatically slows down time where things happen fast. Near a close approach it takes
> many small steps; during the long, slow parts of the journey it takes a few big ones.

### 16.4 Velocity-dependent forces: the auxiliary velocity

The PN accelerations depend on velocities (§6). SpaceHub keeps a second, **auxiliary velocity** $\mathbf w$ alongside the real
velocity $\mathbf v$ (Hellström & Mikkola 2010). Each kick is split symmetrically (`archain.hpp`):
1. half-kick the real velocity $\mathbf v$, using PN forces evaluated with $\mathbf w$;
2. full kick of $\mathbf w$, using PN forces evaluated with $\mathbf v$;
3. half-kick $\mathbf v$ again, using $\mathbf w$.

This keeps the step time-symmetric, which the extrapolation needs. The energy variable $B$ is updated by the work done by
the PN forces, $dB/dt = -\sum_i m_i\mathbf v_i\cdot\mathbf a_{{\rm PN},i}$. In this way the 2.5PN radiation-reaction energy loss is tracked.

### 16.5 The Bulirsch–Stoer driver and how rtol is used

- Each big step $H$ (in $s$) is integrated with $n$ leapfrog substeps, $n\in\lbrace 1, 2, 3, 5, 8, 12, 17, 25\rbrace$. The results are
  extrapolated to $n\to\infty$ with an Aitken–Neville polynomial in $h^2$, up to 8 levels (`Bulirsch-Stoer.hpp`).
- **Error estimate** (`worst-offender.hpp`): compare the last two extrapolations, variable by variable, and take the worst
  (the formula follows this list). The step is accepted if $\mathrm{err}\le1$. The driver sets $\mathrm{rtol} = 10^{-13}$
  (SpaceHub's own default is $10^{-14}$) and $\mathrm{atol} = 0$.
- **Next step size:** $H_{\rm new} = H\cdot0.72\thinspace (0.95/\mathrm{err})^{1/(2k+1)}$ for extrapolation order $2k+1$, limited so it
  can't change too abruptly. The extrapolation order itself adapts, to minimize work per unit of $s$.
- **First step:** a tenth of the shortest free-fall time in the system.

The error estimate used above:

$$\mathrm{err} = \max_i\frac{\left|y_i^{(k)}-y_i^{(k-1)}\right|}{\mathrm{atol}+\mathrm{rtol}\cdot\max(\left|y_i^{(k)}\right|,\left|y_i^{(0)}\right|)} .$$

Here $y_i^{(k)}$ and $y_i^{(k-1)}$ are the last two extrapolated values of variable $i$, and $y_i^{(0)}$ is its value at the start of the step.

### 16.6 The forces

- **Newtonian:** $\mathbf a_i = -\sum_{j\ne i}Gm_j(\mathbf r_i-\mathbf r_j)/r_{ij}^3$.
- **PN1, PN2, PN2.5** (`interaction/post-newtonian.hpp`): the **two-body** equations of motion in harmonic coordinates from
  Blanchet (2014, Living Reviews, eq. 203; cited in the code), applied to **every pair** and summed.

For example, the 1PN acceleration of body 1 due to body 2, with $\mathbf n = (\mathbf r_1-\mathbf r_2)/r$ and $\mathbf v_{12} = \mathbf v_1-\mathbf v_2$, is

$$\mathbf a_1^{\rm 1PN} = \frac{Gm_2}{r^2c^2}\left\lbrace\left[-v_1^2-2v_2^2+4\thinspace\mathbf v_1\cdot\mathbf v_2+\tfrac32(\mathbf n\cdot\mathbf v_2)^2+\frac{5Gm_1}{r}+\frac{4Gm_2}{r}\right]\mathbf n+\left[4\thinspace\mathbf n\cdot\mathbf v_1-3\thinspace\mathbf n\cdot\mathbf v_2\right]\mathbf v_{12}\right\rbrace.$$

The 2.5PN (radiation-reaction) term has prefactor $\frac{8}{5}\frac{G^2m_1m_2}{r^3c^5}$, and is what removes orbital energy as gravitational waves.

**Caveat:** summing two-body PN terms pair by pair omits the genuinely three-body PN cross terms of the full
Einstein–Infeld–Hoffmann equations, such as the SMBH's potential modifying the BH1–star interaction. Those terms are of
relative size $\sim r_g/r_p$ ($\sim10^{-3}$ for the example in §6.2). This is a standard approximation in few-body PN codes.

### 16.7 Stellar radii

When no radius is supplied, `stellar::stellar_radius` uses:
- stars: $R = R_\odot(m/M_\odot)^{0.75}$;
- white dwarfs: $R = 0.01\thinspace R_\odot(M_\odot/m)^{1/3}$;
- neutron stars: 10 km;
- black holes: $2Gm/c^2$.

The BH–star campaign passes the radii explicitly ($1\thinspace R_\odot$ for the star; $2r_g$ for BH1).

### 16.8 The energy-error diagnostic

For the regularized integrator, `EnergyErrWriter` records $\log\left|(T+B)/U\right|$, with Newtonian $T$ and $U$. Along the
exact solution, $T+B = -U$, so this quantity is zero. Its deviation from zero measures how well the LogH Hamiltonian is
conserved, and the PN work enters through $B$. It is normalized by the potential energy, not by the initial energy.

### 16.9 Initial-condition helpers

- **`orbit::Hyperbolic(M, m, v_inf, b, ...)`** builds the incoming hyperbola from $v_\infty$ and $b$:
  $\mu = G(M+m)$, $a = -\mu/v_\infty^2$, $e = \sqrt{1+b^2/a^2}$, $p = (bv_\infty)^2/\mu$. It places the body at distance $r$ on the
  **incoming** branch, with true anomaly $\nu = -\arccos[(p-r)/(er)]$, negative because the body is approaching.
- **`orbit::calc_a_e`** recovers $a = -\mu/(v^2-2\mu/r)$ and the eccentricity from the Laplace–Runge–Lenz vector,
  $e = \left|\mathbf r(v^2-\mu/r)-\mathbf v(\mathbf r\cdot\mathbf v)\right|/\mu$.

---

## 17. The driver `smbh_binary.cpp`: setup, stopping conditions, classification

### 17.1 Initial conditions (element mode)

For one encounter, with all quantities converted to code units:
1. **Inner binary.** BH1 and the star are placed on a Kepler orbit with elements $(a,e,i,\Omega,\omega,\nu)$
   (`orbit::Elliptic`) and moved to their own centre-of-mass frame. The inclination is measured relative to the outer orbit's plane.
2. **Outer orbit.** The binary's pericentre is $r_p = r_{t,\rm bin}/\beta$, and the impact parameter follows from
   $b^2 = r_p^2+2GM_{\rm tot}r_p/v_\infty^2$ (§4). The binary is placed on the incoming branch of the hyperbola
   (`orbit::Hyperbolic(..., Hyper::in)`) at distance $r_{\rm start} = \max(50\thinspace r_{t,\rm bin},\ 10\thinspace r_p)$.
   **Why there:** the tidal perturbation scales as $(r_{t,\rm bin}/r)^3$, which is about $10^{-5}$ at $50\thinspace r_{t,\rm bin}$, so the
   binary starts essentially unperturbed. Starting farther out (as far as $b$) would waste thousands of binary orbits on an
   uneventful approach. The encounter geometry is fixed by $(r_p, b)$, and the binary's phase is random, so starting closer
   introduces no bias.
3. **Everything is moved to the centre-of-mass frame of all three bodies.**

### 17.2 Stopping conditions

The integration ends at the first of:
1. **Terminal event** (checked every step, for every pair, by `pair_event`):
   - **TDE:** if the tidal radius of the extended body is the outermost scale (outside both the contact radius
     $R_i+R_j$ and $4r_g$), a TDE fires when the separation drops below $r_t$.
   - **COLLISION:** otherwise, if the bodies' surfaces are the outer scale, a collision fires on contact.
   - **PLUNGE:** for a compact pair (both black holes, or above the Hills mass), the osculating pericentre is tested
     against $4r_g$ (only within $50r_g$, to save time).
2. **Escape:** a body is farther than $r_{\rm far} = 2r_{\rm start}$ from the centre of mass of the other two, *receding*, **and**
   energetically unbound from them ($\varepsilon\gt 0$). The energy clause prevents stopping a body that is merely at the
   apocentre of a wide but bound orbit; such a body will come back for more.
3. **Time backstop:** $t \gt n_{\rm dyn}t_{\rm dyn}(r_{\rm start})$, with $n_{\rm dyn} = 2\times10^4$ in the BH–star campaign. Runs that
   end this way are flagged `exit_mode = backstop`.

### 17.3 Classification

If no terminal event occurred, the final energies decide the outcome (with $\varepsilon$ from §2):
- inner pair bound ($\varepsilon_{12}\lt 0$) and bound to the SMBH → BOUND_TRIPLE; inner pair bound but free of the SMBH → FLYBY;
- inner pair unbound, one member bound to the SMBH → EXCHANGE; both members bound to it → BOUND_TRIPLE; neither → IONIZATION.

The code also records the surviving orbit's $(a,e)$, the escaper's speed, the closest approaches of each pair (in $r_g$),
the achieved penetration $\beta_{\rm actual}$, and the **final positions and velocities of all three bodies**. For a μTDE,
those final values are the disruption-time state that seeds the debris calculations (§12–13), and they can restart a
follow-up run (`--restart`).

---

## 18. Monte Carlo scattering experiments and cross sections

### 18.1 The idea

We want the **cross section** $\Sigma_X$ for outcome $X$: the effective area of the target. If binaries arrive with
number flux $n v_\infty$ (per area per time), outcome $X$ happens at rate $n\Sigma_Xv_\infty$. We measure it by throwing
binaries at random impact parameters within a large disk of radius $b_{\max}$:

$$\Sigma_X = f_X\thinspace\pi b_{\max}^2,\qquad f_X = \frac{N_X}{N}.$$

### 18.2 Sampling the impact parameter

Random arrival points are uniform **in area**, not in radius. The area between $b$ and $b+db$ is $2\pi b\thinspace db$, so
$p(b)\propto b$, i.e. **$b^2$ is uniform** on $[0, b_{\max}^2]$. The sweep draws $b^2$ uniformly, converts to $r_p$ by inverting
$b^2 = r_p^2+r_{\rm foc}r_p$:

$$r_p = \frac{-r_{\rm foc}+\sqrt{r_{\rm foc}^2+4b^2}}{2},$$

and passes $\beta = r_{t,\rm bin}/r_p$ to the code.

### 18.3 Choosing b_max

$b_{\max}$ must be large enough to contain every encounter that can produce outcome $X$, but not so large that most runs
are wasted on distant flybys. The campaign sets the largest pericentre

$$r_{p,\max} = \max(10\thinspace r_{t,\rm bin},\ 8\thinspace r_g),\qquad b_{\max}^2 = r_{p,\max}^2 + r_{\rm foc}\thinspace r_{p,\max}$$

(`pkg/scales.py`; the first campaign used $5\thinspace r_{t,\rm bin}$). Beyond $10\thinspace r_{t,\rm bin}$ the tide is $10^{-3}$ of the binary's
self-gravity, too weak to drive a μTDE. Doubling the cut from 5 to 10 and checking that the μTDE count converged is a
**convergence test** of the kind recommended in Part 0.

**Note:** the *fraction* $f_X$ depends on the arbitrary $b_{\max}$, but the *cross section* $\Sigma_X = f_X\pi b_{\max}^2$
doesn't, as long as $b_{\max}$ encloses all events.

### 18.4 Random orientations

The binary's orientation is random:
- $\cos i$ uniform in [−1, 1] (uniform directions on a sphere; choosing $i$ itself uniformly would over-represent the poles);
- $\Omega$ and $\omega$ uniform in [0, 2π];
- phase uniform in mean anomaly (§3.3).

### 18.5 Statistical uncertainty

Each run either is or isn't a μTDE: a **binomial** trial. With $N$ runs and $N_X$ successes, $\hat f = N_X/N$ has standard error

$$\delta f = \sqrt{\frac{\hat f(1-\hat f)}{N}}.$$

> **Worked example.** With $N = 1500$ runs per cell and $\hat f = 0.04$: $\delta f = \sqrt{0.04\times0.96/1500} = 0.005$, a 13%
> relative error per cell. Pooling over many cells reduces it further; the notebooks quote per-mass errors below 1%.

### 18.6 The campaign grid

| Parameter | Values | Count |
|---|---|---|
| $M_\bullet$ | 10 to 10¹⁰ $M_\odot$, half-decade steps | 19 |
| $a$ | 0.03 to 30 AU | 13 |
| $v_\infty/v_{\rm orb}$ | 0.001 to 10 | 7 |
| $e$ | 0 to 0.99 (survival cut applied) | ≤ 12 |
| runs per cell | | 1500 |

That's about $2.8\times10^7$ simulated encounters, stored in the HDF5 "bundle" (archived on Zenodo).

---

# Part V: From simulations to the Universe

## 19. Event rates: loss-cone theory

### 19.1 From cross section to rate

A single black hole embedded in a population of binaries with number density $n$ and velocity dispersion $\sigma_v$ sees
outcome $X$ at rate $\Gamma\approx n\thinspace\Sigma_X\thinspace v$. In a galactic nucleus, though, the supply of binaries on orbits that
reach the black hole is limited, and that's what loss-cone theory describes.

### 19.2 The loss cone

A binary on an orbit around the massive black hole only gets close enough if its angular momentum is small. The orbits
with $r_p$ inside the critical radius form a narrow cone of velocity directions: the **loss cone**. Its angular size at
radius $r$ is $\theta_{\rm lc}^2\approx\Sigma/\pi r^2$ (Syer & Ulmer 1999 form, as used in `mutde_losscone.ipynb`). Those orbits are
emptied within one orbital period. They're refilled by **two-body relaxation**: random gravitational nudges between stars,
which slowly change orbits on the relaxation time

$$t_r = \frac{0.34\thinspace\sigma_v^3}{G^2m_\star\rho\ln\Lambda}$$

(Spitzer 1987; Binney & Tremaine 2008). Two regimes (Frank & Rees 1976; Lightman & Shapiro 1977):
- **Empty loss cone** (close to the black hole): relaxation is slow and the loss cone is refilled only diffusively.
- **Full loss cone** (far out): relaxation refills the cone faster than it empties.

The rate per galaxy is dominated by the radius $r_{\rm crit}$ where the two regimes meet, and roughly

$$\Gamma\sim\frac{N_{\rm bin}(\lt r_{\rm crit})}{t_r(r_{\rm crit})}.$$

The notebook interpolates between the regimes and integrates over a two-power-law nuclear star cluster.

### 19.3 Rate density in the Universe

Multiplying the per-galaxy rate by the number of galaxies hosting black holes of each mass (the black-hole mass function),
times the occupation fraction, gives a volumetric rate in Gpc⁻³ yr⁻¹. Integrating over cosmic volume gives the
number of double flares per year out to a given redshift. The ordinary star–SMBH TDE rate is computed in the same framework
for comparison, and its ratio to the double-flare rate scales with the binary fraction.

---

## 20. The f_μ fit and the population average

Running 28 million encounters gives $f_\mu$ on a grid. To use it in rate calculations, and to understand it, the analysis
fits a **semi-analytic** form (`analyze_mutde_cross_sections.ipynb`). The physical idea: a μTDE needs the encounter to drive
the star inside BH1's tidal radius. That requires a penetration above some critical value $\beta_{\rm crit}$, which depends on
how close the binary already is to self-disruption:

$$\rho\equiv\frac{a(1-e)}{r_{\rm TDE}}\qquad(\text{binary pericentre in units of the }\mu\text{TDE tidal radius}),$$

and on $e$. A binary with $\rho$ just above 1 needs only a small nudge; a wide binary needs a deep plunge. The fraction is
then the fraction of the $b^2$-disk inside the critical pericentre $r_c = r_{t,\rm bin}/\beta_{\rm crit}$:

$$f_\mu = \frac{r_c^2+r_{\rm foc}r_c}{r_{p,\max}^2+r_{\rm foc}r_{p,\max}},\qquad\log_{10}\beta_{\rm crit} = c_0+c_1L+c_2L^2+c_3e+c_4e^2,\ \ L = \log_{10}\rho .$$

In the strong-focusing limit this reduces to $f_\mu\approx1/(10\beta_{\rm crit})$, a function of $(a,e)$ alone, so $M_\bullet$ and
$v_\infty$ enter only through $b_{\max}^2$ in $\Sigma_\mu = f_\mu\pi b_{\max}^2$. The fitted coefficients live in `fmu_fit.json`. Here $r_{\rm TDE} = R_\star(m_{\rm BH1}/m_\star)^{1/3} = 0.0100$ AU.
The fit reproduces the simulated $f_\mu$ with a median error of about 0.1 dex ($R^2 = 0.91$ on $\log_{10}f_\mu$).

**Population averages** then weight over assumed distributions $p(a,e)$ of binaries in galactic nuclei (three regimes from
Stone & Leigh), giving $\langle\Sigma_\mu(M_\bullet,v)\rangle$ for the rate calculation.

---

## Appendix A: Symbols

| Symbol | Meaning | Typical units |
|---|---|---|
| $M_\bullet$ | massive (central) black hole mass, "BH2" | $M_\odot$ |
| $m_1$, $m_{\rm BH1}$ | stellar-mass black hole in the binary | $M_\odot$ |
| $m_2$, $m_\star$, $R_\star$ | star mass, radius | $M_\odot$, $R_\odot$ |
| $m_{\rm bin}$ | binary mass $m_1+m_2$ | $M_\odot$ |
| $a$, $e$ | binary semi-major axis, eccentricity | AU, — |
| $i,\Omega,\omega,\nu$ | inclination, node, argument of pericentre, true anomaly | rad |
| $\mathcal M$, $E$ | mean, eccentric anomaly | rad |
| $v_\infty$ | approach speed at infinity | km/s |
| $v_{\rm orb}$ | binary orbital speed $\sqrt{Gm_{\rm bin}/a}$ | km/s |
| $b$, $b_{\max}$ | impact parameter, sampling cut | AU |
| $r_p$ | pericentre of the outer orbit | AU |
| $r_{\rm foc}$ | focusing length $2GM/v_\infty^2$ | AU |
| $r_{t,\rm bin}$ | binary tidal radius $a(M_\bullet/m_{\rm bin})^{1/3}$ | AU |
| $\beta$ | penetration factor $r_{t,\rm bin}/r_p$ | — |
| $r_t$, $r_{\rm TDE}$ | stellar tidal radius $R_\star(m_{\rm BH1}/m_\star)^{1/3}$ | AU |
| $r_g$, $r_s$ | gravitational radius $GM/c^2$, Schwarzschild radius $2r_g$ | AU, km |
| $\varepsilon$ | specific orbital energy | (km/s)², erg/g |
| $h$ | specific angular momentum | AU km/s, cm²/s |
| $\Sigma_X$, $f_X$ | cross section, fraction of outcome $X$ | AU², pc² |
| $\sigma_v$ | velocity dispersion of the stars around the black hole | km/s |
| $\Delta t_{\rm delay}$ | time between the two flares | yr |
| $\Gamma$ | event rate per galaxy | yr⁻¹ |

## Appendix B: Constants and conversions

| Quantity | Value |
|---|---|
| $V_1 = \sqrt{GM_\odot/\mathrm{AU}}$ | 29.78 km/s |
| code time unit | 0.159 yr = 5.02×10⁶ s |
| $c$ in code units | 1.006×10⁴ |
| 1 AU | 215 $R_\odot$ = 1.496×10¹³ cm |
| 1 pc | 206 265 AU |
| 1 pc² | 4.25×10¹⁰ AU² |
| $r_g$ per $M_\odot$ | 1.48 km = 9.87×10⁻⁹ AU |
| 1 (km/s)² | 10¹⁰ erg/g |

## Appendix C: References

**Orbits, dynamics, and few-body physics**
- Binney, J. & Tremaine, S. (2008), *Galactic Dynamics*, 2nd ed., Princeton. Orbits, relaxation, loss cones.
- Murray, C. D. & Dermott, S. F. (1999), *Solar System Dynamics*, Cambridge. Orbital elements, Kepler's equation.
- Heggie, D. C. (1975), MNRAS 173, 729. Binary–single encounters, Heggie's law, thermal eccentricities.
- Hills, J. G. (1988), Nature 331, 687. Binary disruption by a massive black hole, hypervelocity stars.
- Brown, W. R. (2015), ARA&A 53, 15. Hypervelocity stars (review).
- Stone, N. C. & Leigh, N. W. C. (2019), Nature 576, 406. Statistical theory of the three-body problem.
- Valtonen, M. & Karttunen, H. (2006), *The Three-Body Problem*, Cambridge.

**Black holes and relativity**
- Bardeen, J. M., Press, W. H. & Teukolsky, S. A. (1972), ApJ 178, 347. Orbits around black holes (ISCO, marginally bound orbit).
- Misner, C. W., Thorne, K. S. & Wheeler, J. A. (1973), *Gravitation*, Freeman.
- Peters, P. C. (1964), Phys. Rev. 136, B1224. Gravitational-wave inspiral.
- Kidder, L. E. (1995), Phys. Rev. D 52, 821. PN two-body equations of motion.
- Blanchet, L. (2014), Living Rev. Relativ. 17, 2. Post-Newtonian theory (review).
- Hills, J. G. (1975), Nature 254, 295. The Hills mass.

**Tidal disruption**
- Rees, M. J. (1988), Nature 333, 523. Tidal disruption of stars by massive black holes.
- Evans, C. R. & Kochanek, C. S. (1989), ApJ 346, L13. Frozen-in approximation, fallback.
- Phinney, E. S. (1989), IAU Symp. 136, 543. The $t^{-5/3}$ law.
- Perets, H. B., Li, Z., Lombardi, J. C. & Milcarek, S. R. (2016), ApJ 823, 113. Micro-TDEs.
- Kremer, K. et al. (2022), ApJ 933, 203. Hydrodynamics of μTDEs.

**Loss cones and rates**
- Frank, J. & Rees, M. J. (1976), MNRAS 176, 633.
- Lightman, A. P. & Shapiro, S. L. (1977), ApJ 211, 244.
- Syer, D. & Ulmer, A. (1999), MNRAS 306, 35.
- Wang, J. & Merritt, D. (2004), ApJ 600, 149.
- Spitzer, L. (1987), *Dynamical Evolution of Globular Clusters*, Princeton.
- Merritt, D. (2013), *Dynamics and Evolution of Galactic Nuclei*, Princeton.

**Numerical methods**
- Press, W. H. et al. (2007), *Numerical Recipes*, 3rd ed., Cambridge. Runge–Kutta, Bulirsch–Stoer.
- Hairer, E., Nørsett, S. P. & Wanner, G. (1993), *Solving Ordinary Differential Equations I*, Springer. DOP853.
- Bulirsch, R. & Stoer, J. (1966), Numer. Math. 8, 1. Extrapolation.
- Kustaanheimo, P. & Stiefel, E. (1965), J. Reine Angew. Math. 218, 204. KS regularization.
- Mikkola, S. & Aarseth, S. J. (1993), Celest. Mech. 57, 439. Chain regularization.
- Mikkola, S. & Tanikawa, K. (1999), MNRAS 310, 745; Preto, M. & Tremaine, S. (1999), AJ 118, 2532. Algorithmic regularization / logarithmic Hamiltonian.
- Mikkola, S. & Merritt, D. (2008), AJ 135, 2398. AR-chain with GBS extrapolation.
- Hellström, C. & Mikkola, S. (2010), Celest. Mech. Dyn. Astron. 106, 143. Velocity-dependent forces in AR.
- Pihajoki, P. (2015), Celest. Mech. Dyn. Astron. 121, 211. Extended phase space.
- Dormand, J. R. & Prince, P. J. (1980), J. Comput. Appl. Math. 6, 19. Embedded Runge–Kutta pairs.
- Wang, Y.-H., Leigh, N. W. C., Liu, B. & Perna, R. (2021), MNRAS 505, 1053. SpaceHub, the integrator used here.

**This project**
- Kritos, K. & Perna, R. (2026), *Tidally Induced Coda: Double Transients from Black Hole–Star Binary Encounters with Massive Black Holes*, in preparation.
- Yu, F. & Lai, D. (2024), ApJ 977, 268. Binary stars approaching supermassive black holes: tidal break-up, double stellar disruptions and stellar collision (the comparison in `analysis/yulai2024/`).

## Appendix D: Changelog

- **2026-11-22.** Corrected the km/AU factor (§1), the IMBH Schwarzschild radius, the surviving $e$ nodes at $a = 1$ AU,
  and the focusing cross section in terms of $v_{\rm esc}$; the Hills ejection speed now includes $\beta$ and the companion
  mass; added the $\Delta t_{\rm delay}$ fit and the $f_\mu$ fit accuracy; cross sections written $\Sigma$ (as in the figures),
  $\sigma_v$ added to Appendix A; references to SpaceHub and Dormand & Prince added; Yu & Lai (2024) reference completed.
- **2026-11-21 (c).** Replaced `\thickspace` and `\negthinspace`, which GitHub's math renderer doesn't support (the energy
  equation in §2.1, the PN expansion in §6.1 and the 1PN force in §16.6 showed errors).
- **2026-11-21 (b).** Fixed equations that rendered incorrectly on GitHub. GitHub's Markdown escaping ate `\,` `\;` `\!` `\{` `\}`
  and double-escaped `<` `>` inside math, so these are now written `\thinspace`, `\lbrace`, `\lt`, etc. Display equations
  inside list items, or directly against text, were moved onto their own lines.
- **2026-11-21 (a).** First version.
