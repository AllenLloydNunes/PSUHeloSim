# Simplified Landing Dynamics Model — Design Plan

**Goal:** a minimal 6-DOF rigid-body landing simulator in MATLAB, built specifically to study
**ML vs. physics friction modeling** at different timesteps (physics `dt=0.001`, ML `dt=0.01`).

Standalone from the UH-60 sim — no GENHEL dependencies. Fresh code, but carrying forward the
lessons from the `LandingGear.m` frame/moment-arm bug (see §7).

---

## 1. Conventions — fix these first, write them in every file header

| choice | recommendation | why |
|---|---|---|
| Units | **SI** (m, kg, N, rad, s) | The parent sim is ft/slug/lb; mixing is how unit bugs happen. This model is about *method*, not matching UH-60 numbers. Convert only if you later cross-compare. |
| Inertial frame | **NED** — x North, y East, **z DOWN**, origin on the ground plane | Matches the parent project so intuition transfers. Ground plane is `z = 0`. **Altitude = −z.** |
| Body frame | x forward, y right, z down, origin at CG | Standard aerospace, matches parent sim |
| Euler sequence | 3-2-1 (yaw ψ → pitch θ → roll φ) | Standard |
| Integrator | **fixed-step RK4**, hand-written | `ode45` hides the timestep; this study is *about* the timestep. You need explicit control. |

> **The one gotcha with NED:** z is positive *downward*, so a foot is penetrating the ground when
> its inertial z is **positive**, and "height" plots need `-z`. Put this in a comment at the top of
> every file that touches geometry.

## 2. State vector (12 states, exactly as requested)

```
x = [ u  v  w    % body-frame velocity        (m/s)
      p  q  r    % body-frame angular rate    (rad/s)
      phi theta psi   % Euler angles          (rad)
      xN yE zD ]      % inertial position     (m)
```

## 3. Equations of motion — `xdot = f(F_b, M_b)`

Given total force `F_b` and moment `M_b` **in body frame** (§5 builds these):

**Translational** (Newton in a rotating frame):
```
udot = Fx/m - (q*w - r*v)
vdot = Fy/m - (r*u - p*w)
wdot = Fz/m - (p*v - q*u)
```
where `F_b = F_contact_b + F_gravity_b`, and gravity in body frame is:
```
F_gravity_b = m*g*[ -sin(theta);  sin(phi)*cos(theta);  cos(phi)*cos(theta) ]
```

**Rotational** (Euler's equations):
```
omegadot = I \ ( M_b - cross(omega, I*omega) )     % omega = [p;q;r]
```

**Euler kinematics** (body rates → Euler rates):
```
phidot   = p + sin(phi)*tan(theta)*q + cos(phi)*tan(theta)*r
thetadot =            cos(phi)*q -            sin(phi)*r
psidot   =   (sin(phi)/cos(theta))*q + (cos(phi)/cos(theta))*r
```

**Position kinematics**:
```
[xNdot; yEdot; zDdot] = R_b2i * [u; v; w]
```

**Rotation matrix** (body → inertial, 3-2-1):
```
R_b2i = [ cT*cS,  sP*sT*cS-cP*sS,  cP*sT*cS+sP*sS;
          cT*sS,  sP*sT*sS+cP*cS,  cP*sT*sS-sP*cS;
           -sT,        sP*cT,           cP*cT     ];   % s=sin, c=cos, P=phi,T=theta,S=psi
```

## 4. Geometry — box body + 3 legs

**Body:** rectangular box, `a` × `b` × `h` (length × width × height), mass `m` at the centroid.
```
Ixx = m*(b^2 + h^2)/12
Iyy = m*(a^2 + h^2)/12
Izz = m*(a^2 + b^2)/12
I   = diag([Ixx Iyy Izz])
```

**Legs:** two forward (splayed outward), one aft — a helicopter-like tricycle stance. For the
simple version the legs are **rigid struts with a compliant foot** (see §5), so each leg is defined
solely by its foot location in body coordinates:
```
p_foot_b(:,1) = [ +2.0; -1.5; +1.8 ];   % front-left   (outboard, below CG)
p_foot_b(:,2) = [ +2.0; +1.5; +1.8 ];   % front-right
p_foot_b(:,3) = [ -3.0;  0.0; +1.8 ];   % aft
```
The outward splay is captured by the ±1.5 m lateral offset (wider than the 2 m body). Positive z
= below the CG, per NED.

**Suggested starting parameters** (round numbers, tuned so the physics is well-conditioned):

| parameter | value | rationale |
|---|---|---|
| `m` | 1000 kg | round; scale later |
| `a, b, h` | 10, 2, 2 m | gives Iyy ≈ Izz ≈ 8667, Ixx ≈ 667 kg·m² — pitch/yaw >> roll, like a real helo |
| `k` (per leg) | 6.54e4 N/m | gives ~5 cm static deflection under gravity across 3 legs |
| `c` (per leg) | 2800 N·s/m | ζ ≈ 0.3 — visibly underdamped so you can *see* bounce/settle |
| `mu` | 0.6 | dry rubber-on-concrete-ish |
| `v_bl` | **the study variable** — see §6 | this is the knob that sets stiffness |

Identical `k` and `c` on all three legs to start, exactly as requested.

## 5. Contact and friction — the swappable core

Per leg `i`, every derivative evaluation:

```
1.  P_foot_i  = P_cg_i + R_b2i * p_foot_b(:,i)          % foot position, inertial
2.  V_foot_i  = R_b2i*( [u;v;w] + cross(omega, p_foot_b(:,i)) )   % foot velocity, inertial
3.  delta     = P_foot_i(3)                             % >0 means penetrating (NED!)
4.  if delta <= 0 : no contact, all forces zero, skip
5.  deltadot  = V_foot_i(3)
6.  N         = max( k*delta + c*deltadot , 0 )         % CLAMP: ground pushes, never pulls
7.  v_t       = V_foot_i(1:2)                           % tangential (horizontal) velocity
8.  F_t_i     = frictionFcn(v_t, N, params)             % <-- SWAPPABLE: physics or ML
9.  F_foot_i  = [ F_t_i ; -N ]                          % inertial: friction horizontal, N pushes up (-z)
10. F_foot_b  = R_b2i' * F_foot_i                       % CONVERT TO BODY *BEFORE* SUMMING
11. M_foot_b  = cross( p_foot_b(:,i), F_foot_b )        % moment arm = CONTACT POINT
```
Then `F_contact_b = sum(F_foot_b, 2)` and `M_b = sum(M_foot_b, 2)`.

**Step 6's `max(...,0)` clamp is not optional** — without it the damper term generates suction as
the leg rebounds and the aircraft gets glued to the ground.

**Physics friction** (`friction_physics.m`), regularized Coulomb:
```
speed = norm(v_t);
F_t   = -mu * N * v_t / max(speed, v_bl);
```
Magnitude saturates at `mu*N`; below `v_bl` it becomes linear viscous (this is the regularization
that makes it integrable — and the thing that costs you a small timestep).

**ML friction** (`friction_ml.m`) — identical signature, so swapping is one line:
```
function F_t = friction_ml(v_t, N, params)
```

## 6. Why `dt` matters here — the arithmetic that makes this study well-posed

> **CORRECTION (verified in step 3, `test_step3_friction_slope.m`).** The table below estimates
> the friction pole from *translation only*. That is wrong. The feet sit 1.8 m below the CG, so
> friction also resists **rotation** — and roll inertia is small, so the effective mass seen by
> lateral friction is ~171 kg rather than 1000 kg. Linearizing the actual model (5° slope, v_bl
> 0.01) gives the fastest pole at **−3,454 rad/s** (lateral/roll), not −589. Consequences:
>
> | claim below | actual |
> |---|---|
> | v_bl 0.01 stable at dt = 0.001 | **unstable** — needs dt ≤ ~0.0008, or v_bl ≥ 0.0124 |
> | v_bl 0.05 stable at dt = 0.01 | **unstable** — needs v_bl ≥ 0.130 |
>
> Minimum creep physics friction can deliver: **1.81 mm/s at dt = 0.001, 18.99 mm/s at dt = 0.01.**
> Below the stability limit the integrator doesn't blow up — friction saturates and it chatters,
> creeping at the limit value regardless of the v_bl requested. A laterally symmetric test never
> seeds the unstable mode and will falsely look fine; seed it. The original table is kept below
> for the record.
>
> **Resolution:** the default `v_bl` is now 0.1 ft/s = 0.03048 m/s, matching the UH-60 model's
> `LG_VBL`. That puts the fastest pole at −1,143 rad/s, stable at dt = 0.001 (and still unstable
> at dt = 0.01, as the experiment needs).

**The normal spring is NOT the stiff part.** Heave frequency is
`sqrt(3k/m) = sqrt(196200/1000) ≈ 14 rad/s`. At `dt = 0.01`, `|λ|·dt ≈ 0.14` — comfortably stable.
This is a deliberate design choice: by *not* modeling unsprung tire mass in the simple version, you
avoid the ~620 rad/s pole that forces `dt=0.001` in the real UH-60 model.

**Friction regularization IS the stiff part.** Inside the boundary layer, friction acts as a
viscous damper with effective coefficient `c_eff = mu*N/v_bl`. With `mu=0.6`, `N≈3270 N` per leg:

| `v_bl` | `c_eff` (3 legs) | `λ = -c_eff/m` | RK4 stable to | verdict |
|---|---|---|---|---|
| 0.01 m/s | 588,600 N·s/m | **−589 rad/s** | dt < 4.7 ms | **needs dt=0.001**; unstable at 0.01 |
| 0.05 m/s | 117,720 | −118 rad/s | dt < 23 ms | stable at 0.01, but 5 cm/s of false creep |

**This is the whole experiment in one table.** Tight `v_bl` = accurate stick, needs fine `dt`.
Loose `v_bl` = runs at coarse `dt`, but the aircraft creeps when it should be stuck. The question
your ML model has to answer is whether it can deliver column-1 fidelity at column-2 cost.

## 7. Sanity checks to build in from day one

These are cheap, and they are exactly the checks that would have caught the `LandingGear.m` bug:

1. **Static equilibrium:** at rest on the ground, `sum(N_i)` must equal `m*g` to ~1e-6. Assert it.
2. **Level drop stays level:** drop with zero attitude/rates → `phi, theta, psi` must remain ~0. Any
   drift means a spurious moment, i.e. a geometry or frame error.
3. **Energy conservation:** set `c=0`, `mu=0`, drop it — total energy (KE + PE + spring PE) must be
   conserved. This catches integrator and sign errors instantly.
4. **Frame discipline:** convert every force to body frame *before* summing (step 10 above), and
   take moments about the **contact point** (step 11). This is the bug you just fixed — don't
   rebuild it.
5. **No-adhesion check:** assert `N >= 0` always.

## 8. File structure

```
SimpleLandingModel/
  landing_params.m        % one struct: mass, inertia, geometry, k, c, mu, v_bl, dt
  landing_eom.m           % xdot = landing_eom(t, x, P, frictionFcn)   <- §3
  contact_forces.m        % per-leg N, friction, force/moment in body   <- §5
  friction_physics.m      % regularized Coulomb
  friction_ml.m           % learned model, IDENTICAL signature
  rk4_step.m              % fixed-step RK4
  run_landing_sim.m       % driver: init, loop, log
  plot_landing_results.m  % xyz vs t, phi/theta/psi vs t, N per leg, contact flags
  run_experiment.m        % the 3-run comparison (§10)
  check_sanity.m          % the §7 assertions
```

The only coupling between the friction model and everything else is the function handle passed into
`landing_eom`. That is what makes the ML/physics comparison clean.

## 9. Build order — simplest working thing first

| step | what | pass criterion |
|---|---|---|
| **1** | EOM + gravity only, **no contact**. Drop from 5 m. | Free-fall matches `z = z0 + ½gt²` analytically |
| **2** | Add normal spring-damper, **no friction** (`mu=0`). Vertical drop, level. | Settles to 5 cm static deflection; stays level; sanity checks 1-3 pass |
| **3** | Add friction, physics model, `dt=0.001`. Drop with `u = 3 m/s`. | Slides, decelerates, **stops and stays stopped** (no creep, no chatter) |
| **4** | Asymmetric case: initial `phi = 5°`. | One front leg contacts first, body rolls level, settles |
| **5** | Wire in `run_experiment.m` and the 3-run comparison | See §10 |

Do not move to the next step until the pass criterion holds. Step 2 is where frame/geometry bugs
surface; step 3 is where the stiffness shows up.

## 10. The experiment design — the part most likely to be gotten wrong

**Do not** just run "physics at 0.001" vs "ML at 0.01" and compare. That conflates *friction model
error* with *integration timestep error*, and you won't know which one you measured.

Run **three** cases from identical initial conditions:

| run | friction | dt | role |
|---|---|---|---|
| **A** | physics, tight `v_bl` | 0.001 | **reference / ground truth** |
| **B** | physics, same `v_bl` | 0.01 | **control** — shows what coarse dt alone costs you (expect instability or gross creep). *This run is the justification for the whole project.* |
| **C** | **ML** | 0.01 | **the proposal** — does it recover A's fidelity at B's cost? |

Compare A↔B and A↔C on: final resting position, stopping distance, peak leg load, settle time,
attitude history. **C is a success only if it beats B and approaches A.**

**What the ML model should learn — important.** Do *not* train it to predict instantaneous friction
force from instantaneous velocity. That reproduces the same stiff dynamics and buys you nothing at
coarse `dt` (this is the same conclusion the UH-60 profiling work reached). Instead train it to
predict the **effective friction over a 0.01 s step** given the state at the start of that step —
i.e. a coarse-grained flow map that absorbs the sub-step stick-slip behavior. Training data comes
from run A: integrate at 0.001, aggregate friction impulse over each 0.01 window, and learn
`(v_t, N, ...) → effective F_t` over that window.

## 11. Scale-up path

Roughly in order of value per unit of effort:

1. **Sloped/uneven terrain** — replace the `z=0` plane with `z_ground(x,y)` and a terrain normal.
   Contact and friction then need the terrain frame, which is exactly where the UH-60 bug lived.
   Reuse the §7.4 frame discipline.
2. **Per-leg compression DOF** — give each leg its own state (`s`, `sdot`) and unsprung mass. This
   is what makes the system genuinely stiff (as in the real model) and is the honest test of whether
   ML can handle stiffness, not just regularization.
3. **Nonlinear leg** — gas-spring force law (`F = Cpv/(Cpv/F0 - s)`), stroke limits, preload —
   matching the UH-60's actual model.
4. **Better friction physics** — LuGre or brush model with internal bristle state. This is also the
   scenario where a friction surrogate becomes *worth it on cost grounds*, since these models are
   50-100× more expensive than regularized Coulomb.
5. **Anisotropic tires** — different longitudinal vs lateral `mu`, rolling resistance.
6. **Rotor thrust input** — a simple vertical force + moment at the hub, so you can fly approaches
   rather than only drop-testing.
7. **Per-leg `k`/`c` variation**, CG offsets, mass variation — parameter sweeps for ML training data
   diversity.

Steps 1 and 2 together get you a model that reproduces the qualitative behavior of the real UH-60
gear while staying small enough to reason about — that's the natural stopping point for this study.
