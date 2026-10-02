# Simplified Landing Model — Steps 1–3 Complete

Minimal 6-DOF rigid-body landing simulator, built per `PLAN.md`. Steps 1–3 are
implemented and verified: gravity/EOM, spring-damper contact, and regularized Coulomb
friction on planar (optionally inclined) terrain. The ML friction model is not started.

**Run everything:** `matlab -batch "run_all_tests"`

## Status

| step | scope | status |
|---|---|---|
| 1 | EOM + gravity, no contact | **PASS** |
| 2 | Normal spring-damper contact, no friction | **PASS** |
| 3 | Friction on a 5° incline: hold, slide-to-stop, integrator stability | **PASS** |
| 4+ | ML friction comparison, scale-up | not started |

## Verification results

**Step 1 — free-fall vs analytic parabola.** RK4 integrates constant-acceleration
dynamics exactly, so this is a machine-precision test, not an engineering tolerance.

| case | position error | attitude drift | body rates |
|---|---|---|---|
| 1a level, at rest | 1.33e-13 m | 0 | 0 |
| 1b level, translating | 1.33e-13 m | 0 | 0 |
| 1c tilted, translating | 2.78e-13 m | 0 | 0 |

Case 1c is the meaningful one: with gravity resolved through a rotated body frame and
no moments acting, inertial acceleration must still be exactly `[0;0;g]` and attitude
must stay frozen. That validates `rot_b2i` and the `R'*[0;0;mg]` gravity transform
together.

**Step 2 — contact.**

| check | result |
|---|---|
| 2A settle matches exact equilibrium | dN = 6.2e-6 N, dθ = 1.8e-11 rad, dz = 1.6e-11 m |
| 2B force balance | \|Σ N − mg\| = 1.9e-6 N |
| 2B moment balance | \|M\| = 1.8e-5 N·m |
| 2C hold at equilibrium, 2 s | drift 2.3e-11 m / 2.6e-11 rad |
| 2D roll and yaw under symmetric drop | **identically 0** |
| 2E energy, c=0, continuous contact | rel drift 1.9e-13 |
| 2F energy, c=0, through bounce | rel drift 2.2e-6 |
| no adhesion (N ≥ 0) | holds in every run |

Steps 1–2 reproduce **bit-identically** after the terrain generalization in step 3.

**Step 3 — friction on a 5° incline** (μ = 0.6, **v_bl = 0.1 ft/s = 0.03048 m/s to match
the UH-60 model's `LG_VBL`**, terrain rises toward North). Validation runs use the planned
physics timestep dt = 0.001, gated on the linearization confirming it stable (max|R| = 0.998).

| check | result |
|---|---|
| 3A creep speed vs `v_bl·tanα/μ` | 4.444424e-3 m/s, rel err 2.4e-8 |
| 3A Σ N vs `mg cos α` | 9772.6700 N, err 4.4e-7 N |
| 3A friction vs `mg sin α`, uphill | 854.9978 N, err 1.4e-6 N |
| 3A force / moment balance | 1.4e-6 N / 6.3e-6 N·m |
| 3A holding, not sliding | max \|F_t\|/(μN) = 0.146 |
| 3B slide downhill from 3 m/s | 0.8984 vs 0.8975 m analytic (0.10%) |
| 3B slide uphill from 3 m/s | 0.6662 vs 0.6690 m analytic (0.42%) |
| 3C stability prediction vs seeded runs | matches in all 5 configurations |
| Coulomb bound \|F_t\| ≤ μN | holds in every run |

## The step-3 stability finding — read this before designing the ML experiment

`PLAN.md` §6 sized the friction stiffness from translation alone (−589 rad/s) and was
wrong. Friction acts at feet 1.8 m below the CG, so it also resists rotation, and the roll
inertia is small (effective mass ~171 kg for lateral friction). At the original
v_bl = 0.01 the fastest pole is **−3,454 rad/s**, which made the planned dt = 0.001 reference
unstable. **Resolved by raising v_bl to the UH-60 value (0.03048 m/s)**; the pole scales as
1/v_bl and drops to −1,143 rad/s.

| v_bl | dt | role | predicted | seed growth | creep vs analytic |
|---|---|---|---|---|---|
| **0.0305** | **0.001** | **physics reference** | **stable** | 0.06× | exact |
| 0.0305 | 0.01 | control run B | unstable | 298× | 4.8× too fast |
| 0.01 | 0.001 | old value | unstable | 8× | — |
| 0.01 | 0.0005 | old value | stable | 0.07× | — |
| 0.15 | 0.01 | loosened until stable | stable | 0.03× | +0.2% |

- **The physics reference is now valid**: v_bl 0.0305 at dt 0.001 creeps at exactly the
  analytic 4.44 mm/s.
- **Control run B is still unstable, even with the UH-60 value.** At dt = 0.01, physics
  friction cannot creep slower than ~19 mm/s (smallest stable v_bl = 0.130). Below that it
  doesn't diverge — friction saturates and it chatters at ~21 mm/s *whatever v_bl you asked
  for*. So the ML target is concrete: **4.4 mm/s of creep at dt = 0.01, where physics
  delivers ~21.**
- The creep column is informational for the `v_bl = 0.01` rows: they start from the
  `v_bl = 0.0305` equilibrium and are still converging after 2 s.
- **A symmetric test hides this.** Run unseeded, v_bl 0.01 / dt 0.001 creeps at exactly
  the analytic value with zero lateral motion, because symmetry never excites the unstable
  mode. `test_step3` seeds a 1e-4 m/s lateral velocity on purpose.

### Largest stable timestep (`sweep_dt_stability.m`)

dt swept upward from 1 ms on the 5° hold case, current parameters (μ 0.6, v_bl 0.1 ft/s), then
bisected by simulation:

| seed | largest stable dt (simulated) | linear prediction |
|---|---|---|
| lateral 1e-4 m/s — **the real limit** | **2.437 ms** | 2.436 ms |
| along-slope only — what a symmetric test reports | 10.482 ms | 10.479 ms |

Margin at the 1 ms reference is 2.44×; the 10 ms ML target is 4.1× beyond the real limit. A
symmetric test would wrongly call 10 ms stable. Past the limit, creep grows roughly linearly with
dt (~2.1 mm/s per ms): the integrator, not v_bl, sets the effective boundary layer.

### Drop-test timestep sweep (`sweep_dt_drop.m`)

A realistic landing: nose-first **downhill** at 3 m/s onto the 5° incline, released with a **3° roll
offset** so one front wheel touches first (the asymmetry that excites the lateral/roll friction
mode). 4 s run: fall → touchdown 0.16 s → slide → rocking → hold. Reference dt = 0.5 ms, itself
within 1.8e-5 m / 0.016% of a 0.25 ms run.

| | result |
|---|---|
| **Largest stable dt** | **2.4367 ms** (bisection by simulation); linear prediction at the landed state 2.4364 ms (0.01%) |
| Accuracy for dt ≤ 2.4 ms | stop distance ≤ 1.9e-4 m of 1.020 m, peak leg load ≤ 0.39%, final position ≤ 1.1 mm |
| dt = 10 ms (ML target) | slide and touchdown still accurate (stop 4.3e-4 m, peak load 0.88%), but hold creeps at 21.3 mm/s and ends 71 mm off |

**The landing is resolved at any dt up to 12 ms; the grip is not.** Friction is only stiff when foot
slip is below v_bl — during the hold, and briefly at each rocking reversal. That is precisely the
regime an ML friction model at dt = 10 ms has to get right. Same limit as the steady-hold sweep,
confirming that sweep was a valid predictor *provided* it is laterally seeded.

Note: at 4 s the reference is still settling (creep 4.85 mm/s, rotating at 7.8e-4 rad/s). Run to
20 s it converges to 4.4444 mm/s, the analytic value to 3e-12. The late creep in this sweep is
therefore compared against the reference, not the analytic value.

## Earlier findings worth carrying forward

**1. `PLAN.md` Sec.9's "stays level" criterion was wrong.** Feet sit at body
x = +2, +2, −3, so `sum(x_i) = +1 m` — the leg triangle's centroid is not under the CG.
Equal loads would leave a 3270 N·m pitch moment. The true equilibrium loads unevenly:

- front legs 2937.62 N (δ = 0.04492 m)
- aft leg 3934.77 N (δ = 0.06016 m)
- nose-up pitch θ = 0.003049 rad (0.175°)
- **mean δ = 0.05000 m**, so `k` is correct as specified

Geometry *is* laterally symmetric, so roll and yaw are exactly zero and are enforced to
1e-10. Pitch is legitimately nonzero.

**2. A first hand-derived equilibrium check disagreed with the sim by ~11 N — the sim
was right.** The exact per-leg pitch moment is `N_i·(z_i·sin θ + x_i·cos θ)`; the crude
check used only `sum(x_i·N_i) = 0`, dropping the `z_i·sin θ` term. Because the feet sit
1.8 m below the CG, a pitch tilt gives the vertical normal forces a moment arm about y,
worth 53.8 N·m here. `test_step2_contact.m` now solves the exact 2-unknown equilibrium
by Newton iteration, independently of the model code, and agrees to 6e-6 N.

## Files

| file | role |
|---|---|
| `rot_b2i.m` | body→inertial rotation — **the single definition**, never duplicate it |
| `landing_params.m` | all parameters in one struct |
| `landing_eom.m` | 12-state EOM, `xdot = f(F, M)` |
| `make_terrain.m` | planar terrain (slope, uphill heading) → normal vector |
| `contact_forces.m` | per-leg contact relative to the terrain plane → body-frame force/moment |
| `friction_physics.m` | regularized Coulomb, in the terrain tangent plane |
| `linearize_eom.m` | numerical Jacobian + eigenvalues, for integrator stability prediction |
| `rk4_step.m` | fixed-step RK4 |
| `run_landing_sim.m` | driver with full logging |
| `landing_energy.m` | mechanical energy, for the conservation checks |
| `check_sanity.m` | `PLAN.md` Sec.7 safety checks |
| `test_step1_freefall.m` / `test_step2_contact.m` / `test_step3_friction_slope.m` | verification suites |
| `run_all_tests.m` | runs all three, writes figures |
| `plot_landing_results.m` | PNG output to `figures/` |

## Frame conventions — read before editing

- **NED inertial**: x North, y East, **z DOWN**. **Altitude = −z.**
- **Terrain**: a plane with unit normal `P.terrain.n_down` (pointing *into* the ground).
  Penetration is measured along that normal, not along z: a foot is penetrating when
  `n_down·(P_foot − p0) > 0`. Flat ground is `n_down = [0;0;1]`.
- **Body**: x forward, y right, z down, origin at CG. Euler 3-2-1.
- Every per-leg force is rotated into body frame **before** summing.
- Moments are taken about the **contact point**.

The last two rules are the direct lesson from the UH-60 `LandingGear.m` frame/moment-arm
bug. Preserve them when adding the leg compression DOF in scale-up, at which point
`p_foot_b` becomes state-dependent.

## Next: the ML experiment

The three-run design in `PLAN.md` §10 is now valid as written: physics at dt = 0.001
(v_bl = 0.1 ft/s, matching the UH-60) is the reference, physics at dt = 0.01 is the
control, ML at dt = 0.01 is the proposal. Seed lateral motion in every
stability-relevant run — a symmetric run hides the instability.
