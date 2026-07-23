# AirfoilTableGenerator

Airfoil polars (CL, CD, CM vs angle of attack) from [SU2](https://github.com/su2code/SU2),
one command per airfoil. Meshing is a fork of [FALCON](https://github.com/Prisha22/FALCON)'s
gmsh C-mesh generator.

```bash
bash INSTALL.sh          # ~30 min, mostly the SU2 build. No sudo.
source env.sh
python tune_np.py        # once per machine: measures the right MPI rank count
python polar.py --airfoil n0012 --re 1e6 --mach 0.15 --aoa -4:16:2
```

Outputs land in `runs/<airfoil>/`: `airfoil.su2`, one `history_<aoa>.csv` +
`aoa_<aoa>.log` per angle, and `polar.csv` / `polar.png`.

## Picking the number of cores

More ranks is not better. A 2D airfoil case is small, and past a point the halo
exchange costs more than the extra cores buy. `tune_np.py` times a fixed number
of SU2 iterations at several rank counts, three repeats each, and writes the
best to `machine.conf`; `polar.py` then uses it as the default `--np`. Run it
once per machine — it is a few minutes, and `machine.conf` is gitignored because
the answer is machine-specific.

Measured on a 32-core box, 112k-point mesh:

| ranks | seconds | speedup | efficiency |
|---|---|---|---|
| 2 | 52.4 | 1.00x | 100% |
| 4 | 35.3 | 1.49x | 74% |
| 8 | 26.7 | 1.96x | 49% |
| 12 | 30.0 | 1.75x | 29% |
| 16 | 27.2 | 1.93x | 24% |
| 24 | 23.0 | 2.28x | 19% |
| 32 | — | fails | — |

Scaling is poor past 8 ranks: 24 is the fastest in raw wall clock but buys only
16% over 8 ranks for three times the cores, and 12 is reproducibly *slower* than
8. So the picker takes the **fewest ranks within 25% of the fastest time**
(`--tolerance`), which chooses 8 here — the knee, and it leaves the machine
usable. Drop to `--tolerance 0` to chase wall clock regardless of core cost, or
pass `--np` to override `machine.conf` entirely.

SU2 built against OpenMPI aborts in `MPI_Win_create` on a single rank, so 2 is
the practical minimum and the study starts there.

## Warm starts

Each angle restarts from the previous angle's converged flow field
(`restart_flow.dat` is copied to `solution_flow.dat` and `RESTART_SOL=YES` is
set), so only the first point of a sweep starts cold. Neighbouring angles differ
by a couple of degrees, so this cuts most of the sweep's cost. Delete
`runs/<airfoil>/` to start over from scratch, including the mesh.

## Options

| flag | default | note |
|---|---|---|
| `--airfoil` | — | a path to a `.dat`, or a name from `opt/FALCON/Airfoil_DAT_Selig` (1624 airfoils, Selig format). Misses print near matches. |
| `--re` / `--mach` | `1e6` / `0.15` | chord is 1 m; sets viscosity (`inc`) or the SU2 `REYNOLDS_NUMBER`/`MACH_NUMBER` (`comp`). |
| `--aoa` | `-4:16:2` | `lo:hi:step` inclusive, or `0,2,4`. |
| `--regime` | `comp` | `comp` = compressible `RANS` (Roe + low-Mach preconditioning), the default and consistent across a whole Mach sweep; `inc` = `INC_RANS`, only for strictly incompressible cases. |
| `--np` | from `machine.conf` | MPI ranks. |
| `--iters` | `10000` | max iterations per angle; the run stops earlier when the coefficient Cauchy criterion is met. |
| `--yplus` | `1.0` | sets the wall spacing and the number of normal layers. |
| `--farfield` | `15` | farfield radius in chords. Push it well past 15 for transonic. |
| `--transition` | `none` | `lm` adds Langtry-Menter laminar-turbulent transition on top of SA. |
| `--tu` | `0.001` | freestream turbulence intensity the transition model keys off. |

## Convergence

Each angle converges on the **force coefficients**, not on a residual field:
`CONV_FIELD= ( LIFT, DRAG, MOMENT_Z )` with a Cauchy tolerance of
`CONV_CAUCHY_EPS= 1E-4` over the last `100` iterations — i.e. the run stops once
lift, drag and pitching moment have each settled to within **1E-4 (one drag count
in C_D)**. This is what a coefficient sweep actually cares about; a residual can
reach its floor while the integrated lift is still drifting (which produced the
earlier frozen-lift LM sweep). If SU2 hits the `--iters` budget first, the angle
is flagged `converged=0` in `polar.csv`.

## Transition

Turbulence is always Spalart–Allmaras. `--transition lm` adds the Langtry-Menter
γ-Reθ laminar-turbulent transition model on top of it (two extra transport
equations) rather than replacing it:

- `none` (default) — fully turbulent from the leading edge. Overpredicts drag at
  low Re and cannot produce a laminar drag bucket.
- `lm` — Langtry-Menter γ-Reθ. Resolves natural transition and the drag bucket;
  runs at a gentler CFL ceiling (15 vs 50) so it is slower per angle.

`lm` reads `--tu`, the freestream turbulence intensity (`0.001` = 0.1%, a
low-turbulence wind tunnel). Transition location is sensitive to it, so match it
to the experiment you are comparing against.

## What comes from where

- `mesh.py` — structured C-mesh, forked from FALCON's `meshing.py`. See the
  module docstring for what changed; the important one is that an open trailing
  edge is now closed at the midpoint of the two TE points rather than onto the
  upper one, which was skewing symmetric airfoils.
- `opt/FALCON` — cloned by `INSTALL.sh` for the Selig `.dat` database and
  `read_airfoil.py`. Its GUI is unused.
- `opt/su2` — SU2 built from source with MPI by `INSTALL.sh`.

## Validation

NACA 0012, Re 1e6, M 0.15, y+ = 1, converged to `rms[P] = -6`:

| | value | expected |
|---|---|---|
| CL at 0° | 6.0e-05 | 0 by symmetry |
| CMz at 0° | -4.2e-05 | 0 by symmetry |
| CD at 0° | 0.0107 | ~0.011 fully turbulent at this Re |
| dCL/dα | 5.9 /rad | 2π thin-airfoil |

Compressible regime, NACA 0012, AoA 0, `--farfield 15`:

| case | CL | CD | rms[Rho] reached |
|---|---|---|---|
| M 0.15, Re 1e6, 3000 iters | -3.0e-04 | 0.0098 | -6.9 |
| M 0.80, Re 9e6, 4000 iters | -2.4e-04 | 0.0169 | -5.1 |

Symmetry holds in both, and M 0.8 picks up the expected shock drag rise. Two
caveats. The compressible solver gives CD 0.0098 at M 0.15 where the
incompressible one gives 0.0107 — an 8% gap at conditions where the two should
very nearly agree, most likely the ROE vs FDS dissipation, so do not mix regimes
within one study. Both compressible cases reached the `-6` convergence
target only marginally (M 0.8 at −5.1 still misses it); transonic
especially needs more `--iters` and a farfield much larger than the 15-chord
default.

### NACA 0015 against experiment (Ames 1×3.5 ft tunnel, 1945)

Validated against NACA Report 832 / the Army Datcom sheet — M 0.30, Re 1.5e6,
compressible, α = 0:16:2. `validate.py` does the comparison and reports CL slope,
zero-lift angle, CLmax, CD at cl = 0.22 (the one drag point figure 32 resolves),
and a CM symmetry check:

```
python validate.py runs/naca0015_M030/polar.csv reference/naca0015_cl_alpha_M030.csv 0.30
```

| quantity | fully turbulent | `--transition lm` | experiment |
|---|---|---|---|
| lift slope /deg | 0.1056 (+5.4%) | ~0.099 (−1%) | 0.1002 |
| CD at cl = 0.22 | 0.0117 (+50%) | 0.0061 (−22%) | 0.0078 |
| CLmax | 1.30 at 14° | 1.16 at 12° | 1.08 at 12° |
| CM at α = 0 | −0.0005 | −0.0005 | 0 (symmetry) |

The fully-turbulent run is the trustworthy one — every angle converges to
`rms[Rho] = -6`. It sits **+5.4% high on lift slope** and **+50% high on drag**,
both the expected sign: fully turbulent from the leading edge, against a 1945
tunnel model that carried a long laminar run, and with no wind-tunnel-wall or
blockage correction applied to the reference. Turbulent skin friction alone is
~0.012 at this Re, which is essentially the computed CD, so that drag is a
physics-of-the-model result, not a mesh or convergence artefact.

Three things were tested as possible fixes for the +5.4% and **none helped** —
recorded here so they are not re-tried:

- **Farfield.** 15 → 30 chords made the slope *worse* (+5.4% → +7.6%), the
  opposite of the confinement argument. `--farfield 15` stays the default.
- **Grid.** 81k → 115k → 190k cells moved the slope inside a ±0.4% band,
  non-monotonically — grid-converged. The default 401/301/180 mesh is adequate.
- **Wall spacing.** y+ 1.0 → 0.5 moved both slope and CD *further* from
  experiment, confirming the CD is a turbulence-model ceiling. y+ = 1.0 stays.

What *did* move everything the right way is transition. `--transition lm` pulls
the slope to within ~1%, brings CLmax and the stall angle down onto the
experimental 12°, and drops CD through the experimental value (it overshoots low
because LM at Tu = 0.1% laminarises more of the surface than the real model had).
LM's two transport equations can hunt at low α on this case — the transition
front moves by a cell or two — so `--transition lm` runs at a lower CFL ceiling
(15 vs 50 turbulent). Convergence is judged on the coefficients themselves
(`CONV_FIELD= (LIFT, DRAG, MOMENT_Z)`, Cauchy tolerance 1E-4 = one drag count),
so an angle only reports `converged=1` once lift, drag and moment have actually
settled — not merely when a residual floors.

### Stall

NACA 0012, Re 1e6, M 0.15, fully turbulent, warm-started through the sweep:

| α | 8 | 10 | 12 | 14 | 16 | 18 | 20 |
|---|---|---|---|---|---|---|---|
| CL | 0.825 | 1.019 | 1.197 | 1.352 | **1.455** | 1.404 | 1.117 |
| CD | 0.0179 | 0.0225 | 0.0286 | 0.0367 | 0.0471 | 0.0624 | 0.1114 |

CLmax lands at 1.455 at 16°, against 1.30 at 12° in the experiment
(`reference/naca0012_stall_landmarks.csv`) — **12% high and 4° late**. That is
the textbook fully-turbulent steady-RANS result, not a bug: no transition model
means no laminar separation to trigger stall early, and steady RANS cannot
represent the unsteady separated wake that limits real CLmax. Treat the
post-stall points as qualitative. Every point reports `converged=1`, so
convergence is not what is wrong — the physics model is.

### Transition models compared

NACA 0012, Re 1e6, M 0.15, α = 0:

| `--transition` | CL | CD | converged |
|---|---|---|---|
| `none` | +0.00006 | 0.01074 | yes |
| `lm` | −0.00138 | 0.00485 | yes |

Transition roughly halves the drag, as it should at this Reynolds number, and
keeps CL symmetric (≈0 at α=0) as it must.

`python mesh.py`, `python polar.py --selftest` and `python tune_np.py --selftest`
run the unit checks (mesh block symmetry, AoA parsing, SU2 history-CSV parsing,
config generation, rank-choice rule). `INSTALL.sh` runs all three at the end.

## Not included

XFOIL cross-checks, CST/PARSEC parametrisation, and shape optimisation — those
live in FALCON's GUI. Mesh independence has not been studied: the defaults are
401 chordwise / 301 wake / 180 leading-edge nodes with a 15-chord farfield, all
exposed as keyword arguments on `mesh.generate_mesh`.

## License

MIT — see [LICENSE](LICENSE), which also carries FALCON's MIT notice and explains
the SU2 (LGPL-2.1) and gmsh (GPL) situation. Neither is redistributed here;
`INSTALL.sh` fetches and builds them on your machine.
