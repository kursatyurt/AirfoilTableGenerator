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

Scaling is poor past 8 ranks: 24 is the fastest but buys 16% over 8 ranks for
three times the cores, and 12 is reproducibly *slower* than 8. Use
`--tolerance 0.25` to bias the choice toward fewer ranks, or pass `--np`
explicitly to override `machine.conf` entirely.

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
| `--regime` | `inc` | `inc` = `INC_RANS` for M < 0.3; `comp` = compressible `RANS`, use for transonic. |
| `--np` | from `machine.conf` | MPI ranks. |
| `--iters` | `2000` | raise it if `polar.csv` shows `converged=0`. |
| `--yplus` | `1.0` | sets the wall spacing and the number of normal layers. |

Turbulence is Spalart–Allmaras, fully turbulent from the leading edge — no
transition model, so drag at low Re is overpredicted and this will not capture a
laminar bucket. Stall angles are not trustworthy either; steady RANS on a 2D
airfoil stops being meaningful once the flow separates.

## What comes from where

- `mesh.py` — structured C-mesh, forked from FALCON's `meshing.py`. See the
  module docstring for what changed; the important one is that an open trailing
  edge is now closed at the midpoint of the two TE points rather than onto the
  upper one, which was skewing symmetric airfoils.
- `opt/FALCON` — cloned by `INSTALL.sh` for the Selig `.dat` database and
  `read_airfoil.py`. Its GUI is unused.
- `opt/su2` — SU2 built from source with MPI by `INSTALL.sh`.

## Validation

NACA 0012, Re 1e6, M 0.15, y+ = 1, converged to `rms[P] = -8`:

| | value | expected |
|---|---|---|
| CL at 0° | 6.0e-05 | 0 by symmetry |
| CMz at 0° | -4.2e-05 | 0 by symmetry |
| CD at 0° | 0.0107 | ~0.011 fully turbulent at this Re |
| dCL/dα | 5.9 /rad | 2π thin-airfoil |

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
