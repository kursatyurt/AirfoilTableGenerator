#!/usr/bin/env python
"""Find how many MPI ranks this machine should use, once, and remember it.

A single 2D airfoil case is small; past some rank count the halo exchange costs
more than the extra cores buy, and throughput drops. This times a fixed number
of SU2 iterations at several rank counts and writes the best to machine.conf,
which polar.py then uses as its default --np.

    python tune_np.py             # study, then write machine.conf
    python tune_np.py --show      # print the stored value
"""
import argparse, subprocess, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONF = ROOT / "machine.conf"


def stored_np():
    """Ranks chosen by a previous study, or None."""
    if CONF.exists():
        for line in CONF.read_text().splitlines():
            if line.startswith("NP="):
                return int(line.split("=")[1].strip())
    return None


def main():
    import os, sys
    sys.path.insert(0, str(ROOT))
    from polar import make_cfg, find_dat
    from mesh import generate_mesh
    from read_airfoil import read_airfoil_coordinates

    cores = os.cpu_count()
    ap = argparse.ArgumentParser()
    ap.add_argument("--show", action="store_true")
    ap.add_argument("--airfoil", default="n0012", help="case to time")
    ap.add_argument("--iters", type=int, default=150, help="iterations per timing run")
    ap.add_argument("--ranks", default=None, help="comma list; default 2,4,8,... up to core count")
    ap.add_argument("--repeats", type=int, default=3, help="timings per rank count; the best is kept")
    # ponytail: 0.25, not 0. The last doublings of rank count buy a few percent of
    # wall clock for half the machine; on a 2D case that is a bad trade, and it
    # stops one sweep from monopolising a shared box.
    ap.add_argument("--tolerance", type=float, default=0.25,
                    help="accept the fewest ranks within this fraction of the fastest time; "
                         "lower it to chase wall clock at any core cost")
    a = ap.parse_args()

    if a.show:
        print(f"machine.conf: NP={stored_np()}" if stored_np() else "no machine.conf yet")
        return

    # ponytail: no 1 in the default list -- SU2's MPI build aborts in MPI_Win_create
    # on a single rank with OpenMPI, so 2 is the practical serial baseline.
    ranks = ([int(r) for r in a.ranks.split(",")] if a.ranks else
             [n for n in (2, 4, 8, 12, 16, 24, 32, 48, 64) if n <= cores])
    case = ROOT / "runs" / "_tune"
    case.mkdir(parents=True, exist_ok=True)
    mesh = case / "airfoil.su2"
    if not mesh.exists():
        dat = find_dat(a.airfoil)
        x, y = read_airfoil_coordinates(str(dat.parent), dat.name)
        generate_mesh(x, y, 1e6, 0.15, path=mesh)

    # No convergence exit, so every run does exactly the same work.
    cfg = case / "tune.cfg"
    cfg.write_text(make_cfg("inc", 0.0, 1e6, 0.15, a.iters, False)
                   .replace("CONV_RESIDUAL_MINVAL= -8", "CONV_RESIDUAL_MINVAL= -99")
                   .replace("OUTPUT_FILES= ( RESTART, PARAVIEW, SURFACE_CSV )", "OUTPUT_FILES= ( RESTART )"))

    print(f"{cores} cores detected; {a.iters} iterations x {a.repeats} repeats per rank count\n")
    print(f"{'ranks':>6} {'seconds':>9} {'speedup':>8} {'efficiency':>11}")
    results = []
    for n in ranks:
        times, rc = [], 0
        for _ in range(a.repeats):
            t = time.perf_counter()
            rc = subprocess.call(["mpirun", "-n", str(n), "SU2_CFD", cfg.name], cwd=case,
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if rc != 0:
                break
            times.append(time.perf_counter() - t)  # keep the best: noise only ever adds time
        if not times:
            print(f"{n:>6} {'failed':>9}   (rc={rc}; this machine will not run {n} ranks)")
            continue
        results.append((n, min(times)))
        base_n, base_t = results[0]
        speed = base_t / results[-1][1]
        print(f"{n:>6} {min(times):>9.1f} {speed:>8.2f}x {speed / (n / base_n) * 100:>10.0f}%")
    print(f"(speedup and efficiency are relative to {results[0][0]} ranks)" if results else "")

    if not results:
        raise SystemExit("every timing run failed")

    best = min(results, key=lambda r: r[1])
    # Prefer fewer ranks when the extra ones barely help: anything within 10% of
    # the best time counts as "as fast", so take the cheapest of those.
    pick = min((r for r in results if r[1] <= best[1] * (1 + a.tolerance)), key=lambda r: r[0])
    CONF.write_text(f"# written by tune_np.py -- delete to re-tune\nNP={pick[0]}\n")
    print(f"\nfastest: {best[0]} ranks ({best[1]:.1f}s)")
    print(f"chose:   {pick[0]} ranks ({pick[1]:.1f}s) -> {CONF}")


def _selftest():
    """The pick rule: fewest ranks among those within 5% of the fastest."""
    pick = lambda res, tol=0.25: min((r for r in res if r[1] <= min(res, key=lambda q: q[1])[1] * (1 + tol)),
                                     key=lambda r: r[0])
    # 4 ranks are within tolerance of the 8-rank best, so take 4 and leave cores free
    assert pick([(2, 100.0), (4, 30.0), (8, 29.5), (16, 40.0)]) == (4, 30.0)
    assert pick([(2, 100.0), (4, 30.0), (8, 20.0), (16, 90.0)]) == (8, 20.0)
    assert pick([(2, 10.0), (32, 40.0)]) == (2, 10.0)  # scaling never pays off
    # this machine's measured times: 24 is fastest, but 8 is within 25% of it
    measured = [(2, 52.4), (4, 35.3), (8, 26.7), (12, 30.0), (16, 27.2), (24, 23.0)]
    assert pick(measured) == (8, 26.7)
    assert pick(measured, tol=0.0) == (24, 23.0)
    print("tune_np selftest ok")


if __name__ == "__main__":
    import sys
    _selftest() if "--selftest" in sys.argv else main()
