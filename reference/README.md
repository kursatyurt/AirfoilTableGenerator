# Reference data

`ADA033425_helicopter_design_datcom_vol1_airfoils.pdf` — *US Army Helicopter
Design Datcom, Volume I: Airfoils*, USAAMRDL CR 76-2, Boeing Vertol, September
1976 (DTIC AD-A033425). A US Government work, "Distribution Statement A:
approved for public release, distribution unlimited", so it is redistributable.
23 MB, which is most of this repository's size.

Sheets of interest:

| sheet | PDF page | content |
|---|---|---|
| 1.3.20 | 153–155 | NACA 0012, NPL transonic tunnel, M 0.30–0.85, Re 1.7–3.75e6 |
| 1.3.30 | 156–158 | NACA 0012 through 180° AoA, Re 1.8e6 (Critzos, NACA TN 3361) |
| 1.3.60 | 165–167 | NACA 0015, Ames 1×3.5 ft tunnel, M 0.30–0.825, Re 1–2e6 |

## What was extracted, and how far to trust it

**`naca0012_coordinates.csv`, `naca0015_coordinates.csv` — trustworthy.**
Transcribed from printed numeric tables, so they are exact to the precision
shown. `check_reference.py` verifies both against the analytic NACA 4-digit
thickness distribution: worst deviation is 0.00011 chord for the 0012 and
0.00001 for the 0015, which confirms no digit was misread.

**`naca0012_stall_landmarks.csv` — approximate, with stated bounds.**
The polar data in this report exists only as plotted curves, and the scan is
poor: the halftone screen printed as dark as the ink, and the lower half of each
figure is nearly solid black. Values were read against a pixel grid calibrated
to the printed axis labels — the calibration was verified by overlaying the grid
and confirming it lands on the printed "1.2", "0.8", "0.4", "0" and "0", "20.",
"40." labels — but they are still eye-read from a curve. The uncertainty columns
bound the reading, not the experiment.

## NACA Report 832 — the original source, and a much better scan

The Datcom cites **NACA Report 832** (Graham, Nitzberg & Olson, 1945) for the
0015 data. NASA's scan of it ([NTRS
19930091909](https://ntrs.nasa.gov/citations/19930091909)) is crisp line art on
white, not a halftone photocopy, and its figure 27 carries all fourteen Mach
curves. `digitize_tr832.py` fetches that PDF (56 MB, gitignored) and extracts:

- `naca0015_lift_slope_vs_mach.csv` — lift slope for **all 14 Mach numbers**,
  0.300 to 0.825
- `naca0015_clmax_vs_mach.csv` — CLmax and stall angle for M = 0.300, 0.400,
  0.500 only

Three independent checks back this, all run by `check_reference.py`:

1. **Symmetry.** Every curve must cross CL = 0 at its own α = 0. All fourteen
   fitted zero-lift angles land within 0.21°, confirming the axis calibration,
   the 4° stagger, and each curve's identity.
2. **Compressibility.** Nothing in the digitiser knows about Mach — each curve
   is fitted independently — yet the slopes track Prandtl-Glauert to within 3%
   up to M = 0.675 and then fall away (−19% at 0.70, −44% at 0.825), exactly the
   force-divergence behaviour a real airfoil must show.
3. **Cross-source.** The M = 0.30 slope from this scan (0.1000/deg) and from the
   Datcom scan (0.1003/deg) agree to **0.3%** — two independent digitisations of
   two different reproductions of the same 1945 experiment.

CLmax stops at M = 0.500 because above it the curves crowd together and the peak
tracer jumps to a neighbour. The lift slopes are unaffected: they come from the
straight lower part of each curve, which stays well separated, and the overlay
was checked visually for all fourteen.

### Drag and moment: attempted, not delivered

Figure 32 (page 50) is the NACA 0015 drag polar and figure 37 the moment. Both
were attempted and neither is delivered, for a specific and checkable reason.

These figures stack each Mach curve on its **own cd = 0 baseline** rather than a
shared axis, so every value depends on knowing which baseline belongs to which
Mach. The evidence does not settle that:

- the y-axis carries **11** labelled "0" ticks, evenly spaced 158.7 px apart
  (= 0.02 in cd, fixed by the graduated 0.02–0.12 scale above the top curve)
- the plot frame spans 223 → 3005 px, which leaves room for only **12** such
  baselines
- but the figure carries **14** Mach labels, 0.300 through 0.825

Twelve slots, fourteen curves. Guessing wrong by one baseline shifts every cd by
0.02 — about **twice** the actual drag of this section at low Mach — and the
result would still look entirely plausible. That is the worst kind of error, so
no drag numbers are recorded.

### Worth a second look: figures 40, 41 and 42

Page 54 carries three figures on ordinary, unstaggered axes, which makes them far
more tractable than the carpet plots:

| figure | content | why it matters |
|---|---|---|
| 40 | lift-curve slope vs Mach at cl = 0.20, all five airfoils | an independent check on the slopes extracted here |
| 41 | angle of zero lift vs Mach, all five | should read 0 for the 0015 at every Mach |
| 42 | **section drag coefficient vs Mach at cl = 0.20**, all five | the absolute drag reference figure 32 could not give |

The obstacle is different here and may be solvable: the five airfoils are drawn
as five dash patterns that nearly coincide below M = 0.7, so the work is telling
the 0015 curve from its neighbours rather than calibrating an axis.

**Not extracted from the Datcom sheets: every Mach number other than 0.30, and
all CD/CM polars.** This was attempted seriously and abandoned. Six
substantively different methods were tried:

1. intensity-following trace of a single curve — wandered onto background ink
2. trace restricted to solid black runs — same, the lower scan region is solid
3. slope-predicting tracker over all 14 Mach curves — merged adjacent curves
   (the giveaway: two Mach numbers fitting to identical intercepts)
4. multi-target tracking with exclusive assignment — worse; tracks lost their
   curve and coasted along gridlines
5. row-crossing detection using a predicted slope — **circular and invalid**: the
   search window assumed 0.10/deg, so every Mach "measured" 0.09–0.10/deg and
   flatly contradicted Prandtl-Glauert. Caught precisely because the extracted
   slopes did not rise with Mach as physics requires
6. row-crossing with gridlines discriminated by their lack of motion — the
   halftone screen makes row thresholding useless; the same figure yielded
   17 to 79 "crossings" depending on which row was sampled

The M=0.30 curve is extractable only because it is the leftmost of a staggered
carpet plot with clear space to its left. The remaining thirteen sit in a fan
2° apart, cross each other past stall, and are interleaved with gridlines of the
same darkness and the same 2° spacing.

CLmax for the next two or three Mach numbers is visually readable, but assigning
a peak to the correct Mach means counting curves through a cluttered region, and
**there is no independent physical check on that assignment** — unlike the
zero-lift test that caught a 2° axis error on the M=0.30 curve. An off-by-one
there would mislabel the Mach number, so those values are not recorded either.

For anything beyond M=0.30, go to the original sources the Datcom cites: NACA
Report 832 (Graham, Nitzberg & Olson) for the 0015, NACA TN 3361 for the 0012
180° data, ARC C.P. 1261 (Gregory & Wilby) for the NPL 0012.

## Using the landmarks

The 180° sheet is the useful one for a stall study: single Reynolds number
(1.8e6), fully turbulent, one curve per plot. Compare a computed sweep against
it with the caveat that steady RANS is not expected to reproduce post-stall
behaviour — see the stall notes in the top-level README.

```
python reference/check_reference.py    # verifies the coordinate tables
```
