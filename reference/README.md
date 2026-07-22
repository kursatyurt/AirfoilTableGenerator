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

### Drag: the baseline problem, and how it was settled

Figure 32 (page 50) is the NACA 0015 drag polar and figure 37 the moment. Both
stack each Mach curve on its **own cd = 0 baseline**, so every value depends on
knowing which baseline belongs to which Mach. An earlier attempt stopped here,
and correctly: it counted 11 labelled "0" ticks against 14 Mach curves, and
guessing wrong by one baseline shifts every cd by 0.02 -- more than twice the
true drag of this section at low Mach -- while still looking plausible.

**The count was right; the inference from it was wrong.** Re-measured at 600
dpi, the y axis carries a labelled tick on *every odd grid row*, seventeen in
all -- verified as a clean binary signal, 0.20 to 0.76 ink on every odd row and
0.00 to 0.06 on every even one. Six of those seventeen are the graduated scale
`.12 .10 .08 .06 .04 .02`; the other eleven are all labelled "0". Counting
curves up from the bottom axis, one per zero, reaches M = 0.725 and runs out.
The remaining four -- 0.750, 0.775, 0.800, 0.825 -- **share** the zero at the
top of the graduated scale. That is what the graduated 0-to-0.12 scale is
printed for, and it is forced rather than chosen: giving M = 0.775 its own zero
one row up puts its minimum at cd = -0.002, and no assignment starting anywhere
but the bottom axis for M = 0.300 fits inside the frame.

Two independent checks, both run as assertions by `digitize_tr832_drag.py`:

1. **Figure 42, page 54.** It plots cd against Mach at cl = 0.20 for all five
   sections on an ordinary shared axis, so it cannot suffer this error at all.
   All fourteen values read out of figure 32 land inside its envelope, across a
   factor of thirteen in cd -- 0.0078 at M = 0.300 up to 0.100 at M = 0.825. An
   off-by-one baseline would read 0.028 at M = 0.300 and miss by fourfold.
2. **Symmetry.** The 0015 is symmetric, so cd at cl = +0.22 and at cl = -0.22
   must agree. They do, to 0.0004 or better for every Mach up to 0.650.

`naca0015_cd_vs_mach_tr832.csv` -- cd at |cl| = 0.22 for **all fourteen Mach
numbers**, 0.300 to 0.825, with the symmetry residual as an error bar.

**Minimum drag is not recorded, and neither are the cd-vs-cl polars.** The
bottom of the drag bucket is printed underneath the Mach label, which sits
directly above its own curve inside the same band; reading around it gave
values that were not monotone in Mach at the 0.001 level, which is the size of
the effect. Traced polars failed the symmetry check by up to 0.005 at scattered
cl. Ten of the fourteen were close to usable and none were checkable point by
point, so none are recorded. The moment figure was not attempted.

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

## TR-832 contains only one symmetric section

The title of Report 832 — "five representative NACA low-drag and conventional
airfoil sections" — refers to the NACA 65-215, 66,2-215, 0015, 23015 and 4415,
**not** to a 0006/0009/0012/0015/0018 thickness family. The lift-vs-alpha carpet
plots are figures 25-29 (PDF pages 46-48): fig 25 = 65-215, fig 26 = 66,2-215,
fig 27 = 0015, fig 28 = 23015, fig 29 = 4415. A full-text search of the scan
finds no occurrence of "0006", "0009", "0012" or "0018" anywhere in the report.

So `reference/digitize_tr832.py` already covers every symmetric section this
report has. The other four figures are cambered sections: the zero-lift-angle
symmetry assertion that validates the 0015 extraction does not apply to them,
and no substitute check of comparable strength was available, so they were not
digitised. Thickness-ordering checks (slope and force-divergence Mach vs t/c)
are likewise impossible from a single thickness — they need a second source.

## NACA 23015 (figure 28)

`reference/digitize_tr832_23015.py` extends the figure-27 method to figure 28,
the cambered NACA 23015. The substitute for the symmetry check that the section
above says was missing is thin-airfoil theory on the 230 mean line: it gives
alpha_L0 = -1.09 deg analytically (published experiment about -1.2 deg), and all
fourteen fitted curves land at -0.98 to -1.27 deg, a spread of 0.29 deg. Because
Prandtl-Glauert does not move alpha_L0, that mutual agreement is an independent
test of the 4 deg stagger; a one-gridline axis error would be 2 deg.

Lift slope follows 1/sqrt(1-M^2) to within 2.6 percent up to M = 0.60 and falls
away past it, and the M = 0.30 slope is 0.1000/deg, matching the 0015 to four
figures as camber should.

Only `naca23015_lift_slope_vs_mach.csv` is emitted. No CLmax file: the 23015
curves peak higher and cross one another, and the peak tracer demonstrably jumps
to a neighbour above M = 0.300, so those numbers were dropped rather than shipped
with caveats.

## NACA 4415 (figure 29)

`reference/digitize_tr832_4415.py` applies the same method to figure 29 (bottom
of PDF page 48), the strongly cambered NACA 4415. Calibration is re-derived, not
reused: this figure is tilted -0.20 deg rather than page 47's +0.51 deg, and the
pixel constants are now fitted by least squares off the printed grid itself (23
horizontal rules at 0.1 cl, 36 vertical rules at 2 deg) instead of being
eyeballed. The M = 0.300 origin is the 4th vertical rule and the stagger is
exactly two rules, i.e. 4.000 deg.

The substitute for the symmetry check is thin-airfoil theory on the 4-digit mean
line (m = 0.04, p = 0.4), computed in the script rather than quoted: alpha_L0 =
-4.15 deg, against a published experimental value near -4.0 deg. The seven
accepted curves land at -4.21 to -3.93 deg, mean -4.03, a spread of 0.28 deg.
Lift slope follows 1/sqrt(1-M^2) to within 2.8 percent over M = 0.30 to 0.65, and
the M = 0.30 slope is 0.0994/deg against the 0015's 0.1000 -- 0.6 percent apart,
as it should be, since camber shifts the curve but not its slope.

Two things are deliberately not emitted:

* **Only M <= 0.650.** From M = 0.675 up the shock-induced kink leaves no
  straight segment inside the fitting band and the self-check notices: the
  extrapolated alpha_L0 walks off to -4.40, -4.74, -5.83 deg. Those slopes are
  rejected. The price is that the post-divergence fall-off cannot be shown from
  this figure, only the rise to divergence.
* **No CLmax file at all.** The 4415 knee is far more rounded than the 0015's, so
  the peak tracer flattens near the top, latches onto the stubs the gridline
  stripper leaves behind and runs horizontally into the next curve; overlaying
  the traced points shows all three low-Mach traces converging on one pixel
  (cl = 1.344). Peaks are legible to a human eye (roughly 1.25 to 1.35 for
  M <= 0.55, higher than the 0015 as expected) but nothing in the script can
  check them, so no number is claimed.

The off-by-one-curve hazard is sharper here than for the 23015: alpha_L0 is about
-4 deg and the stagger is exactly 4 deg, so a one-curve index error would produce
a plausible-looking 0 or -8. The alpha_L0 window (-5, -3) excludes both, but on
its own it is not sufficient -- the fitter latches onto whichever curve is nearest
its seed, so shifting the whole carpet by one stagger relabels every curve and
leaves every alpha_L0 unchanged (verified: shifting by +-1 just rotates the list).
What pins the labelling is that M = 0.300 is the leftmost curve, asserted by
re-running the fit one stagger further left, where it finds no curve at all --
only 1574 px of y-axis tick label against 3584 px for the weakest real curve.

## Figures 40, 41 and 42 (printed page 49, PDF page 54)

`digitize_tr832_fig40_42.py` extracts the three Mach-sweep summary figures:
lift-curve slope at c_l=0.20 (fig 40), angle of zero lift (fig 41) and section
drag coefficient at c_l=0.20 (fig 42).  All three plot the same five sections
listed above as five dash patterns on one set of axes, and telling those apart
is the whole job: below M~0.66 the curves lie within a line width of each other.

Only what could be pinned by evidence independent of the dash pattern is
emitted:

| file | contents |
|---|---|
| `tr832_fig41_zero_lift_angle_vs_mach.csv` | NACA 0015 (M 0.30-0.83) and NACA 4415 (M 0.31-0.75) only |
| `tr832_fig40_lift_slope_vs_mach.csv` | five-curve envelope, M 0.305-0.795, no per-airfoil curve |
| `tr832_fig42_drag_vs_mach.csv` | five-curve envelope, M 0.305-0.570, no per-airfoil curve |

Figure 41 is the good one: the 0015 is the only symmetric section and so the
only curve that can sit at zero, and the 4415 is alone in the lower half of the
plot at -4 deg, which fixes both identities without reading a dash pattern.  The
other three sections sit inside a 0.5 deg band around -1.3 deg and are omitted.
For figures 40 and 42 no individual curve could be identified, so the bundle
envelope -- the highest and lowest of the five at each Mach -- is given instead.
Figure 42's post-divergence drag rise is not extracted: past M~0.57 the two
envelope edges merge into a single drawn line and cannot be separated again
without guessing.

Self-checks, all asserted in the script and printed on every run: fig 41 reads
+0.19 deg for the 0015 (symmetry demands 0) and -4.05 deg for the 4415
(published -4.0); fig 40 follows Prandtl-Glauert to within 2.6 % between M=0.30
and M=0.60, peaks at M=0.680, and its envelope contains all seven points of
`naca0015_lift_slope_vs_mach.csv` below M=0.68 with no slack; fig 42's band is
0.0048-0.0086 below M=0.60.
