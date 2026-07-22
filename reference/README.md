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

**Not extracted: full CL/CD/CM polars.** Three attempts at automated curve
tracing failed. On the multi-Mach sheets (1.3.20, 1.3.60) eight to twelve curves
overlap on one axis and cannot be separated at this scan quality; on the
single-curve 180° sheet (1.3.30) the negative-CL branch sits in a region of the
scan that is uniformly black. Anything presented as a dense digitised polar from
this document would be invented, so none is provided. For a precision benchmark,
go to the original sources the Datcom cites — NACA TN 3361 for the 180° data,
ARC C.P. 1261 (Gregory & Wilby) for the NPL NACA 0012.

## Using the landmarks

The 180° sheet is the useful one for a stall study: single Reynolds number
(1.8e6), fully turbulent, one curve per plot. Compare a computed sweep against
it with the caveat that steady RANS is not expected to reproduce post-stall
behaviour — see the stall notes in the top-level README.

```
python reference/check_reference.py    # verifies the coordinate tables
```
