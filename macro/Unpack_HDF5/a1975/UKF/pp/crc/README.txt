CRC curves digitised from the theory report (report.pdf, "Results of 16C(p,p') and (d,d') Analysis",
dated 2025-09-05) and its figures gs.pdf / 0.74.pdf / 3.1.pdf.

These are COUPLED REACTION CHANNEL calculations, not the DWBA I built in 16C_dp/fresco. The report
fitted the <16C|15C+n> spectroscopic amplitudes TO the measured 16C(p,d)15C angular distributions,
so the curves are ABSOLUTE -- they must not be rescaled when overlaid on data. "Full" is the
complete coupling scheme; "Reduced" drops some couplings.

  gs_*.dat    15C g.s.   1/2+
  e074_*.dat  15C 0.740  5/2+
  e310_*.dat  15C 3.103  1/2-   (our fit puts this structure at 3.38, see the note below)

Digitised by digitize.py from the VECTOR content of the PDFs -- the curve coordinates are in the
file, so no pixel-picking was involved. Calibration comes from each plot's own axis furniture:
angles from the labels below the frame, the log scale from the 10^n decade labels (or, where only
one decade is labelled, from the 2..9 minor-tick ladder, which agreed to 0.21 %).

The report's own data points are ALSO in those figures and were NOT digitised here; only the curves
were. The report analyses the same 11.5 MeV/nucleon 16C data set, so its points are an independent
analysis of our experiment and would be worth extracting as a cross-check.

NOTE ON THE THIRD STATE: the report couples to the 15C 3.103 MeV 1/2- level; our spectrum fit puts
that structure at 3.380 +- 0.008. The 277 keV difference is not accounted for.
