# PtolemyCpp cross-checks of the a1954 14C FRESCO calculations

[PtolemyCpp](https://github.com/goluckyryan/PtolemyCpp) is a C++ reimplementation of the Argonne
Ptolemy code. It builds with nothing but `g++ -std=c++17` (`make`, then `make test` -> 128/128
here), which makes it a cheap independent check on the FRESCO decks in `../fresco/`.

```
git clone https://github.com/goluckyryan/PtolemyCpp.git && cd PtolemyCpp && make
PTOLEMY=/path/to/PtolemyCpp/ptolemy ./run_ptolemy.sh
```

`run_ptolemy.sh` runs every deck in `inputs/` and parses the angular distributions into `dat/`,
which is what the ROOT macros read. Only `dat/` is tracked; `outputs/` is regenerable.

## What it was used for

**1. Validating the FRESCO elastic deck.** Given byte-identical KD03 input, the two codes agree to
**0.05 %** over the whole angular range including through the diffraction minimum; the residual is
Ptolemy printing 4 significant figures. Ptolemy uses `R = r0 * A_target^(1/3)`, the same convention
as the FRESCO decks, so the potentials really are identical (`R0TARGET` is a no-op for these).
This matters because a wrong FRESCO deck had already cost a day once -- see the `p(1:5)` multipole
trap documented in `../fresco/README.md`.

**2. Validating the FRESCO inelastic decks -- but only after picking the right coupling.**
Ptolemy has three collective treatments and they are *not* interchangeable:

| coupling | rms ln(Ptolemy/FRESCO), 7.012 2+ , 25-135 deg |
|----------|-----------------------------------------------|
| INELOCA1 | **0.087**  (within 2 % everywhere forward of 105 deg) |
| INELOCA2 | 0.131 |
| INELOCA3 | 0.210 |

**Use INELOCA1.** INELOCA3 is flatter by a factor 1.7 and disagreeing with it means nothing.
The one genuine difference that survives is the backward tail: at 135 deg FRESCO is 29 % below
Ptolemy and still falling while Ptolemy flattens. Since the data run to 135 deg that is a real
model systematic on the backward end of every extracted deformation length.

Ptolemy's inelastic engine keys off `BELX` = B(EL), not a deformation length, so its normalisation
is not directly comparable to FRESCO's `&pot type=11` delta -- only the shape is. Set `LMAX=30`:
the default fails to converge (60 of 160 angles above 5 % high-L truncation error) and `LMAX >= 60`
falls into the extrapolation failure and returns all zeros. The compact `.dwba` one-liner format
expands but does not run to a cross section for inelastic; use the native deck.

**3. The optical-potential sweep.** This is what Ptolemy adds over FRESCO -- a built-in OMP library,
so a sweep is one line per potential. See `../pp/elastic_omp_C14.C` and `../pp/elastic_dip_C14.C`.

## The validity caveat that frames all of it

Every global proton OMP in the library is **outside its stated validity** for p + 14C at 11.58 MeV:

| potential | key | stated range | p+14C at 11.58 MeV |
|-----------|-----|--------------|--------------------|
| Koning-Delaroche (KD03) | `K` | 24 < A < 209   | A = 14, below |
| Varner (CH89)           | `V` | 4 < A < 209, 16 < E < 65 | E = 11.6, below |
| Becchetti-Greenlees     | `G` | 40 < A         | below |
| Perey                   | `P` | 30 < A < 100   | below |
| Menet                   | `M` | 40 < A, 30 < E < 60 | both below |

So the spread across them is a spread among five extrapolations, not a menu with a right answer
in it. Quote it as a systematic; do not "choose the best potential" from it.
