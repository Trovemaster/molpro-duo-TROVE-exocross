6LiH MARVEL transition and segment files - 6 October 2026

TRANSITIONS: LiH6_MARVEL_transitions_audited.txt
SEGMENTS:    LiH6_MARVEL_segments.txt

2818 measured records: 2352 active and 466 excluded.
Active: 2333 primary X-X/A-X, one added microwave X-X, eight B-X, ten C-A.
The segment file contains all 11 source prefixes, each with unit cm-1.

MARVEL INPUT SETTINGS
Headerless ASCII, whitespace separated, with 12 columns:
wavenumber  initial_u  working_u  state' v' J' ef'  state" v" J" ef"  tag
Set FOUR quantum labels per level: electronic state, v, J, e/f.
Upper level comes first. X/A/C are singlet Sigma+ states with e parity.
B is singlet Pi: e for P/R and f for Q branches. Keep e and f distinct.
Wavenumbers and both uncertainty columns are in cm-1.
Negative wavenumbers retain excluded observations; do not change their signs.
This is the same input layout as the previously supplied 7LiH release.

AUDIT AND PROVENANCE
Source catalogue: LiH6_refinement/audit/transitions.json.
Primary workbook: LiH_MARVEL_MARVEL 2025-11.xlsx, sheet ^6LiH.
All 2799 primary rows are retained, including all 466 original sign exclusions.
The four PDF-verified frequency corrections at rows 125, 126, 151 and 207
are applied. Primary assignments and both uncertainty columns are preserved.
No new fit-based exclusions, calibration shifts or uncertainty rescaling were
introduced by this export. These are observed frequencies, not Duo predictions.

The additional Plummer (1984) X(0,1)-X(0,0) line is 453160.28 MHz converted
using 1 cm-1 = 29979.2458 MHz, with the audited 50 kHz uncertainty.
The eight Luh (1988) B-X observations retain the audited 0.1 nm wavelength
accuracy propagated to about 11.66-11.86 cm-1, rather than the unsupported
1 cm-1 value in the older workbook. This is not a formal one-sigma error.
The ten Chen (1999) C-A lines retain 0.1 cm-1 accuracy and the audited absolute
C vibrational labels 24,25,26,28,29; all upper J=5.

UNCERTAINTIES AND REMAINING LIMITS
Initial uncertainty means the supplied primary initial column where present;
it does not assert original published precision. For the additional 19 lines,
both columns use the already audited source accuracy/uncertainty scale.
The broad Iowa working uncertainties remain 0.8-1 cm-1. The previously noted
possible 0.23 cm-1 calibration correction was not applied without evidence.
The audit has not verified every original published line individually.

SOURCE TAGS
The 42 primary records sharing 21 reused 90MaOlTh tags receive fresh unused
numeric suffixes. The added microwave tag is normalized to 84PlHeDe.1.
All other tags are preserved. The line-by-line mapping JSON retains old tags,
workbook rows, original frequencies, corrections, PDF references and reasons.

VALIDATION
Every exported field was read back and checked against the curated catalogue.
All tags are unique, all source prefixes have exactly one segment entry, and
all assignments satisfy the electric-dipole J/e-f selection rules.
The active graph has 2 components and 658 levels;
656 levels are connected to X v=0, J=0.
Disconnected components have no determined zero relative to the ground level.
See LiH6_MARVEL_export_verification.json for their tags, counts and file hashes.
No MARVEL inversion or uncertainty optimization was run in this export.
Original workbooks, earlier exports, and Duo files remain unchanged.
