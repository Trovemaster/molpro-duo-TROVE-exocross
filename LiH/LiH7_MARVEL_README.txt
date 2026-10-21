7LiH MARVEL transition input - 5 October 2026

MAIN FILE: LiH7_MARVEL_transitions_audited.txt
3731 measured records: 3241 active, 490 excluded/provisionally held.
Active: 3158 X-X/A-X, 73 B-X, 10 C-A.
Headerless ASCII, whitespace separated, with 12 fields per line:

wavenumber  initial_u  working_u  state' v' J' ef'  state" v" J" ef"  tag

All three numerical spectral quantities are in cm-1. Set the number of quantum
labels per state to FOUR: electronic state, v, J, e/f parity. Upper state comes
first. X, A and C are singlet Sigma+ states (e); B is singlet Pi, with e for P/R
and f for Q branches. Do not merge B e/f levels or run this as a 3-label input.
Negative wavenumbers retain excluded observations without activating them.
The companion LiH7_MARVEL_segments.txt declares cm-1 for every source prefix.

PROVENANCE
This uses the latest curated 7LiH catalogue in LiH_refinement/audit/transitions.json,
combining the primary November workbook audit with subsequently verified B/C
tables. These are observed frequencies, not values calculated in the Duo fit.
The original primary workbook was LiH_MARVEL_MARVEL 2025-11.xlsx.
All 3634 primary records, 87 Velasco B-X entries and 10 Chen C-A entries remain.
Nine PDF-verified primary frequency corrections are applied. The five corrected
Velasco Q(1-0) J assignments and the audited C absolute vibrational numbering
are retained. No unverified later-export assignment or blanket sign changes
have been adopted. No new uncertainty adjustment was made by this exporter.

UNCERTAINTIES
Both uncertainty columns are preserved for primary records. 'Initial' means
the primary workbook's initial column, not an assertion of original published
precision: some values were already inflated. For B/C, both columns start from
the audited published accuracy scales (0.04 and 0.1 cm-1 respectively). These
scales are not assumed formal one-standard-deviation uncertainties.

EXCLUSION POLICY
The latest conservative curation is retained: 363 primary-sign exclusions,
108 additional positive Rafi-1993 lines held following Coxon/Dickinson's source
assessment, and five individually unresolved primary rows (542, 598, 3366,
3376, 3501). Row 3366 is a provisional model-discrepancy hold, not a PDF-proven
transcription error; its measured frequency is unchanged. Row 598 retains the
primary assignment and remains negative pending the original source. The
unphysical Chan bridge at row 3376 is negative; no guessed repair was made.
Fourteen B entries remain negative because the original table flags them as
diffuse, blended, or doubtful. Exclusion reasons are in the mapping JSON.

MARVEL VERSUS DUO
Four otherwise accepted B(v=0,J=8) shape-resonance transitions are restored to
positive frequency: a Duo finite-box limitation does not invalidate an observed
transition for MARVEL. Three B(v=2,J=2) lines also remain active; their previous
box-convergence caution concerns calculated energies, not observed frequencies.
No further inclusion decision was inferred from a small Duo residual.

TAGS AND TRACEABILITY
All tags are unique and have numeric source.sequence suffixes. The four records
sharing two duplicated 90MaOlTh identifiers were given new unused numbers.
The 87 descriptive Velasco tags were converted to 57Velasco.1 through .87 in
audited table order. All original tags, workbook rows, PDF references, changes,
uncertainty bases and exclusion reasons are preserved line-by-line in
LiH7_MARVEL_audit_mapping.json.

CHECKS
The exported file was read back and every field checked against the catalogue:
isotope, signed frequency, both uncertainties, four labels per state, unique
tags, electric-dipole J/e-f selection rules, and all nine frequency corrections.
The active graph has 6 component(s), containing 906 levels.
No official MARVEL inversion or uncertainty optimisation was executed. This is
an audited input set for MARVEL, not a claim of newly validated MARVEL energies.
Input/output hashes and counts are in LiH7_MARVEL_export_verification.json.

Format reference: ReSpecTh data descriptor (2025), MARVEL transition input:
https://www.nature.com/articles/s41597-025-05272-6
Original workbooks, previous exports, and Duo results were not overwritten.
