#!/usr/bin/env python3
"""make_trove_inputs.py: TROVE CH3OH inputs for a new PES (potential.chk).

Takes the working TROVE input (R-ALPHA132-THETA-S-TAU, C3v(M)-2, exact KEO
from kinetic.chk) and replaces EQUILIBRIUM, SPECPARAM and POTEN with the
values of the fitted PEF, so that the TROVE potential coordinates
(morse(x-xe) with chi_eq/specparam, x-xe) are identical to the fit
coordinates.  Writes J0.inp (BASIS_SET/CONTRACT/matelem save) and
J<n>.inp (read) files.
"""
import argparse
import re

TEMPLATE = ("/home/user/ExoMol/CH3OH/nb/"
            "CH3OH_PES_CH3OH_2d_alpha-123_C3v_202_pot_full_1609_KEO-10-11-12_pol8_harm_auto_step1.inp")


def replace_block(text, start_kw, new_block):
  """Replace block starting with a line == start_kw up to its 'end' line."""
  lines = text.split("\n")
  out, i = [], 0
  while i < len(lines):
    if lines[i].strip().upper() == start_kw.upper():
      j = i + 1
      while lines[j].strip().lower() != "end":
        j += 1
      out.extend(new_block.rstrip("\n").split("\n"))
      i = j + 1
      continue
    out.append(lines[i])
    i += 1
  return "\n".join(out)


def block(text, kw):
  m = re.search(rf"^{kw}\s*\n(.*?)^(end|END)\s*$", text, re.M | re.S)
  return f"{kw}\n{m.group(1)}{m.group(2)}"


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("--coords", required=True, help="potential.chk.trove_coords file")
  ap.add_argument("--poten", required=True, help=".trove file written by CH3OH_fitpes")
  ap.add_argument("--jmax", type=int, default=3)
  ap.add_argument("--template", default=TEMPLATE)
  ap.add_argument("--outdir", default=".")
  a = ap.parse_args()
  tmpl = open(a.template).read()
  crd = open(a.coords).read()
  eq = block(crd, "EQUILIBRIUM")
  sp = block(crd, "SPECPARAM")
  # POTEN: keep only the 8 non-linear parameters (documentation; the PES is
  # read from potential.chk, POTEN is not evaluated by TROVE in that case)
  pot_lines = [l for l in open(a.poten).read().split("\n")]
  i0 = [i for i, l in enumerate(pot_lines) if l.strip().upper() == "POTEN"][0]
  nl = [l for l in pot_lines[i0:] if l.split() and l.split()[0] in
        ("rCOe", "rOHe", "rCH", "aHOC", "aHCO", "a1", "a2", "b1")]
  poten = ("POTEN\nNPARAM  9\nPOT_TYPE  poten_ch3oh_sym_IV\nCOEFF  powers  (powers or list)\n"
           + "\n".join(nl) + "\n f   " + "   ".join(["0"] * 12) + "   1      0.0\nend")
  txt = replace_block(tmpl, "EQUILIBRIUM", eq)
  txt = replace_block(txt, "SPECPARAM", sp)
  txt = replace_block(txt, "POTEN", poten)
  for J in range(0, a.jmax + 1):
    t = re.sub(r"Jrot\s+\d+", f"Jrot  {J}", txt)
    mode = "save" if J == 0 else "read"
    t = re.sub(r"^BASIS_SET\s+\w+", f"BASIS_SET   {mode}", t, flags=re.M)
    t = re.sub(r"^CONTRACT\s+\w+", f"CONTRACT    {mode}", t, flags=re.M)
    t = re.sub(r"^matelem\s+\w+\s+split", f"matelem     {mode} split", t, flags=re.M)
    with open(f"{a.outdir}/J{J}.inp", "w") as f:
      f.write(t)
  print(f"wrote J0..J{a.jmax}.inp in {a.outdir}")


if __name__ == "__main__":
  main()
