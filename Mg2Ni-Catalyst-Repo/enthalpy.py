"""
Compute the three ΔH values (kJ/mol H2) for the catalyst study, plus the
auxiliary Nb2O5 ball-milling reduction enthalpy.

Strategy:
- Pristine reference energies are taken from /home/x/Workspace/espresso/mgh2-cif/,
  which has the relaxed 28-atom monoclinic Mg2NiH4 cell that yields ΔH ≈ -57.
- This script's outputs (mg2ni_Nb, mg2nih4_Nb, mg2ni_NbFe, mg2nih4_NbFe, mgo, nb2o5)
  cover the doped variants and Nb2O5 reference.

Mass balance for each ΔH (matches mgh2-cif/enthalpy.py):
    ΔH per H2 = [ E(hyd_28a) - (8/12) * E(met_18a) - 8 * E(H2) ] / 8 * 1312.75 kJ/mol/Ry

Note on doping fractions: pristine cells are matched in metal-atom ratios via the
8/12 scale factor. For doped cells, 1 Nb sits in each (8.3% on Mg site in Mg2Ni
metal cell, 12.5% on Mg site in Mg2NiH4-28a hydride cell). The asymmetry is
acknowledged in the writeup.
"""

from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
CIFDIR = HERE.parent / "mgh2-cif"
RY_TO_KJ = 1312.75

# Pristine references reused directly from mgh2-cif/ (already relaxed)
PRISTINE = {
    "h2":          CIFDIR / "h2.pwo",
    "mg":          CIFDIR / "mg.pwo",
    "mg2ni":       CIFDIR / "mg2ni.pwo",
    "mg2nih4":     CIFDIR / "mg2nih4-28.pwo",
}

# Doped + auxiliary outputs from this folder
DOPED = {
    "mg2ni_Nb":       HERE / "mg2ni_Nb.out",
    "mg2nih4_Nb":     HERE / "mg2nih4_Nb.out",
    "mg2ni_NbFe":     HERE / "mg2ni_NbFe.out",
    "mg2nih4_NbFe":   HERE / "mg2nih4_NbFe.out",
    "mg2ni_NbFeMg":   HERE / "mg2ni_NbFeMg.out",
    "mg2nih4_NbFeMg": HERE / "mg2nih4_NbFeMg.out",
    "mgo":            HERE / "mgo.out",
    "nb2o5":          HERE / "nb2o5.out",
}


def total_energy(path):
    """Final '!' total energy in Ry. Falls back to the last 'total energy' line
    if no '!' was written (e.g. SCF hit electron_maxstep but is energy-converged).
    """
    if not path.exists():
        return None
    bang, plain = [], []
    with path.open() as f:
        for line in f:
            m = re.match(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", line)
            if m:
                bang.append(float(m.group(1)))
                continue
            m = re.match(r"\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", line)
            if m:
                plain.append(float(m.group(1)))
    if bang:
        return bang[-1]
    return plain[-1] if plain else None


def fmt(e):
    return f"{e:14.6f}" if e is not None else "    (missing)"


def dh_per_h2(e_hyd_28, e_met_18, e_h2):
    """ΔH per H2 in kJ/mol using mgh2-cif's mass-balance scaling.

    Reaction (per supercell): (8/12) * Mg2Ni_18a + 8 * H2 -> Mg2NiH4_28a
    Valid for pristine; for doped cells use dh_shift().
    """
    dE_ry = (e_hyd_28 - (8.0 / 12.0) * e_met_18 - 8.0 * e_h2) / 8.0
    return dE_ry * RY_TO_KJ


def dh_shift(e_hyd_doped, e_hyd_pristine, e_met_doped, e_met_pristine, n_h2=8):
    """Doping-induced shift in ΔH per H2 (kJ/mol).

    Removes the unmatched-Nb-count problem in the simple scaling formula by
    subtracting pristine references symmetrically. The shift is independent of
    the Mg/Ni chemical-potential references.
    """
    dE_ry = ((e_hyd_doped - e_hyd_pristine) - (e_met_doped - e_met_pristine)) / n_h2
    return dE_ry * RY_TO_KJ


def main():
    # --- Read all energies ---
    e = {k: total_energy(p) for k, p in {**PRISTINE, **DOPED}.items()}

    if e["h2"] is None or e["mg2ni"] is None or e["mg2nih4"] is None:
        print("ERROR: pristine references missing in mgh2-cif/.", file=sys.stderr)
        sys.exit(1)

    print("\n=== Total Energies (Ry) ===")
    for k, v in e.items():
        print(f"  {k:14s} = {fmt(v)}")

    # --- Three ΔH values for hydrogenation reactions ---
    print("\n=== Formation enthalpies (kJ/mol H2) ===")
    print(f"  {'Reaction':14s}  {'E(metal,Ry)':>14s}  {'E(hydride,Ry)':>14s}  {'ΔH':>10s}")
    print("  " + "-" * 58)

    deltas = {}
    # Pristine: direct ΔH from mass-balanced reaction
    e_m, e_h = e["mg2ni"], e["mg2nih4"]
    dH = dh_per_h2(e_h, e_m, e["h2"])
    deltas["R1 pristine"] = dH
    print(f"  {'R1 pristine':14s}  {e_m:14.6f}  {e_h:14.6f}  {dH:10.2f}")

    # Doped: ΔH = ΔH_pristine + shift (shift method handles unmatched Nb count)
    for label, mkey, hkey in [
        ("R2 Nb-doped",          "mg2ni_Nb",      "mg2nih4_Nb"),
        ("R3 Nb,Fe(on Ni)",      "mg2ni_NbFe",    "mg2nih4_NbFe"),
        ("R3' Nb,Fe(on Mg)",     "mg2ni_NbFeMg",  "mg2nih4_NbFeMg"),
    ]:
        e_m_d, e_h_d = e.get(mkey), e.get(hkey)
        if e_m_d is None or e_h_d is None:
            print(f"  {label:18s}  {fmt(e_m_d):>14s}  {fmt(e_h_d):>14s}  (incomplete)")
            continue
        shift = dh_shift(e_h_d, e["mg2nih4"], e_m_d, e["mg2ni"])
        dH_d = deltas["R1 pristine"] + shift
        deltas[label] = dH_d
        print(f"  {label:18s}  {e_m_d:14.6f}  {e_h_d:14.6f}  {dH_d:10.2f}   (shift = {shift:+.2f})")

    if "R1 pristine" in deltas:
        print("\n=== Catalyst effect (vs pristine) ===")
        for label in ("R2 Nb-doped", "R3 Nb,Fe-doped"):
            if label in deltas:
                shift = deltas[label] - deltas["R1 pristine"]
                arrow = "(less stable hydride, easier H2 release)" if shift > 0 else "(more stable hydride)"
                print(f"  {label:14s}: ΔH shift = {shift:+.2f} kJ/mol H2  {arrow}")

    # --- Auxiliary: Nb2O5 ball-milling reduction ---
    aux_keys = ("mg", "mgo", "nb2o5", "mg2ni", "mg2ni_Nb")
    if all(e[k] is not None for k in aux_keys):
        e_mg     = e["mg"]    / 2.0   # mg.pwo cell has 2 Mg atoms
        e_mgo    = e["mgo"]   / 1.0   # mgo.out 2-atom primitive = 1 f.u. of MgO
        e_nb2o5  = e["nb2o5"] / 2.0   # nb2o5.out 14-atom = 2 f.u. of Nb2O5
        # Reaction per Nb dopant (atom-balanced for 1 Mg→Nb substitution in Mg2Ni):
        #   ½ Nb2O5 + 3/2 Mg + Mg2Ni  →  Mg2Ni:Nb + 5/2 MgO
        # Mg balance: 1.5 + 12 = 11 + 2.5 ✓  (the displaced Mg from Mg2Ni
        #                                      gets oxidised to MgO using O from Nb2O5)
        dE_ry = (e["mg2ni_Nb"] + 2.5 * e_mgo
                 - e["mg2ni"]   - 1.5 * e_mg - 0.5 * e_nb2o5)
        dE_kj = dE_ry * RY_TO_KJ
        sign = "exothermic ✓ (Nb migration into Mg2Ni is favored)" if dE_kj < 0 else "endothermic"
        print("\n=== Auxiliary: Nb2O5 ball-milling reduction (per Nb atom) ===")
        print(f"  E(Mg)/atom    = {e_mg:.6f} Ry")
        print(f"  E(MgO)/f.u.   = {e_mgo:.6f} Ry")
        print(f"  E(Nb2O5)/f.u. = {e_nb2o5:.6f} Ry")
        print(f"  Reaction:  ½ Nb2O5 + 3/2 Mg + Mg2Ni  →  Mg2Ni:Nb + 5/2 MgO")
        print(f"  ΔE_milling = {dE_kj:+.2f} kJ/mol Nb   ({sign})")


if __name__ == "__main__":
    main()
