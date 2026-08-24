"""
00_coarsen_datafile.py
======================
Collapses a LigParGen-style LAMMPS data file where every atom has its own
unique type into a molecule-aware element-grouped data file.

Typing scheme:  (mol_id, element)  →  new integer type
  → C in mol 1 ≠ C in mol 2, but all C atoms *within* mol 1 share one type.
  → Bonds/angles/dihedrals are typed by (mol_id, elem-sequence).

Force field coefficients (Pair Coeffs, Bond Coeffs, Angle Coeffs,
Dihedral Coeffs, Improper Coeffs) are preserved: values are carried over
from the original file and written under the remapped new type indices.

Also writes a molecule-topology JSON that the TACF/RgACF scripts use to
identify chemically equivalent dihedrals / chain membership without needing
to re-parse the original bloated file.

Usage:
    python 00_coarsen_datafile.py <original.data> [output_prefix]

Outputs:
    <prefix>_coarse.data      – coarsened LAMMPS data file
    <prefix>_topology.json    – per-molecule atom/bond/angle/dihedral lists
                                 + element-sequence dihedral pattern map
"""

import sys
import json
import math
from collections import defaultdict

# ---------------------------------------------------------------------------
# known masses → element symbol (extend as needed)
# ---------------------------------------------------------------------------
MASS_TO_ELEM = {
    1.008:   "H",
    12.011:  "C",
    14.007:  "N",
    15.999:  "O",
    18.998:  "F",
    32.067:  "S",
    35.453:  "Cl",
    6.941:   "Li",
    6.940:   "Li",
    22.990:  "Na",
    39.098:  "K",
}

def nearest_elem(mass):
    best = min(MASS_TO_ELEM.keys(), key=lambda m: abs(m - mass))
    if abs(best - mass) < 0.05:
        return MASS_TO_ELEM[best]
    return f"M{mass:.3f}"   # fallback: mass-labelled pseudoelement

# ---------------------------------------------------------------------------
# Parse helpers
# ---------------------------------------------------------------------------
def find_section(lines, name):
    """Return index of the blank-then-data line after 'name' header."""
    for i, l in enumerate(lines):
        if l.strip() == name or l.strip().startswith(name + " "):
            return i
    return None

def parse_masses(lines):
    idx = find_section(lines, "Masses")
    masses = {}
    if idx is None:
        return masses
    i = idx + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            break
        p = s.split()
        if len(p) >= 2:
            masses[int(p[0])] = float(p[1])
        i += 1
    return masses

def parse_atoms_full(lines):
    """Returns list of dicts with keys: id, mol, type, q, x, y, z, ix, iy, iz"""
    idx = find_section(lines, "Atoms")
    atoms = []
    if idx is None:
        return atoms
    i = idx + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            if atoms:
                break
            i += 1
            continue
        p = s.split("#", 1)[0].split()   # strip trailing "# comment" before column count
        if len(p) >= 7:
            d = dict(id=int(p[0]), mol=int(p[1]), type=int(p[2]),
                     q=float(p[3]), x=float(p[4]), y=float(p[5]), z=float(p[6]))
            if len(p) >= 10:
                d.update(ix=int(p[7]), iy=int(p[8]), iz=int(p[9]))
            else:
                d.update(ix=0, iy=0, iz=0)
            atoms.append(d)
        i += 1
    return atoms

def parse_bonds(lines):
    idx = find_section(lines, "Bonds")
    bonds = []
    if idx is None:
        return bonds
    i = idx + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            if bonds: break
            i += 1; continue
        p = s.split()
        if len(p) >= 4:
            bonds.append((int(p[0]), int(p[1]), int(p[2]), int(p[3])))
        i += 1
    return bonds

def parse_angles(lines):
    idx = find_section(lines, "Angles")
    angles = []
    if idx is None:
        return angles
    i = idx + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            if angles: break
            i += 1; continue
        p = s.split()
        if len(p) >= 5:
            angles.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4])))
        i += 1
    return angles

def parse_dihedrals(lines):
    idx = find_section(lines, "Dihedrals")
    dihedrals = []
    if idx is None:
        return dihedrals
    i = idx + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            if dihedrals: break
            i += 1; continue
        p = s.split()
        if len(p) >= 6:
            dihedrals.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]),
                               int(p[4]), int(p[5])))
        i += 1
    return dihedrals

# ---------------------------------------------------------------------------
# Coeffs parsers — each returns {old_type: (val, ...)}
# ---------------------------------------------------------------------------
def _parse_coeffs_generic(lines, section_name):
    """Return {type_int: tuple_of_floats} for any *Coeffs section."""
    idx = find_section(lines, section_name)
    coeffs = {}
    if idx is None:
        return coeffs
    i = idx + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            if coeffs:
                break
            i += 1
            continue
        # stop at next section header (no leading digit on a header line)
        p = s.split()
        try:
            t = int(p[0])
        except (ValueError, IndexError):
            break
        coeffs[t] = tuple(float(x) for x in p[1:] if _is_float(x))
        i += 1
    return coeffs

def _is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def parse_pair_coeffs(lines):
    return _parse_coeffs_generic(lines, "Pair Coeffs")

def parse_bond_coeffs(lines):
    return _parse_coeffs_generic(lines, "Bond Coeffs")

def parse_angle_coeffs(lines):
    return _parse_coeffs_generic(lines, "Angle Coeffs")

def parse_dihedral_coeffs(lines):
    return _parse_coeffs_generic(lines, "Dihedral Coeffs")

def parse_improper_coeffs(lines):
    return _parse_coeffs_generic(lines, "Improper Coeffs")

# ---------------------------------------------------------------------------
# Self-check: coarsened Coeffs must faithfully represent the source data
# ---------------------------------------------------------------------------
def verify_coarsening(atoms, bonds, angles, dihedrals,
                       atom_id_to_newtype,
                       pair_coeffs, new_pair_coeffs,
                       bond_key_to_newtype, bond_mol_seq, bond_coeffs, new_bond_coeffs,
                       angle_key_to_newtype, angle_mol_seq, angle_coeffs, new_angle_coeffs,
                       dih_key_to_newtype, dih_mol_seq, dihedral_coeffs, new_dihedral_coeffs):
    """
    Coarsening only relabels types -- it must never average or silently
    pick one of several genuinely different force-field values when many
    old types collapse into one new type. For every level (atom/bond/
    angle/dihedral) this checks:
      (a) every old type that had a Coeffs entry maps to a new type that
          also has one (catches dropped/overwritten entries), and
      (b) every old type collapsed into the same new type carried the
          IDENTICAL coeff values (catches silent loss of real variation).

    Raises RuntimeError on any violation. Always prints a summary.
    """
    errors = []

    def check_level(label, old_new_pairs, old_coeffs, new_coeffs):
        groups   = defaultdict(set)
        had_coeff = {}
        for old_t, new_t in old_new_pairs:
            if old_t in old_coeffs:
                groups[new_t].add(old_coeffs[old_t])
                had_coeff[new_t] = True
            else:
                had_coeff.setdefault(new_t, False)
        for new_t, vals in groups.items():
            if len(vals) > 1:
                errors.append(
                    f"{label} type {new_t}: {len(vals)} distinct coeff values "
                    f"collapsed together: {sorted(vals)}")
        for new_t, had in had_coeff.items():
            if had and new_t not in new_coeffs:
                errors.append(
                    f"{label} type {new_t}: source coeff existed but is "
                    f"missing from the coarsened Coeffs table")

    check_level(
        "Pair",
        ((a["type"], atom_id_to_newtype[a["id"]]) for a in atoms),
        pair_coeffs, new_pair_coeffs)

    check_level(
        "Bond",
        ((b_type, bond_key_to_newtype[bond_mol_seq[(b_id, b_type, i, j)]])
         for b_id, b_type, i, j in bonds),
        bond_coeffs, new_bond_coeffs)

    check_level(
        "Angle",
        ((ag_type, angle_key_to_newtype[angle_mol_seq[(ag_id, ag_type, i, j, k)]])
         for ag_id, ag_type, i, j, k in angles),
        angle_coeffs, new_angle_coeffs)

    check_level(
        "Dihedral",
        ((d_type, dih_key_to_newtype[dih_mol_seq[(d_id, d_type, i, j, k, l)]])
         for d_id, d_type, i, j, k, l in dihedrals),
        dihedral_coeffs, new_dihedral_coeffs)

    print(f"\n[verify] checked {len(atoms)} atoms / {len(bonds)} bonds / "
          f"{len(angles)} angles / {len(dihedrals)} dihedrals against "
          f"coarsened Coeffs")
    if errors:
        print(f"[verify] FAIL -- {len(errors)} issue(s):")
        for e in errors[:40]:
            print(f"    {e}")
        if len(errors) > 40:
            print(f"    ... and {len(errors) - 40} more")
        raise RuntimeError(
            f"coarsening self-check failed with {len(errors)} issue(s) -- "
            f"coarsened Coeffs do not faithfully represent the source data; "
            f"see printed detail above.")
    print("[verify] PASS -- every collapsed type carries a single, "
          "consistent coeff value; no coeffs dropped.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def coarsen(datafile, prefix=None):
    """
    Coarsen a LigParGen-style LAMMPS data file: collapse per-atom types into
    (mol_id, element) types, remap all connectivity and *Coeffs sections to
    new type indices, write <prefix>_coarse.data and <prefix>_topology.json.

    Args:
        datafile: path to original LAMMPS data file
        prefix:   output prefix (default: datafile with .data stripped)

    Returns:
        (coarse_data_path, topology_json_path)
    """
    if prefix is None:
        prefix = datafile.replace(".data", "")

    print(f"Reading {datafile} ...")
    with open(datafile) as f:
        lines = f.readlines()

    # ---- Parse geometry ----
    masses    = parse_masses(lines)
    atoms     = parse_atoms_full(lines)
    bonds     = parse_bonds(lines)
    angles    = parse_angles(lines)
    dihedrals = parse_dihedrals(lines)

    print(f"  {len(atoms)} atoms, {len(bonds)} bonds, "
          f"{len(angles)} angles, {len(dihedrals)} dihedrals")
    print(f"  {len(masses)} original atom types")

    # ---- Parse coeffs (values unchanged, indices will be remapped) ----
    pair_coeffs     = parse_pair_coeffs(lines)
    bond_coeffs     = parse_bond_coeffs(lines)
    angle_coeffs    = parse_angle_coeffs(lines)
    dihedral_coeffs = parse_dihedral_coeffs(lines)
    improper_coeffs = parse_improper_coeffs(lines)

    # ---- Map old type → element ----
    type_to_elem = {t: nearest_elem(m) for t, m in masses.items()}

    # ---- Per-atom maps ----
    atom_id_to_elem = {a["id"]: type_to_elem[a["type"]] for a in atoms}
    atom_mol        = {a["id"]: a["mol"] for a in atoms}

    # ---- New atom type numbering: key = (mol_id, element) ----
    mol_elem_pairs     = sorted(set((a["mol"], type_to_elem[a["type"]]) for a in atoms))
    mol_elem_to_newtype = {pair: n+1 for n, pair in enumerate(mol_elem_pairs)}
    atom_id_to_newtype = {
        a["id"]: mol_elem_to_newtype[(a["mol"], type_to_elem[a["type"]])]
        for a in atoms
    }

    # (old atom type, mol) → new atom type (for Pair Coeffs remapping)
    # Keyed by (old_type, mol), not old_type alone: LigParGen combine can
    # assign the SAME old type to multiple mols (e.g. a shared placeholder
    # type across repeated ion instances), which would otherwise collide
    # in a plain old_type→new_type dict and silently drop entries.
    old_atype_mol_to_new = {}
    for a in atoms:
        old_atype_mol_to_new[(a["type"], a["mol"])] = atom_id_to_newtype[a["id"]]

    n_atom_types = len(mol_elem_pairs)
    print(f"  Collapsed to {n_atom_types} mol-element types "
          f"(was {len(masses)} unique types)")
    print("  Atom type map  (mol, elem) → new_type:")
    for (mol, elem), ntype in sorted(mol_elem_to_newtype.items()):
        print(f"    mol {mol:2d}  {elem:2s}  →  type {ntype}")

    # ---- Per-molecule topology ----
    mol_atoms     = defaultdict(list)
    mol_bonds     = defaultdict(list)
    mol_angles    = defaultdict(list)
    mol_dihedrals = defaultdict(list)
    for a in atoms:
        mol_atoms[a["mol"]].append(a["id"])
    for b_id, b_type, i, j in bonds:
        mol_bonds[atom_mol[i]].append((b_id, i, j))
    for ag_id, ag_type, i, j, k in angles:
        mol_angles[atom_mol[i]].append((ag_id, i, j, k))
    for d_id, d_type, i, j, k, l in dihedrals:
        mol_dihedrals[atom_mol[i]].append((d_id, i, j, k, l))

    # ---- Dihedral element-sequence pattern map ----
    pattern_dihedrals = defaultdict(list)
    for d_id, d_type, i, j, k, l in dihedrals:
        seq = tuple(atom_id_to_elem[x] for x in [i, j, k, l])
        pattern_dihedrals[str(seq)].append([atom_mol[i], i, j, k, l])

    print("  Dihedral element patterns:")
    for pat, entries in sorted(pattern_dihedrals.items(), key=lambda x: -len(x[1])):
        mols = sorted(set(e[0] for e in entries))
        print(f"    {pat}: {len(entries)} dihedrals, mols {mols}")

    # ---- Topology JSON ----
    topo = {
        "atom_id_to_elem":  atom_id_to_elem,
        "atom_id_to_mol":   {str(k): v for k, v in atom_mol.items()},
        "mol_atom_ids":     {str(k): v for k, v in mol_atoms.items()},
        "mol_dihedral_ids": {str(k): [[d[0], d[1], d[2], d[3]] for d in v]
                             for k, v in mol_dihedrals.items()},
        "pattern_dihedrals": dict(pattern_dihedrals),
    }
    topo_path = prefix + "_topology.json"
    with open(topo_path, "w") as f:
        json.dump(topo, f)
    print(f"Saved topology → {topo_path}")

    # ---- New bond types: key = (mol_id, canonical-elem-pair) ----
    bond_mol_seq = {}
    # Also track old_btype → new_btype for Bond Coeffs remapping.
    # One old bond type may appear in multiple mols; each maps to a different
    # new type.  Build mapping per (old_btype, mol) pair instead.
    old_btype_mol_to_new = {}
    for b_id, b_type, i, j in bonds:
        mol = atom_mol[i]
        seq = tuple(sorted([atom_id_to_elem[i], atom_id_to_elem[j]]))
        bond_mol_seq[(b_id, b_type, i, j)] = (mol, seq)
    unique_bond_keys    = sorted(set(bond_mol_seq.values()))
    bond_key_to_newtype = {k: n+1 for n, k in enumerate(unique_bond_keys)}
    for b_id, b_type, i, j in bonds:
        key = bond_mol_seq[(b_id, b_type, i, j)]
        old_btype_mol_to_new[(b_type, atom_mol[i])] = bond_key_to_newtype[key]

    # ---- New angle types: key = (mol_id, canonical-elem-triple) ----
    angle_mol_seq = {}
    old_agtype_mol_to_new = {}
    for ag_id, ag_type, i, j, k in angles:
        mol = atom_mol[i]
        seq = (atom_id_to_elem[i], atom_id_to_elem[j], atom_id_to_elem[k])
        if seq[0] > seq[2]:
            seq = seq[::-1]
        angle_mol_seq[(ag_id, ag_type, i, j, k)] = (mol, seq)
    unique_angle_keys    = sorted(set(angle_mol_seq.values()))
    angle_key_to_newtype = {k: n+1 for n, k in enumerate(unique_angle_keys)}
    for ag_id, ag_type, i, j, k in angles:
        key = angle_mol_seq[(ag_id, ag_type, i, j, k)]
        old_agtype_mol_to_new[(ag_type, atom_mol[i])] = angle_key_to_newtype[key]

    # ---- New dihedral types: key = (mol_id, canonical-elem-quad) ----
    dih_mol_seq = {}
    old_dtype_mol_to_new = {}
    for d_id, d_type, i, j, k, l in dihedrals:
        mol = atom_mol[i]
        seq = tuple(atom_id_to_elem[x] for x in [i, j, k, l])
        if seq[0] > seq[3]:
            seq = seq[::-1]
        dih_mol_seq[(d_id, d_type, i, j, k, l)] = (mol, seq)
    unique_dih_keys    = sorted(set(dih_mol_seq.values()))
    dih_key_to_newtype = {k: n+1 for n, k in enumerate(unique_dih_keys)}
    for d_id, d_type, i, j, k, l in dihedrals:
        key = dih_mol_seq[(d_id, d_type, i, j, k, l)]
        old_dtype_mol_to_new[(d_type, atom_mol[i])] = dih_key_to_newtype[key]

    n_bond_types  = len(unique_bond_keys)
    n_angle_types = len(unique_angle_keys)
    n_dih_types   = len(unique_dih_keys)

    # ---- Remap Pair Coeffs: (old_atype, mol) → new_atype ----
    new_pair_coeffs = {}
    for (old_t, mol), new_t in old_atype_mol_to_new.items():
        if old_t in pair_coeffs:
            new_pair_coeffs[new_t] = pair_coeffs[old_t]

    # ---- Remap Bond Coeffs ----
    # old bond type can appear in multiple mols → multiple new types with same values.
    new_bond_coeffs = {}
    for (old_bt, mol), new_bt in old_btype_mol_to_new.items():
        if old_bt in bond_coeffs:
            new_bond_coeffs[new_bt] = bond_coeffs[old_bt]

    # ---- Remap Angle Coeffs ----
    new_angle_coeffs = {}
    for (old_agt, mol), new_agt in old_agtype_mol_to_new.items():
        if old_agt in angle_coeffs:
            new_angle_coeffs[new_agt] = angle_coeffs[old_agt]

    # ---- Remap Dihedral Coeffs ----
    new_dihedral_coeffs = {}
    for (old_dt, mol), new_dt in old_dtype_mol_to_new.items():
        if old_dt in dihedral_coeffs:
            new_dihedral_coeffs[new_dt] = dihedral_coeffs[old_dt]

    # ---- Self-check before writing anything ----
    verify_coarsening(atoms, bonds, angles, dihedrals,
                       atom_id_to_newtype,
                       pair_coeffs, new_pair_coeffs,
                       bond_key_to_newtype, bond_mol_seq, bond_coeffs, new_bond_coeffs,
                       angle_key_to_newtype, angle_mol_seq, angle_coeffs, new_angle_coeffs,
                       dih_key_to_newtype, dih_mol_seq, dihedral_coeffs, new_dihedral_coeffs)

    # ---- Extract box ----
    xlo_xhi = ylo_yhi = zlo_zhi = None
    for l in lines:
        if "xlo xhi" in l:
            p = l.split(); xlo_xhi = (float(p[0]), float(p[1]))
        if "ylo yhi" in l:
            p = l.split(); ylo_yhi = (float(p[0]), float(p[1]))
        if "zlo zhi" in l:
            p = l.split(); zlo_zhi = (float(p[0]), float(p[1]))

    # ---- Write coarsened data file ----
    out_path = prefix + "_coarse.data"
    with open(out_path, "w") as f:
        f.write("LAMMPS coarsened data file (element-grouped types)\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{n_atom_types} atom types\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{n_bond_types} bond types\n")
        f.write(f"{len(angles)} angles\n")
        f.write(f"{n_angle_types} angle types\n")
        f.write(f"{len(dihedrals)} dihedrals\n")
        f.write(f"{n_dih_types} dihedral types\n\n")
        f.write(f"{xlo_xhi[0]:.6f} {xlo_xhi[1]:.6f} xlo xhi\n")
        f.write(f"{ylo_yhi[0]:.6f} {ylo_yhi[1]:.6f} ylo yhi\n")
        f.write(f"{zlo_zhi[0]:.6f} {zlo_zhi[1]:.6f} zlo zhi\n\n")

        # Masses
        f.write("Masses\n\n")
        for (mol, elem), ntype in sorted(mol_elem_to_newtype.items(),
                                         key=lambda x: x[1]):
            mass = next(m for t, m in masses.items() if type_to_elem[t] == elem)
            f.write(f"{ntype} {mass:.4f}  # mol{mol}-{elem}\n")

        # Pair Coeffs
        if new_pair_coeffs:
            f.write("\nPair Coeffs\n\n")
            for ntype in sorted(new_pair_coeffs):
                vals = "  ".join(f"{v:.6f}" for v in new_pair_coeffs[ntype])
                f.write(f"{ntype:8d}  {vals}\n")

        # Bond Coeffs
        if new_bond_coeffs:
            f.write("\nBond Coeffs\n\n")
            for ntype in sorted(new_bond_coeffs):
                vals = "  ".join(f"{v:.4f}" for v in new_bond_coeffs[ntype])
                f.write(f"{ntype:8d}  {vals}\n")

        # Angle Coeffs
        if new_angle_coeffs:
            f.write("\nAngle Coeffs\n\n")
            for ntype in sorted(new_angle_coeffs):
                vals = "  ".join(f"{v:.4f}" for v in new_angle_coeffs[ntype])
                f.write(f"{ntype:8d}  {vals}\n")

        # Dihedral Coeffs
        if new_dihedral_coeffs:
            f.write("\nDihedral Coeffs\n\n")
            for ntype in sorted(new_dihedral_coeffs):
                vals = "  ".join(f"{v:.4f}" for v in new_dihedral_coeffs[ntype])
                f.write(f"{ntype:8d}  {vals}\n")

        # Improper Coeffs — no type remapping needed here (impropers share
        # the same old type numbering as dihedrals in OPLS-AA; if present,
        # write as-is since combine_lammps doesn't retype them).
        if improper_coeffs:
            f.write("\nImproper Coeffs\n\n")
            for old_t in sorted(improper_coeffs):
                vals = "  ".join(f"{v:.4f}" for v in improper_coeffs[old_t])
                f.write(f"{old_t:8d}  {vals}\n")

        # Atoms
        f.write("\nAtoms # full\n\n")
        for a in sorted(atoms, key=lambda x: x["id"]):
            ntype = atom_id_to_newtype[a["id"]]
            f.write(f"{a['id']:6d} {a['mol']:4d} {ntype:3d} "
                    f"{a['q']:12.6f} {a['x']:15.8f} {a['y']:15.8f} {a['z']:15.8f} "
                    f"{a['ix']:2d} {a['iy']:2d} {a['iz']:2d}\n")

        # Bonds
        f.write("\nBonds\n\n")
        for b_id, b_type, i, j in bonds:
            ntype = bond_key_to_newtype[bond_mol_seq[(b_id, b_type, i, j)]]
            f.write(f"{b_id:6d} {ntype:3d} {i:6d} {j:6d}\n")

        # Angles
        f.write("\nAngles\n\n")
        for ag_id, ag_type, i, j, k in angles:
            ntype = angle_key_to_newtype[angle_mol_seq[(ag_id, ag_type, i, j, k)]]
            f.write(f"{ag_id:6d} {ntype:3d} {i:6d} {j:6d} {k:6d}\n")

        # Dihedrals
        f.write("\nDihedrals\n\n")
        for d_id, d_type, i, j, k, l in dihedrals:
            ntype = dih_key_to_newtype[dih_mol_seq[(d_id, d_type, i, j, k, l)]]
            f.write(f"{d_id:6d} {ntype:3d} {i:6d} {j:6d} {k:6d} {l:6d}\n")

    print(f"Saved coarsened data → {out_path}")
    print(f"\nSummary of type collapse:")
    print(f"  Atom types:    {len(masses)} → {n_atom_types}  (mol × element)")
    print(f"  Bond types:    (original) → {n_bond_types}  (mol × elem-pair)")
    print(f"  Angle types:   (original) → {n_angle_types}  (mol × elem-triple)")
    print(f"  Dihedral types:(original) → {n_dih_types}  (mol × elem-quad)")
    print(f"\nNew dihedral type map:")
    for (mol, seq), ntype in sorted(dih_key_to_newtype.items(), key=lambda x: x[1]):
        print(f"  type {ntype:3d}: mol {mol:2d}  {'-'.join(seq)}")

    return out_path, topo_path


def main():
    if len(sys.argv) < 2:
        print("Usage: python 00_coarsen_datafile.py <datafile> [output_prefix]")
        sys.exit(1)
    datafile = sys.argv[1]
    prefix   = sys.argv[2] if len(sys.argv) > 2 else None
    coarsen(datafile, prefix)

if __name__ == "__main__":
    main()
