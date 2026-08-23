"""
mol2 reading for the mol2-native pipeline entry point.

mol2 already carries everything PSF+PDB carries between the two of them --
coordinates, elements (via SYBYL atom type), connectivity, and partial
charges -- in one file, and (for a packed system with several
@<TRIPOS>MOLECULE blocks concatenated back to back, as most packing tools
emit) it can represent a whole multi-molecule system directly. No PSF/PDB
synthesis, no obabel round-trip: species detection and coordinates are read
straight out of mol2's own sections.
"""
import re

# SYBYL atom type -> element. Dot-suffixed hybridization/geometry variants
# (C.3, C.2, C.ar, N.pl3, O.co2, S.O2, Cr.th, Cr.oh, Co.oh, ...) all collapse
# to the element before the dot. A handful of non-elemental "types" exist in
# the SYBYL spec (dummy atoms, lone pairs, wildcards) and are rejected below
# rather than silently mapped to something wrong.
_NOT_REAL_ATOMS = {'Du', 'Du.C', 'LP', 'Any', 'Hal', 'Het', 'Hev'}

# A few SYBYL types don't already look like a bare element after the dot-split
# (case or naming quirks); everything else round-trips through
# ``type.split('.')[0]`` correctly (C.3 -> C, N.ar -> N, Cl -> Cl, Na -> Na).
_SYBYL_ELEMENT_OVERRIDES = {
    'S.O': 'S', 'S.O2': 'S',
}


def sybyl_to_element(sybyl_type):
    if sybyl_type in _NOT_REAL_ATOMS:
        raise ValueError(
            f"SYBYL atom type {sybyl_type!r} is not a real atom (dummy/lone "
            f"pair/wildcard) -- mol2 input must be a concrete atomistic "
            f"structure, not a query/template mol2.")
    if sybyl_type in _SYBYL_ELEMENT_OVERRIDES:
        return _SYBYL_ELEMENT_OVERRIDES[sybyl_type]
    element = sybyl_type.split('.')[0]
    return element[0].upper() + element[1:].lower() if len(element) > 1 else element.upper()


_BOND_ORDER_TO_INT = {'1': 1, '2': 2, '3': 3, 'am': 1, 'ar': 1, 'du': 1, 'un': 1, 'nc': 0}


def parse_mol2(mol2_path):
    """Parse a mol2 file, possibly containing multiple concatenated
    @<TRIPOS>MOLECULE blocks (the usual way packing tools emit a packed
    multi-molecule system as mol2).

    Returns (natom, bonds, elements, coords, charges):
      natom   : total atom count across the whole file, globally renumbered
                1..natom in file order (block order, then atom order
                within each block)
      bonds   : [(a, b), ...] global 1-based atom-index pairs; bonds are
                intra-block only (mol2 has no cross-molecule bond records)
      elements: {serial: 'C'/'H'/'O'/...} from each atom's SYBYL type
      coords  : {serial: (x, y, z)} in Angstrom, as written in the file
      charges : {serial: float} per-atom partial charge column (0.0 if
                the file has no charge column for that block)
    """
    with open(mol2_path) as f:
        lines = f.readlines()

    blocks = []  # list of (atom_lines, bond_lines)
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip() == '@<TRIPOS>MOLECULE':
            i += 1
            # MOLECULE header: name, then counts line "natoms nbonds ..."
            i += 1  # skip mol name
            counts = lines[i].split()
            n_atoms_blk = int(counts[0])
            n_bonds_blk = int(counts[1]) if len(counts) > 1 else 0
            i += 1
            atom_lines = []
            bond_lines = []
            while i < n and lines[i].strip() != '@<TRIPOS>MOLECULE':
                if lines[i].strip() == '@<TRIPOS>ATOM':
                    i += 1
                    for _ in range(n_atoms_blk):
                        atom_lines.append(lines[i])
                        i += 1
                elif lines[i].strip() == '@<TRIPOS>BOND':
                    i += 1
                    for _ in range(n_bonds_blk):
                        bond_lines.append(lines[i])
                        i += 1
                else:
                    i += 1
            blocks.append((atom_lines, bond_lines))
        else:
            i += 1

    if not blocks:
        raise ValueError(f"No @<TRIPOS>MOLECULE block found in {mol2_path}")

    elements = {}
    coords = {}
    charges = {}
    bonds = []
    serial = 0
    for atom_lines, bond_lines in blocks:
        local_to_global = {}
        for line in atom_lines:
            parts = line.split()
            local_id = int(parts[0])
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            sybyl_type = parts[5]
            charge = float(parts[8]) if len(parts) > 8 else 0.0

            serial += 1
            local_to_global[local_id] = serial
            elements[serial] = sybyl_to_element(sybyl_type)
            coords[serial] = (x, y, z)
            charges[serial] = charge

        for line in bond_lines:
            parts = line.split()
            a_local, b_local = int(parts[1]), int(parts[2])
            bonds.append((local_to_global[a_local], local_to_global[b_local]))

    return serial, bonds, elements, coords, charges


def write_pdb_from_mol2(natom, elements, coords, out_pdb_path, resname='SYS'):
    """Write a plain PDB (single residue per atom-count-1000 chunk is not
    attempted here; residue/chain fields are cosmetic placeholders -- the
    rest of the pipeline only reads coordinates, element, and atom order
    from this file) from parsed mol2 data, for combine_lammps.run to source
    real system coordinates from.
    """
    with open(out_pdb_path, 'w') as f:
        for serial in range(1, natom + 1):
            x, y, z = coords[serial]
            el = elements[serial]
            name = el.upper().ljust(4) if len(el) <= 2 else el.upper()[:4]
            f.write(
                f"ATOM  {serial:5d} {name:<4s} {resname:<3s} A{1:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                f"{el.rjust(2)}\n"
            )
        f.write("END\n")
