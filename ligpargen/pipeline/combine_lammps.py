#!/usr/bin/env python3
"""
Step 2 of 2. Reassemble full LAMMPS data file from per-species ligpargen
.lmp outputs + manifest.json (from extract_unique_molecules.py) + the
original combined PDB (for real per-instance coordinates).

Species with n_atoms==1 (bare ions) have no ligpargen .lmp -- these get
a synthesized single-atom type with PLACEHOLDER charge/sigma/epsilon
(sentinel -9999, flagged with a trailing comment) so they're easy to
grep and fill in by hand from oplsaa.par.

Type numbering: each species contributes N_species distinct atom/bond/
angle/dihedral/improper TYPES (ligpargen gives every atom and every bond
its own type -- no reuse within a species). All instances of the same
species SHARE that species' types; only atom/bond/angle/... IDs are
unique per instance.

Usage: python combine_lammps.py manifest.json original.pdb species_dir/ out.lmp [box_from.data]
Looks for species_dir/species_N.lammps.lmp for each polyatomic species.

Box bounds: pass box_from (any LAMMPS data file) to use its real periodic
box. Each atom is then wrapped into it independently (standard LAMMPS/
VMD/OVITO convention: x - L*floor((x-lo)/L) per dimension, per atom) --
a molecule straddling the boundary WILL have atoms stored on opposite
sides, which is expected; bonded distances are still correct because
LAMMPS (and any PBC-aware viewer) uses the minimum image. Without
box_from, the box falls back to bbox(coords)+10A padding, which is NOT
physically meaningful for a PDB with molecules unwrapped across the
periodic boundary (common for EMC-packed systems) -- set it manually
in that case.

After wrapping, does an O(n^2) minimum-image overlap scan and warns on
any non-bonded atom pair closer than 1.0 Angstrom (real clash).
"""
import sys
import os
import json
import math

PLACEHOLDER = -9999.0

# full periodic table (Z=1-118), IUPAC standard atomic weights -- conventional
# value where an interval is given; mass number of longest-lived isotope for
# Tc, Pm, and Z>=84 (no stable isotope)
ATOMIC_MASS = {
    'H': 1.008, 'HE': 4.0026, 'LI': 6.94, 'BE': 9.0122, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'NE': 20.180,
    'NA': 22.990, 'MG': 24.305, 'AL': 26.982, 'SI': 28.085, 'P': 30.974,
    'S': 32.06, 'CL': 35.45, 'AR': 39.948, 'K': 39.098, 'CA': 40.078,
    'SC': 44.956, 'TI': 47.867, 'V': 50.942, 'CR': 51.996, 'MN': 54.938,
    'FE': 55.845, 'CO': 58.933, 'NI': 58.693, 'CU': 63.546, 'ZN': 65.38,
    'GA': 69.723, 'GE': 72.630, 'AS': 74.922, 'SE': 78.971, 'BR': 79.904,
    'KR': 83.798, 'RB': 85.468, 'SR': 87.62, 'Y': 88.906, 'ZR': 91.224,
    'NB': 92.906, 'MO': 95.95, 'TC': 97, 'RU': 101.07, 'RH': 102.91,
    'PD': 106.42, 'AG': 107.87, 'CD': 112.41, 'IN': 114.82, 'SN': 118.71,
    'SB': 121.76, 'TE': 127.60, 'I': 126.90, 'XE': 131.29, 'CS': 132.91,
    'BA': 137.33, 'LA': 138.91, 'CE': 140.12, 'PR': 140.91, 'ND': 144.24,
    'PM': 145, 'SM': 150.36, 'EU': 151.96, 'GD': 157.25, 'TB': 158.93,
    'DY': 162.50, 'HO': 164.93, 'ER': 167.26, 'TM': 168.93, 'YB': 173.05,
    'LU': 174.97, 'HF': 178.49, 'TA': 180.95, 'W': 183.84, 'RE': 186.21,
    'OS': 190.23, 'IR': 192.22, 'PT': 195.08, 'AU': 196.97, 'HG': 200.59,
    'TL': 204.38, 'PB': 207.2, 'BI': 208.98, 'PO': 209, 'AT': 210,
    'RN': 222, 'FR': 223, 'RA': 226, 'AC': 227, 'TH': 232.04,
    'PA': 231.04, 'U': 238.03, 'NP': 237, 'PU': 244, 'AM': 243,
    'CM': 247, 'BK': 247, 'CF': 251, 'ES': 252, 'FM': 257,
    'MD': 258, 'NO': 259, 'LR': 262, 'RF': 267, 'DB': 268,
    'SG': 271, 'BH': 272, 'HS': 270, 'MT': 276, 'DS': 281,
    'RG': 280, 'CN': 285, 'NH': 286, 'FL': 289, 'MC': 290,
    'LV': 293, 'TS': 294, 'OG': 294,
}


def parse_species_lmp(path):
    """Parse a ligpargen-written .lmp file. Format is fixed by
    ligpargen/inout/lammps.py -- this is not a generic LAMMPS parser."""
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f]

    def section(name):
        start = None
        for i, l in enumerate(lines):
            if l.strip() == name:
                start = i + 1
                break
        if start is None:
            return []
        i = start
        while i < len(lines) and lines[i].strip() == '':  # skip blank line(s) after header
            i += 1
        out = []
        while i < len(lines):
            row = lines[i].strip()
            if row == '' or not row.split()[0].lstrip('-').isdigit():
                break  # blank line or next section header -- empty section, stop
            out.append(lines[i].split())
            i += 1
        return out

    masses = [float(r[1]) for r in section('Masses')]
    pair = [(float(r[1]), float(r[2])) for r in section('Pair Coeffs')]  # (epsilon, sigma)
    bond_coeffs = [(float(r[1]), float(r[2])) for r in section('Bond Coeffs')]
    angle_coeffs = [(float(r[1]), float(r[2])) for r in section('Angle Coeffs')]
    dihedral_coeffs = [tuple(float(x) for x in r[1:5]) for r in section('Dihedral Coeffs')]
    improper_coeffs = [(float(r[1]), int(r[2]), float(r[3])) for r in section('Improper Coeffs')]

    atoms = section('Atoms')
    charges = [float(r[3]) for r in atoms]  # local atom order == index into masses/pair

    bonds = [(int(r[2]), int(r[3])) for r in section('Bonds')]
    angles = [(int(r[2]), int(r[3]), int(r[4])) for r in section('Angles')]
    dihedrals = [tuple(int(x) for x in r[2:6]) for r in section('Dihedrals')]
    impropers = [tuple(int(x) for x in r[2:6]) for r in section('Impropers')]

    n_atoms = len(masses)
    if not (len(pair) == len(charges) == n_atoms):
        raise ValueError(f"{path}: atom-type count mismatch "
                          f"(masses={len(masses)} pair={len(pair)} atoms={len(charges)})")

    return dict(masses=masses, pair=pair, charges=charges,
                bond_coeffs=bond_coeffs, angle_coeffs=angle_coeffs,
                dihedral_coeffs=dihedral_coeffs, improper_coeffs=improper_coeffs,
                bonds=bonds, angles=angles, dihedrals=dihedrals, impropers=impropers)


def extract_element(raw_token):
    """PDB element column isn't always a clean symbol (e.g. psfgen-style
    'saltLi', 'LI1'). Match the longest known symbol found in the token,
    case-insensitive, rather than requiring an exact match."""
    token = ''.join(c for c in raw_token.upper() if c.isalpha())
    for sym in sorted(ATOMIC_MASS, key=len, reverse=True):  # 2-letter before 1-letter
        if sym in token:
            return sym
    raise ValueError(f"can't find a known element in '{raw_token}' -- "
                      f"add the element to ATOMIC_MASS in combine_lammps.py")


def synth_ion_species(element):
    el = extract_element(element)
    return dict(masses=[ATOMIC_MASS[el]], pair=[(PLACEHOLDER, PLACEHOLDER)],
                charges=[PLACEHOLDER],
                bond_coeffs=[], angle_coeffs=[], dihedral_coeffs=[], improper_coeffs=[],
                bonds=[], angles=[], dihedrals=[], impropers=[])


def read_pdb_coords_and_elements(pdb_path):
    """1-indexed by ATOM/HETATM order -- matches extract_unique_molecules.py."""
    coords, elements, serial = {}, {}, 0
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                serial += 1
                coords[serial] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
                elements[serial] = l.split()[-1]
    return coords, elements


def parse_box_bounds(data_path):
    """Read xlo/xhi, ylo/yhi, zlo/zhi from any LAMMPS data file's header."""
    bounds = {}
    with open(data_path) as f:
        for l in f:
            row = l.split('#')[0].split()
            if len(row) >= 4 and row[2:4] in (['xlo', 'xhi'], ['ylo', 'yhi'], ['zlo', 'zhi']):
                bounds[row[2]] = float(row[0])
                bounds[row[3]] = float(row[1])
    missing = {'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi'} - bounds.keys()
    if missing:
        raise ValueError(f"{data_path}: missing box bounds {missing}")
    return bounds


def parse_box_from_cryst1(pdb_path):
    """Read box dimensions from the CRYST1 record of a PDB file.
    CRYST1 format (cols 1-80):
      CRYST1  a  b  c  alpha  beta  gamma  sGroup  Z
    Only orthogonal cells (alpha=beta=gamma=90) are supported here --
    non-orthogonal cells need a tilt-factor extension in the LAMMPS data file
    which combine_lammps.py doesn't currently write.
    Returns (xlo, xhi, ylo, yhi, zlo, zhi) with lo=0."""
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('CRYST1'):
                parts = l.split()
                a, b, c = float(parts[1]), float(parts[2]), float(parts[3])
                alpha, beta, gamma = float(parts[4]), float(parts[5]), float(parts[6])
                if abs(alpha - 90.0) > 0.01 or abs(beta - 90.0) > 0.01 or abs(gamma - 90.0) > 0.01:
                    raise ValueError(
                        f"CRYST1 in {pdb_path} has non-orthogonal cell "
                        f"(alpha={alpha}, beta={beta}, gamma={gamma}) -- "
                        "only orthogonal cells (90/90/90) are supported; "
                        "use box_from with a LAMMPS data file instead.")
                return (0.0, a, 0.0, b, 0.0, c)
    raise ValueError(f"No CRYST1 record found in {pdb_path}")


def wrap_coord(v, lo, hi):
    """Standard periodic wrap: map a single coordinate into [lo, hi).
    Applied per-atom, independently in each dimension -- this is the
    LAMMPS/VMD/OVITO convention. A molecule straddling the boundary WILL
    end up with atoms on opposite sides in the stored coordinates (looks
    "torn" in a naive plot) -- that's expected. Bonded interactions in
    LAMMPS (and any PBC/topology-aware viewer) use the minimum image, so
    the true bond length is still correct regardless of raw coordinates."""
    L = hi - lo
    return v - L * math.floor((v - lo) / L)


def find_overlaps(atom_ids, atom_coords, bonded_pairs, box, threshold=1.0):
    """Brute-force O(n^2) minimum-image overlap scan -- flags any
    non-bonded atom pair closer than `threshold` (real clash, not just a
    short bond). Fine into the low thousands of atoms; a cell list would
    be needed for much larger systems."""
    xlo, xhi, ylo, yhi, zlo, zhi = box
    Lx, Ly, Lz = xhi - xlo, yhi - ylo, zhi - zlo
    n = len(atom_ids)
    hits = []
    for i in range(n):
        xi, yi, zi = atom_coords[i]
        for j in range(i + 1, n):
            pair = frozenset((atom_ids[i], atom_ids[j]))
            if pair in bonded_pairs:
                continue
            xj, yj, zj = atom_coords[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj
            dx -= Lx * round(dx / Lx)
            dy -= Ly * round(dy / Ly)
            dz -= Lz * round(dz / Lz)
            d = (dx * dx + dy * dy + dz * dz) ** 0.5
            if d < threshold:
                hits.append((atom_ids[i], atom_ids[j], d))
    return hits


def run(manifest_path, pdb_path, species_dir, out_path, box_from=None,
         box_from_pdb=False):
    with open(manifest_path) as f:
        manifest = json.load(f)

    coords, elements = read_pdb_coords_and_elements(pdb_path)

    box = None
    if box_from_pdb:
        box = parse_box_from_cryst1(pdb_path)
        print(f"Box from CRYST1 in {pdb_path}: "
              f"x=[{box[0]:.4f},{box[1]:.4f}] "
              f"y=[{box[2]:.4f},{box[3]:.4f}] "
              f"z=[{box[4]:.4f},{box[5]:.4f}]")
    elif box_from is not None:
        b = parse_box_bounds(box_from)
        box = (b['xlo'], b['xhi'], b['ylo'], b['yhi'], b['zlo'], b['zhi'])

    species_data, placeholder_species = {}, []
    for sid_str, sp in manifest['species'].items():
        sid = int(sid_str)
        if sp['n_atoms'] == 1:
            rep_serial = sp['representative_atom_serials'][0]
            species_data[sid] = synth_ion_species(elements[rep_serial])
            placeholder_species.append(sid)
        else:
            lmp_path = os.path.join(species_dir, f"species_{sid}.lammps.lmp")
            if not os.path.exists(lmp_path):
                raise FileNotFoundError(
                    f"missing {lmp_path} -- run ligpargen on species_{sid}.pdb first")
            species_data[sid] = parse_species_lmp(lmp_path)

    # --- assign global type offsets, one pass, fixed species order ---
    sid_order = sorted(species_data)
    type_off = {}  # sid -> (atom, bond, angle, dihedral, improper) offsets
    off_a = off_b = off_ang = off_d = off_i = 0
    all_masses, all_pair = [], []
    all_bond_coeffs, all_angle_coeffs, all_dihedral_coeffs, all_improper_coeffs = [], [], [], []
    for sid in sid_order:
        d = species_data[sid]
        type_off[sid] = (off_a, off_b, off_ang, off_d, off_i)
        all_masses += d['masses']
        all_pair += d['pair']
        all_bond_coeffs += d['bond_coeffs']
        all_angle_coeffs += d['angle_coeffs']
        all_dihedral_coeffs += d['dihedral_coeffs']
        all_improper_coeffs += d['improper_coeffs']
        off_a += len(d['masses'])
        off_b += len(d['bond_coeffs'])
        off_ang += len(d['angle_coeffs'])
        off_d += len(d['dihedral_coeffs'])
        off_i += len(d['improper_coeffs'])

    # --- replicate per instance ---
    atom_id = bond_id = angle_id = dihedral_id = improper_id = 0
    atoms_out, bonds_out, angles_out, dihedrals_out, impropers_out = [], [], [], [], []
    xs, ys, zs = [], [], []
    bonded_pairs = set()

    for inst_id_str, inst in manifest['instances'].items():
        inst_id = int(inst_id_str)
        sid = inst['species_id']
        d = species_data[sid]
        a_off, b_off, ang_off, dih_off, imp_off = type_off[sid]
        serials = inst['atom_serials']  # real absolute serials; local order == template order

        raw_coords = [coords[s] for s in serials]
        if box is not None:
            xlo, xhi, ylo, yhi, zlo, zhi = box
            wrapped_coords = [(wrap_coord(x, xlo, xhi), wrap_coord(y, ylo, yhi), wrap_coord(z, zlo, zhi))
                               for x, y, z in raw_coords]
        else:
            wrapped_coords = raw_coords

        local_to_global = {}
        for local_idx, (real_serial, (x, y, z)) in enumerate(zip(serials, wrapped_coords), start=1):
            atom_id += 1
            local_to_global[local_idx] = atom_id
            xs.append(x); ys.append(y); zs.append(z)
            atype = a_off + local_idx
            charge = d['charges'][local_idx - 1]
            line = (f"{atom_id:8d} {inst_id:6d} {atype:6d} {charge:10.6f} "
                    f"{x:10.6f} {y:10.6f} {z:10.6f}")
            if charge == PLACEHOLDER:
                line += "  # PLACEHOLDER ion charge -- fill from oplsaa.par"
            atoms_out.append(line)

        for local_bidx, (a, b) in enumerate(d['bonds'], start=1):
            bond_id += 1
            ga, gb = local_to_global[a], local_to_global[b]
            bonded_pairs.add(frozenset((ga, gb)))
            bonds_out.append(f"{bond_id:6d} {b_off+local_bidx:6d} {ga:6d} {gb:6d}")

        for local_aidx, (a, b, c) in enumerate(d['angles'], start=1):
            angle_id += 1
            angles_out.append(f"{angle_id:6d} {ang_off+local_aidx:6d} "
                               f"{local_to_global[a]:6d} {local_to_global[b]:6d} {local_to_global[c]:6d}")

        for local_didx, (a, b, c, e) in enumerate(d['dihedrals'], start=1):
            dihedral_id += 1
            dihedrals_out.append(f"{dihedral_id:6d} {dih_off+local_didx:6d} "
                                  f"{local_to_global[a]:6d} {local_to_global[b]:6d} "
                                  f"{local_to_global[c]:6d} {local_to_global[e]:6d}")

        for local_iidx, (a, b, c, e) in enumerate(d['impropers'], start=1):
            improper_id += 1
            impropers_out.append(f"{improper_id:6d} {imp_off+local_iidx:6d} "
                                  f"{local_to_global[a]:6d} {local_to_global[b]:6d} "
                                  f"{local_to_global[c]:6d} {local_to_global[e]:6d}")

    # --- finalize box + overlap check ---
    if box is not None:
        xlo, xhi, ylo, yhi, zlo, zhi = box
        # sanity check only -- wrap_coord guarantees every atom lands in
        # [lo, hi) by construction, so this should never actually trip
        assert min(xs) >= xlo and max(xs) < xhi, "x coord escaped the box after wrapping"
        assert min(ys) >= ylo and max(ys) < yhi, "y coord escaped the box after wrapping"
        assert min(zs) >= zlo and max(zs) < zhi, "z coord escaped the box after wrapping"
    else:
        # atom-coordinate bbox + padding is NOT a real periodic cell -- it's only
        # sound if every atom is already wrapped into one periodic image, which an
        # unwrapped/whole-molecule PDB (e.g. chains kept continuous across PBC) will
        # violate. Prefer box_from with a real reference data file when you have one.
        pad = 10.0
        xlo, xhi = min(xs) - pad, max(xs) + pad
        ylo, yhi = min(ys) - pad, max(ys) + pad
        zlo, zhi = min(zs) - pad, max(zs) + pad

    overlaps = []
    if box is not None:
        atom_ids = list(range(1, atom_id + 1))
        atom_coords = list(zip(xs, ys, zs))
        overlaps = find_overlaps(atom_ids, atom_coords, bonded_pairs, box)

    with open(out_path, 'w') as f:
        f.write("LAMMPS data file -- combined via combine_lammps.py\n")
        if placeholder_species:
            f.write(f"# PLACEHOLDER ion species (need manual params from oplsaa.par): "
                    f"{placeholder_species}\n")
        if box_from is not None:
            f.write(f"# box bounds copied from {os.path.basename(box_from)}; atoms wrapped "
                    f"into it independently per-atom (standard convention) -- a molecule "
                    f"straddling the boundary will have atoms stored on opposite sides; "
                    f"bonded distances are still correct via minimum image\n")
        else:
            f.write("# box bounds are an auto bbox + padding, NOT a real periodic cell -- "
                    "likely wrong if the source PDB has unwrapped chains -- set manually\n")
        if overlaps:
            f.write(f"# WARNING: {len(overlaps)} non-bonded atom pair(s) closer than "
                    f"1.0 Angstrom after wrapping -- see stdout for the list\n")
        f.write("\n")

        f.write(f"{len(atoms_out):8d} atoms\n")
        f.write(f"{len(bonds_out):8d} bonds\n")
        f.write(f"{len(angles_out):8d} angles\n")
        f.write(f"{len(dihedrals_out):8d} dihedrals\n")
        f.write(f"{len(impropers_out):8d} impropers\n\n")

        f.write(f"{len(all_masses):8d} atom types\n")
        f.write(f"{len(all_bond_coeffs):8d} bond types\n")
        f.write(f"{len(all_angle_coeffs):8d} angle types\n")
        f.write(f"{len(all_dihedral_coeffs):8d} dihedral types\n")
        if all_improper_coeffs:
            f.write(f"{len(all_improper_coeffs):8d} improper types\n")
        f.write("\n")

        f.write(f"{xlo:12.6f} {xhi:12.6f} xlo xhi\n")
        f.write(f"{ylo:12.6f} {yhi:12.6f} ylo yhi\n")
        f.write(f"{zlo:12.6f} {zhi:12.6f} zlo zhi\n")

        f.write("\nMasses\n\n")
        for i, m in enumerate(all_masses, start=1):
            f.write(f"{i:8d} {m:10.3f}\n")

        f.write("\nPair Coeffs\n\n")
        for i, (eps, sig) in enumerate(all_pair, start=1):
            line = f"{i:8d} {eps:11.3f} {sig:11.7f}"
            if eps == PLACEHOLDER:
                line += "  # PLACEHOLDER"
            f.write(line + "\n")

        f.write("\nBond Coeffs\n\n")
        for i, (k0, r0) in enumerate(all_bond_coeffs, start=1):
            f.write(f"{i:8d} {k0:11.4f} {r0:11.4f}\n")

        f.write("\nAngle Coeffs\n\n")
        for i, (k0, a0) in enumerate(all_angle_coeffs, start=1):
            f.write(f"{i:8d} {k0:11.3f} {a0:11.3f}\n")

        f.write("\nDihedral Coeffs\n\n")
        for i, (v1, v2, v3, v4) in enumerate(all_dihedral_coeffs, start=1):
            f.write(f"{i:8d} {v1:11.3f} {v2:11.3f} {v3:11.3f} {v4:11.3f}\n")

        if all_improper_coeffs:
            f.write("\nImproper Coeffs\n\n")
            for i, (k, d0, n0) in enumerate(all_improper_coeffs, start=1):
                f.write(f"{i:8d} {k:11.3f} {d0:8d} {n0:8.1f}\n")

        f.write("\nAtoms\n\n")
        f.write("\n".join(atoms_out) + "\n")

        f.write("\nBonds\n\n")
        f.write("\n".join(bonds_out) + "\n")

        f.write("\nAngles\n\n")
        f.write("\n".join(angles_out) + "\n")

        f.write("\nDihedrals\n\n")
        f.write("\n".join(dihedrals_out) + "\n")

        if impropers_out:
            f.write("\nImpropers\n\n")
            f.write("\n".join(impropers_out) + "\n")

    print(f"Wrote {out_path}: {len(atoms_out)} atoms, {len(placeholder_species)} "
          f"placeholder ion species ({placeholder_species})")
    if placeholder_species:
        print(f"grep PLACEHOLDER {out_path}   # to find what needs manual fill-in")
    if overlaps:
        print(f"\nWARNING: {len(overlaps)} non-bonded atom pair(s) < 1.0 Angstrom apart "
              f"after wrapping (global atom ids, min-image distance):")
        for a, b, dist in overlaps[:20]:
            print(f"  {a:6d} -- {b:6d}   {dist:.3f} A")
        if len(overlaps) > 20:
            print(f"  ... and {len(overlaps) - 20} more")


def cli(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    box_from_pdb = '--box-from-pdb' in args
    if box_from_pdb:
        args.remove('--box-from-pdb')
    if len(args) not in (4, 5):
        print("Usage: python combine_lammps.py manifest.json original.pdb species_dir/ out.lmp "
              "[box_from.data] [--box-from-pdb]")
        print()
        print("  box_from.data   : copy box bounds from this LAMMPS data file")
        print("  --box-from-pdb  : read box from CRYST1 record in original.pdb "
              "(orthogonal cells only)")
        print()
        print("  Without either flag, box falls back to atom bbox + 10A padding,")
        print("  which is wrong for unwrapped/whole-molecule PDBs.")
        return 1
    run(*args[:4],
        box_from=(args[4] if len(args) == 5 else None),
        box_from_pdb=box_from_pdb)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
