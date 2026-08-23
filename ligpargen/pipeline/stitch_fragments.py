#!/usr/bin/env python3
"""
Step 1c (runs after 1b_split + ligpargen on each fragment).
Merges per-fragment ligpargen .lmp outputs back into a single
species_N.lammps.lmp in the same format parse_species_lmp() in
2_combine_lammps.py reads -- so 2_combine sees it identically to
a whole-molecule ligpargen run.

Usage:
    python 1c_stitch_fragments.py manifest.json fragment_manifest_N.json \\
           fragment_dir/ species_dir/ [--charge-tol F]

    manifest.json              -- from 1_extract_unique_molecules.py
    fragment_manifest_N.json   -- from 1b_split_large_species.py
    fragment_dir/              -- directory containing fragment_N_K.lammps.lmp
                                  (ligpargen run: ligpargen -i fragment_N_K.pdb
                                                            -n fragment_N_K
                                                            -c <charge> -o 0)
    species_dir/               -- where to write species_N.lammps.lmp

    --charge-tol F             -- warn (don't fail) if residual charge after
                                  redistribution exceeds F (default 0.5)

Output:
    species_dir/species_N.lammps.lmp

Algorithm (see spec in conversation):
  1. Parse each fragment_N_K.lammps.lmp via parse_species_lmp().
  2. For each fragment, from fragment_manifest: local_to_species, core_species_atoms,
     cap_local_ids.
  3. Drop every bond/angle/dihedral/improper that touches any cap atom.
  4. Keep only atoms/terms whose species-level atom(s) are in core_species_atoms
     (ownership partition from 1b is already unambiguous -- no heuristic detection).
  5. Reindex merged numbering to 1..n_atoms (original species_N.pdb order).
  6. Type-remap: each fragment's ligpargen assigns fresh type ids -- accumulate offsets
     fragment by fragment, exactly as 2_combine_lammps.py does across species.
  7. Sum real-atom charges, compare to expected_net_charge from manifest,
     distribute residual equally across all atoms.
  8. Sanity asserts before write.
  9. Write in parse_species_lmp's expected section order.
"""
import sys
import os
import json

from .combine_lammps import parse_species_lmp


# ── helpers ──────────────────────────────────────────────────────────────────

def load_fragment_lmp(fragment_dir, species_id, fragment_id):
    name = f"fragment_{species_id}_{fragment_id}.lammps.lmp"
    path = os.path.join(fragment_dir, name)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Missing {path} -- run ligpargen on fragment_{species_id}_{fragment_id}.pdb first.\n"
            f"  ligpargen -i fragment_{species_id}_{fragment_id}.pdb "
            f"-n fragment_{species_id}_{fragment_id} -c <charge> -o 0")
    return parse_species_lmp(path)


def read_species_pdb_coords(pdb_path):
    """Return {1-based local index: (x,y,z)} from a species_N.pdb.
    Order matches the atom order used throughout this pipeline."""
    coords = {}
    idx = 0
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                idx += 1
                coords[idx] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return coords


def term_touches_cap(indices, cap_local_ids_set):
    """Return True if any atom index in a bond/angle/dihedral/improper tuple
    is a cap atom (caps are never carried into the merged file)."""
    return bool(set(indices) & cap_local_ids_set)


# ── main ─────────────────────────────────────────────────────────────────────

def run(manifest_path, frag_manifest_path, fragment_dir, species_dir,
        charge_tol=0.5):

    # ── load manifests ───────────────────────────────────────────────────────
    with open(manifest_path) as f:
        manifest = json.load(f)
    with open(frag_manifest_path) as f:
        frag_manifest = json.load(f)

    species_id  = frag_manifest['species_id']
    n_atoms     = frag_manifest['n_atoms']
    sid_str     = str(species_id)

    sp = manifest['species'].get(sid_str)
    if sp is None:
        raise KeyError(f"species {species_id} not found in {manifest_path}")

    expected_net_charge = sp.get('expected_net_charge')
    if expected_net_charge is None:
        print(f"WARNING: 'expected_net_charge' missing from manifest for species "
              f"{species_id} -- charge redistribution check will be skipped.\n"
              f"  Re-run molpargen extract with the updated version "
              f"to add PSF-derived net charges to the manifest.")

    frag_ids = sorted(int(k) for k in frag_manifest['fragments'])
    print(f"species {species_id}: stitching {len(frag_ids)} fragment(s), "
          f"n_atoms = {n_atoms}")

    # ── per-fragment parse + filter ──────────────────────────────────────────
    # After filtering we accumulate into species-indexed dicts/lists.
    # Type ids are remapped fragment-by-fragment with running offsets
    # (same pattern as 2_combine_lammps.py's off_a/off_b/... loop).

    # Accumulators (species-indexed):
    #   masses[species_local_idx] = float
    #   pair[species_local_idx]   = (epsilon, sigma)
    #   charges[species_local_idx]= float
    #   atom_type[species_local_idx] = merged_type_id (1-based)
    #
    # These are filled as each fragment's core atoms are claimed.
    # "merged type id" is assigned per fragment block via offset accumulation.

    masses_by_sidx  = {}   # species_local_idx (1..n_atoms) -> mass
    pair_by_sidx    = {}   # species_local_idx -> (eps, sig)
    charges_by_sidx = {}   # species_local_idx -> charge

    # Topology terms + their coeffs are collected only for kept (owned) terms.
    # Type IDs are sequential 1..N with no gaps so parse_species_lmp's
    # (n_types == n_atoms for atoms; n_bond_types == n_bonds; etc.) holds.
    bonds_species     = []   # (sidx_a, sidx_b)
    angles_species    = []   # (sidx_a, sidx_b, sidx_c)
    dihedrals_species = []
    impropers_species = []
    all_bond_coeffs   = []
    all_angle_coeffs  = []
    all_dihedral_coeffs  = []
    all_improper_coeffs  = []

    # Boundary terms (bond/angle/dihedral straddling a cut point) exist as
    # real, uncapped atoms in TWO adjacent fragments (trust_radius overlap)
    # but are core-owned by only one atom's-worth in each -- so ownership
    # can't gate keeping them, only cap-touching can. Dedup by canonical
    # species-index form so an overlap-visible term is written once.
    seen_bonds, seen_angles, seen_dih, seen_imp = set(), set(), set(), set()

    def canon_angle(a, b, c):
        return (a, b, c) if a < c else (c, b, a)

    def canon_torsion(a, b, c, d):
        t = (a, b, c, d)
        return t if t <= t[::-1] else t[::-1]

    for fid in frag_ids:
        finfo = frag_manifest['fragments'][str(fid)]
        n_real          = finfo['n_real']
        cap_local_ids   = set(finfo['cap_local_ids'])
        # JSON stores keys as strings; convert to int.
        local_to_species = {int(k): int(v) for k, v in finfo['local_to_species'].items()}
        core_species_set = set(finfo['core_species_atoms'])

        lmp = load_fragment_lmp(fragment_dir, species_id, fid)
        n_frag_types = len(lmp['masses'])  # = n_atoms_in_fragment (ligpargen: 1 type/atom)

        print(f"  fragment {fid}: {finfo['n_atoms']} atoms "
              f"({n_real} real + {finfo['n_caps']} caps), "
              f"{n_frag_types} atom types, owns {len(core_species_set)} species atoms")

        # -- Atom params for this fragment's core atoms only -----------------
        # ligpargen: atom i has type i. Caps and non-owned overlap atoms are
        # skipped so we never create orphan type entries.
        for local_id in range(1, n_real + 1):
            if local_id in cap_local_ids:
                continue
            sidx = local_to_species[local_id]
            if sidx not in core_species_set:
                continue

            frag_type = local_id
            masses_by_sidx[sidx]  = lmp['masses'][frag_type - 1]
            pair_by_sidx[sidx]    = lmp['pair'][frag_type - 1]
            charges_by_sidx[sidx] = lmp['charges'][local_id - 1]

        # -- Topology terms (no caps; ownership no longer gates keeping --
        # only whether a term touches a cap atom does. A boundary term is
        # real/uncapped in whichever fragment(s) cover it via overlap; the
        # seen_* dedup below keeps it exactly once regardless of how many
        # fragments see it. Type IDs assigned later as sequential 1..N so
        # n_bond_types == n_bonds (ligpargen convention). --------------------
        for bidx, (a, b) in enumerate(lmp['bonds']):
            if term_touches_cap((a, b), cap_local_ids):
                continue
            sa, sb = local_to_species[a], local_to_species[b]
            key = frozenset((sa, sb))
            if key in seen_bonds:
                continue
            seen_bonds.add(key)
            bonds_species.append((sa, sb))
            all_bond_coeffs.append(lmp['bond_coeffs'][bidx])

        for aidx, (a, b, c) in enumerate(lmp['angles']):
            if term_touches_cap((a, b, c), cap_local_ids):
                continue
            sa, sb, sc = local_to_species[a], local_to_species[b], local_to_species[c]
            key = canon_angle(sa, sb, sc)
            if key in seen_angles:
                continue
            seen_angles.add(key)
            angles_species.append((sa, sb, sc))
            all_angle_coeffs.append(lmp['angle_coeffs'][aidx])

        for didx, (a, b, c, e) in enumerate(lmp['dihedrals']):
            if term_touches_cap((a, b, c, e), cap_local_ids):
                continue
            sa, sb, sc, se = (local_to_species[a], local_to_species[b],
                               local_to_species[c], local_to_species[e])
            key = canon_torsion(sa, sb, sc, se)
            if key in seen_dih:
                continue
            seen_dih.add(key)
            dihedrals_species.append((sa, sb, sc, se))
            all_dihedral_coeffs.append(lmp['dihedral_coeffs'][didx])

        for iidx, (a, b, c, e) in enumerate(lmp['impropers']):
            if term_touches_cap((a, b, c, e), cap_local_ids):
                continue
            sa, sb, sc, se = (local_to_species[a], local_to_species[b],
                               local_to_species[c], local_to_species[e])
            key = canon_torsion(sa, sb, sc, se)
            if key in seen_imp:
                continue
            seen_imp.add(key)
            impropers_species.append((sa, sb, sc, se))
            all_improper_coeffs.append(lmp['improper_coeffs'][iidx])

    # ── sanity checks before writing ────────────────────────────────────────

    # 1. Every species atom was claimed by exactly one fragment.
    missing = [i for i in range(1, n_atoms + 1) if i not in charges_by_sidx]
    duplicate_check = len(charges_by_sidx)
    if missing:
        raise AssertionError(
            f"species {species_id}: {len(missing)} atom(s) never claimed by any fragment: "
            f"{missing[:10]}{'...' if len(missing) > 10 else ''}")
    if duplicate_check != n_atoms:
        raise AssertionError(
            f"species {species_id}: expected {n_atoms} owned atoms, got {duplicate_check} "
            f"-- ownership overlap or count mismatch")

    # 2. No bond/angle/dihedral/improper references an out-of-range index.
    all_term_atoms = set()
    for a, b in bonds_species:                all_term_atoms |= {a, b}
    for a, b, c in angles_species:            all_term_atoms |= {a, b, c}
    for a, b, c, e in dihedrals_species:      all_term_atoms |= {a, b, c, e}
    for a, b, c, e in impropers_species:      all_term_atoms |= {a, b, c, e}
    bad = [i for i in all_term_atoms if i < 1 or i > n_atoms]
    if bad:
        raise AssertionError(
            f"species {species_id}: topology terms reference out-of-range species indices: "
            f"{sorted(bad)[:10]}")

    # ── build sequential type tables (1 type per atom, ligpargen style) ─────
    # parse_species_lmp requires len(masses) == len(pair) == len(charges) == n_atoms.
    # Atom type for species atom i is i (1-based).
    all_masses_list = [masses_by_sidx[i] for i in range(1, n_atoms + 1)]
    all_pair_list   = [pair_by_sidx[i]   for i in range(1, n_atoms + 1)]
    atype_by_sidx   = {i: i for i in range(1, n_atoms + 1)}

    # ── charge redistribution ───────────────────────────────────────────────
    charges_list = [charges_by_sidx[i] for i in range(1, n_atoms + 1)]
    raw_sum = sum(charges_list)

    if expected_net_charge is not None:
        residual = expected_net_charge - raw_sum
        residual_per_atom = residual / n_atoms
        charges_list = [q + residual_per_atom for q in charges_list]
        adj_sum = sum(charges_list)
        print(f"  charge: raw sum = {raw_sum:.6f}, target = {expected_net_charge}, "
              f"residual = {residual:.6f} ({residual:.4f}e distributed equally, "
              f"{residual_per_atom:.6f}e/atom)")
        if abs(residual) > charge_tol:
            print(f"  WARNING: residual |{residual:.4f}|e > tolerance {charge_tol}e -- "
                  f"fragment charge assignments may be inconsistent. "
                  f"Check -c flags passed to ligpargen for each fragment.")
        print(f"  adjusted sum = {adj_sum:.8f} (should be {expected_net_charge})")
    else:
        print(f"  charge: raw sum = {raw_sum:.6f} (no expected_net_charge in manifest, "
              f"skipping redistribution)")

    # ── write species_N.lammps.lmp ──────────────────────────────────────────
    # Single-structure LAMMPS data file. Atom types are 1..n_atoms (one per
    # atom). Bond/angle/dihedral/improper types are 1..N for kept terms only
    # so parse_species_lmp's n_types == n_terms invariant holds.
    # Coords from species_N.pdb so the file is a real molecule, not zeros.
    os.makedirs(species_dir, exist_ok=True)
    out_path = os.path.join(species_dir, f"species_{species_id}.lammps.lmp")

    pdb_name = sp.get('pdb', f"species_{species_id}.pdb")
    pdb_path = os.path.join(species_dir, pdb_name)
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(
            f"Cannot write real coordinates: species PDB missing: {pdb_path}\n"
            f"  Re-run 1_extract so species_{species_id}.pdb sits in species_dir.")
    pdb_coords = read_species_pdb_coords(pdb_path)
    if len(pdb_coords) != n_atoms:
        raise ValueError(
            f"{pdb_path} has {len(pdb_coords)} atoms, fragment_manifest says n_atoms={n_atoms}")

    xs = [pdb_coords[i][0] for i in range(1, n_atoms + 1)]
    ys = [pdb_coords[i][1] for i in range(1, n_atoms + 1)]
    zs = [pdb_coords[i][2] for i in range(1, n_atoms + 1)]
    pad = 5.0
    xlo, xhi = min(xs) - pad, max(xs) + pad
    ylo, yhi = min(ys) - pad, max(ys) + pad
    zlo, zhi = min(zs) - pad, max(zs) + pad

    n_atom_types     = len(all_masses_list)       # == n_atoms
    n_bond_types     = len(all_bond_coeffs)       # == len(bonds_species)
    n_angle_types    = len(all_angle_coeffs)
    n_dihedral_types = len(all_dihedral_coeffs)
    n_improper_types = len(all_improper_coeffs)

    with open(out_path, 'w') as f:
        f.write(f"LAMMPS data file -- stitched by molpargen stitch "
                f"(species {species_id}, {len(frag_ids)} fragments)\n")
        f.write("# Single molecular structure. Coords from species PDB; "
                "combine replaces them from the full system.\n\n")

        f.write(f"{n_atoms:8d} atoms\n")
        f.write(f"{len(bonds_species):8d} bonds\n")
        f.write(f"{len(angles_species):8d} angles\n")
        f.write(f"{len(dihedrals_species):8d} dihedrals\n")
        f.write(f"{len(impropers_species):8d} impropers\n\n")

        f.write(f"{n_atom_types:8d} atom types\n")
        f.write(f"{n_bond_types:8d} bond types\n")
        f.write(f"{n_angle_types:8d} angle types\n")
        f.write(f"{n_dihedral_types:8d} dihedral types\n")
        if n_improper_types:
            f.write(f"{n_improper_types:8d} improper types\n")
        f.write("\n")

        f.write(f"{xlo:12.6f} {xhi:12.6f} xlo xhi\n")
        f.write(f"{ylo:12.6f} {yhi:12.6f} ylo yhi\n")
        f.write(f"{zlo:12.6f} {zhi:12.6f} zlo zhi\n")

        f.write("\nMasses\n\n")
        for i, m in enumerate(all_masses_list, start=1):
            f.write(f"{i:8d} {m:10.3f}\n")

        f.write("\nPair Coeffs\n\n")
        for i, (eps, sig) in enumerate(all_pair_list, start=1):
            f.write(f"{i:8d} {eps:11.3f} {sig:11.7f}\n")

        f.write("\nBond Coeffs\n\n")
        for i, (k0, r0) in enumerate(all_bond_coeffs, start=1):
            f.write(f"{i:8d} {k0:11.4f} {r0:11.4f}\n")

        f.write("\nAngle Coeffs\n\n")
        for i, (k0, a0) in enumerate(all_angle_coeffs, start=1):
            f.write(f"{i:8d} {k0:11.3f} {a0:11.3f}\n")

        f.write("\nDihedral Coeffs\n\n")
        for i, (v1, v2, v3, v4) in enumerate(all_dihedral_coeffs, start=1):
            f.write(f"{i:8d} {v1:11.3f} {v2:11.3f} {v3:11.3f} {v4:11.3f}\n")

        if n_improper_types:
            f.write("\nImproper Coeffs\n\n")
            for i, (k, d0, n0) in enumerate(all_improper_coeffs, start=1):
                f.write(f"{i:8d} {k:11.3f} {int(d0):8d} {n0:8.1f}\n")

        f.write("\nAtoms\n\n")
        for sidx in range(1, n_atoms + 1):
            atype  = atype_by_sidx[sidx]   # == sidx
            charge = charges_list[sidx - 1]
            x, y, z = pdb_coords[sidx]
            f.write(f"{sidx:8d} {1:6d} {atype:6d} {charge:10.6f} "
                    f"{x:10.6f} {y:10.6f} {z:10.6f}\n")

        f.write("\nBonds\n\n")
        for bidx, (a, b) in enumerate(bonds_species, start=1):
            f.write(f"{bidx:6d} {bidx:6d} {a:6d} {b:6d}\n")  # type == index

        f.write("\nAngles\n\n")
        for aidx, (a, b, c) in enumerate(angles_species, start=1):
            f.write(f"{aidx:6d} {aidx:6d} {a:6d} {b:6d} {c:6d}\n")

        f.write("\nDihedrals\n\n")
        for didx, (a, b, c, e) in enumerate(dihedrals_species, start=1):
            f.write(f"{didx:6d} {didx:6d} {a:6d} {b:6d} {c:6d} {e:6d}\n")

        if impropers_species:
            f.write("\nImpropers\n\n")
            for iidx, (a, b, c, e) in enumerate(impropers_species, start=1):
                f.write(f"{iidx:6d} {iidx:6d} {a:6d} {b:6d} {c:6d} {e:6d}\n")

        f.write("\n")

    print(f"\nWrote {out_path}")
    print(f"  {n_atoms} atoms, {len(bonds_species)} bonds, {len(angles_species)} angles, "
          f"{len(dihedrals_species)} dihedrals, {len(impropers_species)} impropers")
    print(f"  {n_atom_types} atom types, {n_bond_types} bond types, "
          f"{n_angle_types} angle types, {n_dihedral_types} dihedral types"
          + (f", {n_improper_types} improper types" if n_improper_types else ""))
    print(f"  coordinates from {pdb_name}; atom types = n_atoms (no orphan types)")
    print(f"Ready for molpargen combine (species_dir = {species_dir})")


def cli(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    charge_tol = 0.5
    if '--charge-tol' in args:
        i = args.index('--charge-tol')
        charge_tol = float(args[i + 1])
        del args[i:i + 2]
    if len(args) != 4:
        print("Usage: molpargen stitch manifest.json "
              "fragment_manifest_N.json fragment_dir/ species_dir/ [--charge-tol F]")
        return 1
    run(args[0], args[1], args[2], args[3], charge_tol=charge_tol)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
