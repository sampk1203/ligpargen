#!/usr/bin/env python3
"""
Step 1 of 2. Given a combined PSF+PDB (multiple molecules, e.g. polymer
chains + ions), find each unique molecule species via PSF connectivity,
and write one representative single-molecule PDB per species -- run
ligpargen on each of these, then feed the .lmp outputs + manifest.json
to combine_lammps.py to reassemble the full system with real coordinates.

Usage: python extract_unique_molecules.py input.psf input.pdb outdir/
Writes: outdir/species_1.pdb, outdir/species_2.pdb, ..., outdir/manifest.json
"""
import sys
import os
import json
from collections import Counter

from .species_graph import (
    connected_components, component_bonds, species_signature,
    group_into_species, estimate_net_charge, ION_ATOM_LIMIT,
)


def parse_psf(psf_path):
    with open(psf_path) as f:
        lines = f.readlines()

    natom = None
    natom_line_idx = None
    for i, l in enumerate(lines):
        if '!NATOM' in l:
            natom = int(l.split()[0])
            natom_line_idx = i
            break
    if natom is None:
        raise ValueError("No !NATOM section found in PSF")

    # PSF atom section: one line per atom, starting the line after !NATOM.
    # Column layout (whitespace-separated):
    #   0:serial  1:segname  2:resid  3:resname  4:atomname  5:atomtype  6:charge  7:mass  8:imove
    psf_charges = {}   # serial (1-based) -> float charge from PSF column 6
    for j in range(natom_line_idx + 1, natom_line_idx + 1 + natom):
        parts = lines[j].split()
        serial = int(parts[0])
        psf_charges[serial] = float(parts[6])

    nbond_idx = None
    for i, l in enumerate(lines):
        if '!NBOND' in l:
            nbond_idx = i
            nbonds = int(l.split()[0])
            break
    if nbond_idx is None:
        raise ValueError("No !NBOND section found in PSF")

    bond_nums = []
    i = nbond_idx + 1
    while len(bond_nums) < nbonds * 2:
        bond_nums.extend(int(x) for x in lines[i].split())
        i += 1
    bonds = [(bond_nums[j], bond_nums[j + 1]) for j in range(0, len(bond_nums), 2)]

    return natom, bonds, psf_charges


def parse_pdb_elements(pdb_path):
    """Element symbol is the last whitespace-separated token on each ATOM line."""
    elements = {}
    serial = 0
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                serial += 1
                elements[serial] = l.split()[-1]
    return elements


def read_pdb_atom_lines(pdb_path):
    with open(pdb_path) as f:
        lines = f.readlines()
    atom_lines = []
    serial = 0
    for l in lines:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            serial += 1
            atom_lines.append((serial, l))
    return atom_lines


def write_species_pdb(out_path, serials_in_order, atom_line_by_serial, local_bonds=()):
    # normalize resname/resSeq to a single uniform residue: this extracted
    # file is always exactly one molecule, so grouping everything under one
    # residue is safe here (unlike the combined multi-molecule file) and
    # avoids the "lone atom in its own residue" proximity-bonding bug that
    # broke ligpargen on EMC's per-monomer/per-terminal-atom residue naming
    #
    # local_bonds (from PSF, ground truth) are written as explicit CONECT
    # records for the same reason split_species.py's write_fragment_pdb
    # does it: without them RDKit falls back to proximityBonding, which on
    # F-dense species (e.g. PVDF-HFP) misreads close nonbonded F...F/F...C
    # contacts as real bonds, overfills F's valence, and MolFromPDBFile
    # returns None -- surfacing downstream as a cryptic MolToPDBFile
    # Boost.Python.ArgumentError. This was previously patched ad hoc per
    # element (e.g. carbon); writing real CONECT records here removes the
    # underlying cause for every element at once.
    with open(out_path, 'w') as f:
        for new_i, serial in enumerate(serials_in_order, start=1):
            l = atom_line_by_serial[serial]
            resfixed = l[:17] + 'MOLA' + l[21:22] + '   1' + l[26:]
            renumbered = resfixed[:6] + '%5d' % new_i + resfixed[11:]
            f.write(renumbered)

        neighbors = {}
        for a, b in local_bonds:
            neighbors.setdefault(a, []).append(b)
            neighbors.setdefault(b, []).append(a)
        for atom in sorted(neighbors):
            nbrs = sorted(neighbors[atom])
            # PDB CONECT caps at 4 partners per line; chunk if ever exceeded
            for k in range(0, len(nbrs), 4):
                chunk = nbrs[k:k + 4]
                f.write("CONECT" + f"{atom:5d}" + "".join(f"{n:5d}" for n in chunk) + "\n")

        f.write("END\n")


def run(psf_path, pdb_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    natom, bonds, psf_charges = parse_psf(psf_path)
    elements = parse_pdb_elements(pdb_path)

    atom_lines = read_pdb_atom_lines(pdb_path)
    if len(atom_lines) != natom:
        raise ValueError(f"PDB has {len(atom_lines)} atoms, PSF has {natom} -- files don't match")
    atom_line_by_serial = {serial: line for serial, line in atom_lines}

    comp_atoms, comp_to_species, species_list = group_into_species(natom, bonds, elements)

    print(f"Found {len(comp_atoms)} molecule(s), {len(species_list)} unique species:")
    manifest = {"species": {}, "instances": {}}
    for species_id, rep_comp_id, sig in species_list:
        serials = comp_atoms[rep_comp_id]
        n_atoms = sig[0]
        elem_counts = Counter(elements.get(s, '?') for s in serials)
        elem_str = ' '.join(f"{el}{ct}" for el, ct in elem_counts.items())
        n_bonds = len(component_bonds(bonds, set(serials)))
        print(f"  species {species_id}: {n_atoms} atoms ({elem_str}), {n_bonds} bonds "
              f"-- representative = molecule {rep_comp_id}")

        out_pdb = os.path.join(outdir, f"species_{species_id}.pdb")
        serial_to_local = {s: i for i, s in enumerate(serials, start=1)}
        local_bonds = [[serial_to_local[a], serial_to_local[b]]
                       for a, b in component_bonds(bonds, set(serials))]
        write_species_pdb(out_pdb, serials, atom_line_by_serial, local_bonds)

        # Derive expected net charge from PSF column 6 charges (set by EMC/psfgen).
        # IMPORTANT: PSF per-atom charges are OPLS force-field type charges, not
        # formal charges. For large polymer chains their partial charges can
        # accidentally sum close to an integer (e.g. 326-atom PVDF-HFP -> 4.08)
        # even though the molecule is formally neutral. See species_graph.py's
        # estimate_net_charge for the small-molecule-ion-vs-polymer-neutral rule.
        psf_charge_sum = sum(psf_charges.get(s, 0.0) for s in serials)
        expected_net_charge = estimate_net_charge(species_id, n_atoms, psf_charge_sum)

        manifest["species"][species_id] = {
            "pdb": os.path.basename(out_pdb),
            "n_atoms": n_atoms,
            "representative_molecule_id": rep_comp_id,
            "representative_atom_serials": serials,  # PSF/PDB serial order == local order fed to ligpargen
            "representative_bonds": local_bonds,  # local (1..n_atoms) index pairs, matches species_N.pdb order
            "expected_net_charge": expected_net_charge,  # integer; ligpargen -c flag and 1c charge target
        }

    for comp_id, serials in comp_atoms.items():
        manifest["instances"][comp_id] = {
            "species_id": comp_to_species[comp_id],
            "atom_serials": serials,  # order matches species template local order
        }

    with open(os.path.join(outdir, "manifest.json"), 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\nWrote {len(species_list)} representative PDB(s) and manifest.json to {outdir}/")
    print("Run ligpargen on each species_N.pdb, then run combine_lammps.py")


def cli(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 3:
        print("Usage: python extract_unique_molecules.py input.psf input.pdb outdir/")
        return 1
    run(args[0], args[1], args[2])
    return 0


if __name__ == "__main__":
    sys.exit(cli())
