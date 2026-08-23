#!/usr/bin/env python3
"""
mol2-native equivalent of extract_unique.py. Given a single mol2 file
(optionally several @<TRIPOS>MOLECULE blocks concatenated -- a packed
multi-molecule system, the usual output of packing tools), finds each
unique molecule species via mol2 connectivity/SYBYL types, writes one
representative single-molecule PDB per species plus a full-system PDB
(so downstream stages -- split, stitch, combine -- are unchanged and
format-agnostic from here on), and a manifest.json identical in schema
to extract_unique.py's.

Usage: python extract_unique_mol2.py input.mol2 outdir/
Writes: outdir/species_1.pdb, ..., outdir/system_from_mol2.pdb, outdir/manifest.json
"""
import sys
import os
import json
from collections import Counter

from .mol2_io import parse_mol2, write_pdb_from_mol2
from .species_graph import component_bonds, group_into_species, estimate_net_charge
from .extract_unique import read_pdb_atom_lines, write_species_pdb


def run(mol2_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    natom, bonds, elements, coords, charges = parse_mol2(mol2_path)

    system_pdb = os.path.join(outdir, "system_from_mol2.pdb")
    write_pdb_from_mol2(natom, elements, coords, system_pdb)

    atom_lines = read_pdb_atom_lines(system_pdb)
    atom_line_by_serial = {serial: line for serial, line in atom_lines}

    comp_atoms, comp_to_species, species_list = group_into_species(natom, bonds, elements)

    print(f"Found {len(comp_atoms)} molecule(s), {len(species_list)} unique species:")
    manifest = {"species": {}, "instances": {}, "source_format": "mol2",
                "system_pdb": os.path.basename(system_pdb)}
    for species_id, rep_comp_id, sig in species_list:
        serials = comp_atoms[rep_comp_id]
        n_atoms = sig[0]
        elem_counts = Counter(elements.get(s, '?') for s in serials)
        elem_str = ' '.join(f"{el}{ct}" for el, ct in elem_counts.items())
        n_bonds = len(component_bonds(bonds, set(serials)))
        print(f"  species {species_id}: {n_atoms} atoms ({elem_str}), {n_bonds} bonds "
              f"-- representative = molecule {rep_comp_id}")

        out_pdb = os.path.join(outdir, f"species_{species_id}.pdb")
        write_species_pdb(out_pdb, serials, atom_line_by_serial)
        serial_to_local = {s: i for i, s in enumerate(serials, start=1)}
        local_bonds = [[serial_to_local[a], serial_to_local[b]]
                       for a, b in component_bonds(bonds, set(serials))]

        charge_sum = sum(charges.get(s, 0.0) for s in serials)
        expected_net_charge = estimate_net_charge(species_id, n_atoms, charge_sum)

        manifest["species"][species_id] = {
            "pdb": os.path.basename(out_pdb),
            "n_atoms": n_atoms,
            "representative_molecule_id": rep_comp_id,
            "representative_atom_serials": serials,
            "representative_bonds": local_bonds,
            "expected_net_charge": expected_net_charge,
        }

    for comp_id, serials in comp_atoms.items():
        manifest["instances"][comp_id] = {
            "species_id": comp_to_species[comp_id],
            "atom_serials": serials,
        }

    with open(os.path.join(outdir, "manifest.json"), 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\nWrote {len(species_list)} representative PDB(s), {system_pdb}, "
          f"and manifest.json to {outdir}/")
    print("Run ligpargen on each species_N.pdb, then combine using "
          f"{system_pdb} as the coordinate source.")
    return system_pdb


def cli(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("Usage: molpargen extract-mol2 input.mol2 outdir/")
        return 1
    run(args[0], args[1])
    return 0


if __name__ == "__main__":
    sys.exit(cli())
