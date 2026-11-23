#!/usr/bin/env python3
"""
mol2-native equivalent of run_pipeline.py. Starts from a single mol2 file
(optionally a packed multi-molecule system as several concatenated
@<TRIPOS>MOLECULE blocks) instead of PSF+PDB. Runs extract_unique_mol2 for
stage 0, then delegates everything past that (ligpargen per species/
fragment, split, stitch, combine) to run_pipeline.run_from_manifest --
those stages are format-independent once species_dir/manifest.json and a
full-system PDB exist, regardless of whether the original input was
PSF+PDB or mol2.

Usage:
    python run_pipeline_mol2.py input.mol2 species_dir/ out.lmp [options]

Same options as run_pipeline.py's `run`, except there's no separate PDB
argument -- extract_unique_mol2 derives species_dir/system_from_mol2.pdb
from the mol2 itself and that's what combine uses for real coordinates.
"""
import sys
import os

from . import extract_unique_mol2
from . import run_pipeline
from .run_pipeline import (
    call_stage, BOSS_CAP, CGEN_DEFAULT, ORBITAL_TARGET_DEFAULT, WORKERS_DEFAULT,
)


def run(mol2_path, species_dir, out_lmp, max_atoms=BOSS_CAP,
        cgen=CGEN_DEFAULT, dry_run=False, force_split=False,
        orbital_target=ORBITAL_TARGET_DEFAULT, box_from=None, box_from_pdb=False,
        coarsen=True, workers=WORKERS_DEFAULT):

    # ── stage 0: extract unique species from mol2 (blocking) ─────────────────
    print(f"── extracting unique species from {mol2_path}")
    call_stage(dry_run, f"extract_unique_mol2.run({mol2_path!r}, {species_dir!r})",
               extract_unique_mol2.run, mol2_path, species_dir)

    manifest_path = os.path.join(species_dir, "manifest.json")
    system_pdb = os.path.join(species_dir, "system_from_mol2.pdb")
    if not dry_run and not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"extract stage ran but {manifest_path} not found")

    if dry_run:
        print(f"   [dry-run: would read {manifest_path} for species list]\n")
        print("DRY RUN: stopping after stage 0 (extraction) -- manifest.json "
              "doesn't exist yet to plan the remaining stages against.")
        return

    run_pipeline.run_from_manifest(
        manifest_path, species_dir, system_pdb, out_lmp,
        max_atoms=max_atoms, cgen=cgen, dry_run=dry_run,
        force_split=force_split, orbital_target=orbital_target,
        box_from=box_from, box_from_pdb=box_from_pdb,
        coarsen=coarsen, workers=workers)


def cli(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)

    dry_run     = '--dry-run'     in args; args = [a for a in args if a != '--dry-run']
    force_split = '--force-split' in args; args = [a for a in args if a != '--force-split']
    coarsen     = '--no-coarsen'  not in args; args = [a for a in args if a != '--no-coarsen']

    max_atoms = BOSS_CAP
    if '--max-atoms' in args:
        i = args.index('--max-atoms')
        max_atoms = int(args[i + 1])
        del args[i:i + 2]

    orbital_target = ORBITAL_TARGET_DEFAULT
    if '--max-orbitals' in args:
        i = args.index('--max-orbitals')
        orbital_target = int(args[i + 1])
        del args[i:i + 2]

    cgen = CGEN_DEFAULT
    if '--cgen' in args:
        i = args.index('--cgen')
        cgen = args[i + 1]
        del args[i:i + 2]

    workers = WORKERS_DEFAULT
    for flag in ('--workers', '-w'):
        if flag in args:
            i = args.index(flag)
            workers = int(args[i + 1])
            del args[i:i + 2]

    box_from_pdb = '--box-from-pdb' in args
    if box_from_pdb:
        args.remove('--box-from-pdb')
    box_from = None
    if '--box-from' in args:
        i = args.index('--box-from')
        box_from = args[i + 1]
        del args[i:i + 2]

    if len(args) != 3:
        print("Usage: molpargen run-mol2 input.mol2 species_dir/ out.lmp [options]")
        print()
        print("Same as 'molpargen run' but starting from a single mol2 file")
        print("instead of PSF+PDB. --box-from-pdb reads CRYST1 from the PDB")
        print("derived from the mol2 (species_dir/system_from_mol2.pdb) -- only")
        print("meaningful if the mol2 itself carries a real box in some form")
        print("your mol2 source wrote into that derived PDB; most mol2 files")
        print("don't carry box info at all, so --box-from a reference LAMMPS")
        print("data file is usually the more reliable option here.")
        print()
        print("Options:")
        print("  --max-atoms N          BOSS/ligpargen atom cap (default: 220)")
        print("  --max-orbitals N       BOSS orbital-count cap (default: 480)")
        print("  --cgen CM1A|CM1A-LBCC  charge model (default: CM1A)")
        print("  --dry-run              print stages without running them")
        print("  --force-split          pass --force to the splitter")
        print("  --box-from FILE        copy box bounds from this LAMMPS data file")
        print("  --box-from-pdb         read box from derived system PDB's CRYST1 record")
        print("  --no-coarsen           skip the post-combine coarsening step (on by default);")
        print("                         collapses per-atom types to (mol × element), writing")
        print("                         <out>_coarse.data + <out>_topology.json")
        print("  --workers N, -w N      run up to N ligpargen calls in parallel (default: 1,")
        print("                         i.e. serial -- same behavior as before this flag existed)")
        return 1

    mol2_path, species_dir, out_lmp = args
    run(mol2_path, species_dir, out_lmp,
        max_atoms=max_atoms, cgen=cgen, dry_run=dry_run,
        force_split=force_split, orbital_target=orbital_target,
        box_from=box_from, box_from_pdb=box_from_pdb,
        coarsen=coarsen, workers=workers)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
