#!/usr/bin/env python3
"""
Step 1d -- full pipeline orchestrator, entry point and exit point.
Starts from the raw PSF+PDB (runs 1_extract_unique_molecules.py to produce
species_dir/manifest.json + species_N.pdb), decides for each species
whether it fits under the BOSS/ligpargen atom cap (default 220), runs
1b_split_large_species.py for those that don't, runs ligpargen on every
fragment/small-species PDB, runs 1c_stitch_fragments.py for every split
species, then runs 2_combine_lammps.py. Each stage blocks on the previous
one (subprocess.run, not backgrounded) -- this script returns only after
out.lmp exists, or raises on the first failing stage.

Species with n_atoms == 1 (bare ions: Li+, Na+, ...) are skipped entirely --
2_combine_lammps.py synthesizes those with placeholder LJ params.

Usage:
    python 1d_run_ligpargen.py input.psf input.pdb species_dir/ out.lmp [options]

Options:
    --max-atoms N       ligpargen/BOSS atom cap (default: 220)
    --cgen {CM1A,CM1A-LBCC}
                        charge model passed to ligpargen (default: CM1A)
    --dry-run           print commands without executing them
    --force-split       pass --force to 1b for species with flagged backbone atoms
    --scripts-dir DIR   directory containing 1b_split_large_species.py and
                        1c_stitch_fragments.py (default: same dir as this script)

Assumes:
    - ligpargen is on PATH
    - 1b_split_large_species.py and 1c_stitch_fragments.py are in --scripts-dir
    - All ligpargen runs happen from within species_dir/ so outputs land there
"""
import sys
import os
import json
import subprocess

from . import extract_unique
from . import split_species
from . import stitch_fragments
from . import combine_lammps
from ..tools.coarsen import coarsen as _coarsen

BOSS_CAP      = 220   # default max atoms per ligpargen run
CGEN_DEFAULT  = 'CM1A'
ORBITAL_TARGET_DEFAULT = 480   # passed through to 1b as --max-orbitals; see that
                                # script for why atom count alone isn't a safe budget

# Mirrors 1b_split_large_species.py's ORBITAL_WEIGHTS -- kept as a small
# duplicate here (rather than importing 1b as a module) so this script can
# sanity-check species that go straight to ligpargen WITHOUT ever going
# through 1b's splitter, which is the only place that currently checks
# BOSS's real orbital cap. A small, F/S-light species can still be under the
# atom cap yet over the orbital cap if it happens to be dense in F/S/P.
ORBITAL_WEIGHTS = {'H': 1, 'HE': 1,
                    'LI': 4, 'BE': 4, 'B': 4, 'C': 4, 'N': 4, 'O': 4, 'F': 4, 'NE': 4,
                    'NA': 4, 'MG': 4, 'AL': 4, 'SI': 4, 'CL': 4, 'BR': 4, 'I': 4}
DEFAULT_ORBITAL_WEIGHT = 9
BOSS_ORBITAL_CAP = 550

_HERE = os.path.dirname(os.path.abspath(__file__))


def estimate_pdb_orbitals(pdb_path):
    """Rough element-weighted orbital estimate for a species PDB, read directly
    (element is the last whitespace token on each ATOM/HETATM line, same
    convention as the rest of this pipeline)."""
    total = 0
    with open(pdb_path) as f:
        for l in f:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                el = l.split()[-1].upper()
                total += ORBITAL_WEIGHTS.get(el, DEFAULT_ORBITAL_WEIGHT)
    return total


def run_cmd(cmd, dry_run, cwd=None):
    """Print and optionally execute an external shell command (ligpargen/BOSS)."""
    print('  $', ' '.join(cmd))
    if not dry_run:
        result = subprocess.run(cmd, cwd=cwd)
        if result.returncode != 0:
            raise RuntimeError(
                f"Command failed (exit {result.returncode}): {' '.join(cmd)}")


def call_stage(dry_run, desc, fn, *args, **kwargs):
    """Print a description and, unless dry_run, call an in-process pipeline
    stage directly (extract/split/stitch/combine -- no more subprocess/sys.path
    tricks to invoke sibling scripts)."""
    print('  ›', desc)
    if not dry_run:
        fn(*args, **kwargs)


def ligpargen_cmd(pdb_name, mol_name, charge, cgen):
    """Build ligpargen invocation. Run from species_dir so output lands there."""
    return [
        'ligpargen',
        '-i', pdb_name,
        '-n', mol_name,
        '-c', str(charge),
        '-o', '0',
        '-cgen', cgen,
    ]


def run(psf_path, pdb_path, species_dir, out_lmp, max_atoms=BOSS_CAP,
        cgen=CGEN_DEFAULT, dry_run=False, force_split=False,
        orbital_target=ORBITAL_TARGET_DEFAULT, box_from=None, box_from_pdb=False,
        coarsen=True):

    # ── stage 0: extract unique species from PSF+PDB (blocking) ──────────────
    print(f"── extracting unique species from {psf_path} / {pdb_path}")
    call_stage(dry_run, f"extract_unique.run({psf_path!r}, {pdb_path!r}, {species_dir!r})",
               extract_unique.run, psf_path, pdb_path, species_dir)

    manifest_path = os.path.join(species_dir, "manifest.json")
    if not dry_run and not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"extract stage ran but {manifest_path} not found")
    original_pdb = pdb_path  # same PDB feeds real coordinates to combine

    if dry_run:
        print(f"   [dry-run: would read {manifest_path} for species list]\n")
        # No manifest exists yet in dry-run -- can't enumerate species below,
        # so dry-run only demonstrates stage 0 and stops here.
        print("DRY RUN: stopping after stage 0 (extraction) -- manifest.json "
              "doesn't exist yet to plan the remaining stages against.")
        return

    run_from_manifest(manifest_path, species_dir, original_pdb, out_lmp,
                       max_atoms=max_atoms, cgen=cgen, dry_run=dry_run,
                       force_split=force_split, orbital_target=orbital_target,
                       box_from=box_from, box_from_pdb=box_from_pdb,
                       coarsen=coarsen)


def run_from_manifest(manifest_path, species_dir, original_pdb, out_lmp,
                       max_atoms=BOSS_CAP, cgen=CGEN_DEFAULT, dry_run=False,
                       force_split=False, orbital_target=ORBITAL_TARGET_DEFAULT,
                       box_from=None, box_from_pdb=False, coarsen=True):
    """Everything past stage 0 (extraction): ligpargen per species/fragment,
    split, stitch, combine. Format-independent -- doesn't care whether
    manifest.json/species_N.pdb came from PSF+PDB (extract_unique) or mol2
    (extract_unique_mol2); only needs a manifest, a species_dir holding the
    species PDBs, and a full-system PDB (original_pdb) for combine's real
    coordinates."""

    with open(manifest_path) as f:
        manifest = json.load(f)

    species_to_stitch = []   # (species_id, fragment_dir) -- need stitch after
    direct_species = []   # species_id -- went straight to ligpargen
    ion_species    = []   # species_id -- skipped (1-atom, handled by combine)

    print(f"manifest: {manifest_path}")
    print(f"species_dir: {species_dir}")
    print(f"ligpargen cap: {max_atoms} atoms, charge model: {cgen}")
    print(f"dry_run: {dry_run}\n")

    for sid_str in sorted(manifest['species'], key=lambda x: int(x)):
        sid = int(sid_str)
        sp  = manifest['species'][sid_str]
        n   = sp['n_atoms']
        charge = sp.get('expected_net_charge', 0)

        print(f"── species {sid}: {n} atoms, expected_net_charge={charge}")

        # ── 1-atom ion: skip ─────────────────────────────────────────────────
        if n == 1:
            print(f"   → bare ion, skipping (combine stage synthesizes placeholder)")
            ion_species.append(sid)
            continue

        pdb_name = sp['pdb']   # e.g. "species_1.pdb"
        mol_name = os.path.splitext(pdb_name)[0]  # "species_1"

        # ── small enough: run ligpargen directly ─────────────────────────────
        if n <= max_atoms:
            est_orb = estimate_pdb_orbitals(os.path.join(species_dir, pdb_name))
            if est_orb > BOSS_ORBITAL_CAP:
                raise ValueError(
                    f"species {sid}: {n} atoms is under the atom cap ({max_atoms}) but its "
                    f"estimated orbital count ({est_orb}) exceeds BOSS's hard cap of "
                    f"{BOSS_ORBITAL_CAP} -- this species is unusually dense in heavy elements "
                    f"(F/S/P/...). Lower --max-atoms so it goes through the 1b splitter instead "
                    f"of running directly.")
            print(f"   → {n} <= {max_atoms} ({est_orb} est. orbitals): running ligpargen directly")
            cmd = ligpargen_cmd(pdb_name, mol_name, charge, cgen)
            run_cmd(cmd, dry_run, cwd=species_dir)
            lmp_out = os.path.join(species_dir, f"{mol_name}.lammps.lmp")
            if not dry_run and not os.path.exists(lmp_out):
                raise FileNotFoundError(
                    f"ligpargen ran but output not found: {lmp_out}\n"
                    f"Check ligpargen output above for errors.")
            direct_species.append(sid)
            print(f"   → wrote {mol_name}.lammps.lmp")
            continue

        # ── too large: split with 1b then ligpargen per fragment ─────────────
        print(f"   → {n} > {max_atoms}: splitting (split_species.run)")
        frag_dir = os.path.join(species_dir, f"fragments_{sid}")

        call_stage(dry_run, f"split_species.run(species {sid} -> {frag_dir})",
                   split_species.run, manifest_path, species_dir, sid, frag_dir,
                   max_atoms, split_species.TRUST_RADIUS, force_split,
                   orbital_target=orbital_target)

        # Load fragment manifest to know which fragments were created
        frag_manifest_path = os.path.join(frag_dir, f"fragment_manifest_{sid}.json")
        if dry_run:
            print(f"   [dry-run: would read {frag_manifest_path} for fragment list]")
            print(f"   [dry-run: would run ligpargen on each fragment_{sid}_K.pdb]")
            species_to_stitch.append((sid, frag_dir))
            continue

        if not os.path.exists(frag_manifest_path):
            # 1b printed "no split needed" -- species fit in one fragment after all
            print(f"   → 1b reported no split needed; running ligpargen directly")
            cmd = ligpargen_cmd(pdb_name, mol_name, charge, cgen)
            run_cmd(cmd, dry_run, cwd=species_dir)
            direct_species.append(sid)
            continue

        with open(frag_manifest_path) as f:
            frag_manifest = json.load(f)

        frag_ids = sorted(int(k) for k in frag_manifest['fragments'])
        print(f"   → {len(frag_ids)} fragments created")

        for fid in frag_ids:
            finfo = frag_manifest['fragments'][str(fid)]
            frag_pdb  = finfo['pdb']                          # fragment_N_K.pdb
            frag_name = os.path.splitext(frag_pdb)[0]         # fragment_N_K
            # Fragment charge: use 0 for polymer (ligpargen -c 0 is always right
            # for a neutral fragment regardless of which slice; cap H's are neutral
            # and the redistribute step in 1c corrects any sum drift).
            frag_charge = 0
            print(f"   fragment {fid}: {finfo['n_atoms']} atoms "
                  f"({finfo['n_real']} real + {finfo['n_caps']} caps)")
            cmd = ligpargen_cmd(frag_pdb, frag_name, frag_charge, cgen)
            run_cmd(cmd, dry_run, cwd=frag_dir)
            lmp_out = os.path.join(frag_dir, f"{frag_name}.lammps.lmp")
            if not os.path.exists(lmp_out):
                raise FileNotFoundError(
                    f"ligpargen ran but output not found: {lmp_out}\n"
                    f"Check ligpargen output above for errors.")
            print(f"   → wrote {frag_name}.lammps.lmp")

        species_to_stitch.append((sid, frag_dir))

    # ── stage: stitch every split species (blocking, one at a time) ──────────
    print("\n" + "═" * 60)
    print(f"DONE with ligpargen runs. {len(species_to_stitch)} species need stitching.\n")

    for sid, frag_dir in species_to_stitch:
        frag_manifest_path = os.path.join(frag_dir, f"fragment_manifest_{sid}.json")
        print(f"── stitching species {sid}")
        call_stage(dry_run, f"stitch_fragments.run(species {sid})",
                   stitch_fragments.run, manifest_path, frag_manifest_path,
                   frag_dir + '/', species_dir + '/')
        stitched = os.path.join(species_dir, f"species_{sid}.lammps.lmp")
        if not dry_run and not os.path.exists(stitched):
            raise FileNotFoundError(f"1c ran but output not found: {stitched}")

    if ion_species:
        print(f"\nIon species skipped by ligpargen (placeholder params): {ion_species}")
        print(f"  combine stage synthesizes these; grep PLACEHOLDER {out_lmp} "
              f"after combine and fill from oplsaa.par")

    # ── stage: combine (blocking) -- this is the last step; 1d only returns
    # after out_lmp actually exists ────────────────────────────────────────
    print(f"\n── combining into {out_lmp}")
    call_stage(dry_run, f"combine_lammps.run(-> {out_lmp})",
               combine_lammps.run, manifest_path, original_pdb, species_dir, out_lmp,
               box_from=box_from, box_from_pdb=box_from_pdb)
    if not dry_run and not os.path.exists(out_lmp):
        raise FileNotFoundError(f"combine stage ran but output not found: {out_lmp}")

    if not dry_run and coarsen:
        print(f"\n── coarsening {out_lmp}")
        _coarsen(out_lmp)

    print("\n" + "═" * 60)
    print(f"PIPELINE COMPLETE: {out_lmp}" if not dry_run else "DRY RUN COMPLETE (nothing written)")


def cli(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)

    dry_run     = '--dry-run'     in args; args = [a for a in args if a != '--dry-run']
    force_split = '--force-split' in args; args = [a for a in args if a != '--force-split']
    coarsen     = '--no-coarsen'  not in args; args = [a for a in args if a != '--no-coarsen']

    max_atoms   = BOSS_CAP
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

    box_from_pdb = '--box-from-pdb' in args
    if box_from_pdb:
        args.remove('--box-from-pdb')
    box_from = None
    if '--box-from' in args:
        i = args.index('--box-from')
        box_from = args[i + 1]
        del args[i:i + 2]

    if len(args) != 4:
        print("Usage: molpargen run input.psf input.pdb species_dir/ out.lmp [options]")
        print()
        print("Runs the WHOLE pipeline end to end: extracts unique species from the")
        print("PSF+PDB, ligpargen per species/fragment, stitches every split species,")
        print("then combines -- returns only once out.lmp exists (or on error).")
        print()
        print("Options:")
        print("  --max-atoms N          BOSS/ligpargen atom cap (default: 220)")
        print("  --max-orbitals N       BOSS orbital-count cap passed to the splitter")
        print("                         (default: 480; hard BOSS limit is 550)")
        print("  --cgen CM1A|CM1A-LBCC  charge model (default: CM1A)")
        print("  --dry-run              print stages without running them (stops after")
        print("                         stage 0 -- manifest.json doesn't exist yet to plan against)")
        print("  --force-split          pass --force to the splitter (skip unsaturated backbone warnings)")
        print("  --box-from FILE        copy box bounds from this LAMMPS data file (passed to combine)")
        print("  --box-from-pdb         read box from input.pdb's CRYST1 record (passed to combine)")
        print("  --no-coarsen           skip the post-combine coarsening step (on by default);")
        print("                         collapses per-atom types to (mol × element), writing")
        print("                         <out>_coarse.data + <out>_topology.json")
        return 1

    psf_path, pdb_path, species_dir, out_lmp = args
    run(psf_path, pdb_path, species_dir, out_lmp,
        max_atoms=max_atoms, cgen=cgen,
        dry_run=dry_run, force_split=force_split,
        orbital_target=orbital_target,
        box_from=box_from, box_from_pdb=box_from_pdb,
        coarsen=coarsen)
    return 0


if __name__ == "__main__":
    sys.exit(cli())
